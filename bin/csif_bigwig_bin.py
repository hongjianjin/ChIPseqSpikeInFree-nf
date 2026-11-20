#!/usr/bin/env python3
"""
csif_bigwig_bin.py: per-bin signal from a bigWig, in the format ChIPseqSpikeInFree reads.

Replaces the BAM counting step (CountRawReads) of ChIPseqSpikeInFree. For one bigWig it

  1. builds the same genome-wide bins as ChIPseqSpikeInFree::GenerateBins() (overlap 0),
  2. sums the bigWig signal over every bin (signal mass = sum of per-base values over the bin;
     bases absent from a sparse bigWig count as 0),
  3. rescales the masses to pseudo read counts that add up to a chosen library size, and
  4. writes <prefix>_rawCounts.txt (columns: bin, <ID>) plus <prefix>_bigwig_qc.json.

Why rescaling is harmless: ParseReadCounts() divides every sample by its own total before
anything else (CPM = bin / sum(bins) * 1e6 * 1000 / binSize), so a per-sample constant in the
bigWig (RPM, RPGC, spike-in factor) cancels exactly. The library size only has to exceed the 1e6
cut-off of ParseReadCounts().

The bigWig must be a non-negative coverage track (read pileup, any constant scaling). Input-
subtracted, log-ratio or otherwise transformed tracks are not valid input and are rejected when
they contain negative values.

http(s) URLs are read with HTTP range requests (no full download).
"""
import argparse, collections, io, json, os, re, sys, time, warnings

import numpy as np

warnings.filterwarnings("ignore")


# --------------------------------------------------------------------------------------
# bins: faithful port of ChIPseqSpikeInFree::GenerateBins(overlap = 0)
# --------------------------------------------------------------------------------------
def read_chrom_sizes(spec, search_dirs=()):
    """spec: path to a chrom.sizes file, or a genome name (hg19, hg38, mm9, mm10, ...)."""
    cands = [spec] + [os.path.join(d, spec + ".chrom.sizes") for d in search_dirs]
    for c in cands:
        if c and os.path.isfile(c):
            rows = []
            with open(c) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t") if "\t" in line else line.split()
                    if len(f) >= 2 and f[0]:
                        rows.append((f[0], int(f[1])))
            return rows, c
    raise FileNotFoundError("chrom sizes not found for %r (looked in %s)" % (spec, cands))


def generate_bins(chrom_rows, bin_size=1000, with_chr=True, min_chr_size=None):
    """Return list of (label_chrom, source_chrom, start1, end1, chrom_len) in GenerateBins order.

    Reproduces GenerateBins() quirks: contigs matching random/Un/alt/hap are dropped, contigs
    shorter than 2 * bin_size are dropped, chrM is labelled MT (chrM with chr), and a chromosome
    whose length is an exact multiple of bin_size gets a trailing empty bin [L+1, L].
    """
    if bin_size < 200 or bin_size > 10000:
        raise ValueError("Recommended binSize range 200 ~ 10000 (bp); got %s" % bin_size)
    min_chr_size = bin_size * 2 if min_chr_size is None else min_chr_size
    out = []
    for name, length in chrom_rows:
        if any(s in name for s in ("random", "Un", "alt", "hap")):
            continue
        lab = name.replace("chrM", "MT").replace("chr", "")
        if with_chr:
            lab = ("chr" + lab).replace("chrMT", "chrM")
        if length < min_chr_size:
            continue
        starts = np.arange(0, length + 1, bin_size) + 1
        ends = np.concatenate([np.arange(bin_size, length + 1, bin_size), [length]])
        starts = starts[: len(ends)]
        out.append((lab, name, starts.astype(np.int64), ends.astype(np.int64), length))
    return out


# --------------------------------------------------------------------------------------
# input: local path or HTTP(S) URL through range requests
# --------------------------------------------------------------------------------------
class RangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP(S) using Range requests and a small block cache.

    Uses `requests`, so HTTP(S)_PROXY settings are honoured; pybigtools' own URL support does not
    go through a proxy.
    """

    def __init__(self, url, block=1 << 20, max_blocks=64):
        import requests

        self._requests = requests
        self.block, self.max_blocks = block, max_blocks
        self.s = requests.Session()
        r = self.s.head(url, allow_redirects=True, timeout=60)
        r.raise_for_status()
        self.url = r.url
        self.size = int(r.headers["Content-Length"])
        if r.headers.get("Accept-Ranges", "").lower() != "bytes":
            raise IOError("server does not advertise byte-range support: %s" % url)
        self.pos = 0
        self.cache = collections.OrderedDict()
        self.n_requests = 0
        self.bytes_fetched = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def _get(self, b):
        if b in self.cache:
            self.cache.move_to_end(b)
            return self.cache[b]
        lo = b * self.block
        hi = min(self.size, lo + self.block) - 1
        for attempt in range(6):
            try:
                r = self.s.get(self.url, headers={"Range": "bytes=%d-%d" % (lo, hi)}, timeout=120)
                if r.status_code == 206:
                    break
                raise IOError("unexpected HTTP %s for a range request" % r.status_code)
            except (self._requests.RequestException, IOError):
                if attempt == 5:
                    raise
                time.sleep(1 + 2 * attempt)
        self.n_requests += 1
        self.bytes_fetched += len(r.content)
        self.cache[b] = r.content
        if len(self.cache) > self.max_blocks:
            self.cache.popitem(last=False)
        return r.content

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        out, end = [], self.pos + n
        while self.pos < end:
            b, o = divmod(self.pos, self.block)
            chunk = self._get(b)[o : o + (end - self.pos)]
            out.append(chunk)
            self.pos += len(chunk)
        return b"".join(out)

    def readinto(self, buf):
        d = self.read(len(buf))
        buf[: len(d)] = d
        return len(d)


def open_bigwig(src):
    import pybigtools as pb

    if re.match(r"^(ftp|http|https)://", src):
        url = re.sub(r"^ftp://", "https://", src)  # GEO lists ftp:// links that are also served over https
        handle = RangeFile(url)
        return pb.open(handle), handle
    if os.name == "nt" and re.match(r"^[A-Za-z]:[\\/]", src):
        # pybigtools reads "C:\..." as a URL scheme; hand it an open file instead
        return pb.open(open(src, "rb")), None
    return pb.open(src), None


def norm_chrom(x):
    x = re.sub(r"^chr", "", x)
    return "MT" if x in ("M", "MT") else x


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bigwig", required=True, help="bigWig path or http(s) URL")
    ap.add_argument("--id", required=True, help="column name written for this sample (the ID used in the metadata file)")
    ap.add_argument("--genome", required=True, help="hg19 | hg38 | mm9 | mm10 or a chrom.sizes file")
    ap.add_argument("--prefix", required=True, help="output prefix (<prefix>_rawCounts.txt, <prefix>_bigwig_qc.json)")
    ap.add_argument("--bin-size", type=int, default=1000)
    ap.add_argument("--reads", type=float, default=None, help="library size (mapped reads) of this sample; sets the pseudo-count total")
    ap.add_argument("--pseudo-depth", type=float, default=20e6, help="pseudo-count total when --reads is not given [%(default).0f]")
    ap.add_argument("--chrom-sizes-dir", action="append", default=[], help="folder with <genome>.chrom.sizes files")
    ap.add_argument("--approx", action="store_true", help="allow zoom-level summaries (faster for remote files, approximate)")
    ap.add_argument("--allow-length-mismatch", action="store_true", help="do not stop when bigWig chromosome lengths differ from the genome")
    args = ap.parse_args(argv)
    t0 = time.time()

    here = os.path.dirname(os.path.abspath(__file__))
    dirs = args.chrom_sizes_dir + [os.path.join(here, "..", "assets", "chrom_sizes"), os.path.join(here, "assets", "chrom_sizes")]
    chrom_rows, chrom_file = read_chrom_sizes(args.genome, dirs)

    bw, handle = open_bigwig(args.bigwig)
    bw_chroms = bw.chroms()
    if not bw_chroms:
        sys.exit("ERROR: bigWig has no chromosomes: %s" % args.bigwig)
    with_chr = any(c.startswith("chr") for c in bw_chroms)

    # a coverage track has no negative values; reject input-subtracted / log-ratio tracks up front
    summ = (bw.info() or {}).get("summary") or {}
    if summ.get("min") is not None and summ["min"] < 0:
        sys.exit("ERROR: bigWig contains negative values (min %.4g). ChIPseqSpikeInFree needs a non-negative "
                 "coverage track (read pileup); input-subtracted or log-ratio tracks are not valid input." % summ["min"])

    # map genome chromosomes to bigWig chromosome names (chr1/1, chrM/MT/M)
    bw_map = {}
    for c in bw_chroms:
        k = norm_chrom(c)
        if k in bw_map:
            sys.exit("ERROR: bigWig chromosome names are ambiguous after chr/M normalisation: %s, %s" % (bw_map[k], c))
        bw_map[k] = c

    bins = generate_bins(chrom_rows, args.bin_size, with_chr)
    labels, masses = [], []
    used, missing, mismatched = [], [], []
    n_neg = 0
    for lab, src_name, starts, ends, length in bins:
        bwname = bw_map.get(norm_chrom(src_name))
        n = len(starts)
        vals = np.zeros(n, dtype=np.float64)
        if bwname is None:
            missing.append(src_name)
        else:
            if bw_chroms[bwname] != length:
                mismatched.append((src_name, length, bw_chroms[bwname]))
                if not args.allow_length_mismatch:
                    continue
            used.append(src_name)
            clen = min(length, bw_chroms[bwname])
            nfull = clen // args.bin_size
            if nfull > 0:
                v = bw.values(bwname, 0, nfull * args.bin_size, bins=nfull, summary="sum",
                              exact=not args.approx, fillna=0)
                v = np.nan_to_num(np.asarray(v, dtype=np.float64), nan=0.0)
                n_neg += int((v < 0).sum())
                vals[:nfull] = v
            if clen > nfull * args.bin_size and n > nfull:  # partial last bin
                v = bw.values(bwname, nfull * args.bin_size, clen, bins=1, summary="sum",
                              exact=not args.approx, fillna=0)
                vals[nfull] = float(np.nan_to_num(np.asarray(v, dtype=np.float64), nan=0.0)[0])
        labels.append(np.char.add(np.char.add(np.char.add(lab + ":", starts.astype(str)), "-"), ends.astype(str)))
        masses.append(vals)

    if mismatched and not args.allow_length_mismatch:
        sys.exit("ERROR: bigWig chromosome lengths differ from genome '%s' (wrong assembly?): %s. "
                 "Use --allow-length-mismatch to override." % (args.genome, ", ".join("%s genome=%d bigWig=%d" % m for m in mismatched[:5])))
    if not used:
        sys.exit("ERROR: no bigWig chromosome matches genome '%s' (bigWig has: %s)" % (args.genome, ", ".join(list(bw_chroms)[:8])))
    if n_neg:
        sys.exit("ERROR: %d bins have negative signal; not a coverage track." % n_neg)

    labels = np.concatenate(labels)
    mass = np.concatenate(masses)
    total = float(mass.sum())
    if not total > 0:
        sys.exit("ERROR: no signal in any bin")
    depth = float(args.reads) if args.reads else float(args.pseudo_depth)
    counts = mass / total * depth

    out = args.prefix + "_rawCounts.txt"
    with open(out, "w", newline="\n") as fh:
        fh.write("bin\t%s\n" % args.id)
        chunk = 500000
        for i in range(0, len(labels), chunk):
            fh.write("".join("%s\t%.6g\n" % (l, c) for l, c in zip(labels[i : i + chunk].tolist(), counts[i : i + chunk].tolist())))

    cpm = counts / counts.sum() * 1e6 * (1000.0 / args.bin_size)
    qc = {
        "id": args.id, "bigwig": args.bigwig, "genome": args.genome, "chrom_sizes": os.path.basename(chrom_file),
        "bin_size": args.bin_size, "n_bins": int(len(labels)), "n_bins_with_signal": int((mass > 0).sum()),
        "fraction_bins_with_signal": float((mass > 0).mean()),
        "total_signal_mass": total, "mean_signal_per_bp": total / (args.bin_size * len(labels)),
        "pseudo_depth": depth, "pseudo_depth_source": "reads" if args.reads else "pseudo_depth_default",
        "fraction_bins_cpm_gt_150": float((cpm > 150).mean()),
        "chroms_used": len(used), "chroms_missing_in_bigwig": missing, "chroms_length_mismatch": mismatched,
        "summary_statistics": "zoom-level allowed (approx)" if args.approx else "exact",
        "remote": handle is not None,
        "http_requests": handle.n_requests if handle else None,
        "http_bytes_fetched": handle.bytes_fetched if handle else None,
        "seconds": round(time.time() - t0, 1),
    }
    with open(args.prefix + "_bigwig_qc.json", "w") as fh:
        json.dump(qc, fh, indent=2)
    if missing:
        print("WARNING: genome chromosomes absent from the bigWig (treated as zero signal): %s" % ", ".join(missing), file=sys.stderr)
    print("%s: %d bins, %.1f%% with signal, %s%s in %.0fs" % (
        args.id, len(labels), 100 * qc["fraction_bins_with_signal"],
        ("%.1f MB fetched over %d requests, " % (handle.bytes_fetched / 1e6, handle.n_requests)) if handle else "",
        "approx" if args.approx else "exact", time.time() - t0))


if __name__ == "__main__":
    main()
