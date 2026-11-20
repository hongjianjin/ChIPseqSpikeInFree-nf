#!/usr/bin/env python3
"""Tests for bin/csif_bigwig_bin.py using small synthetic bigWigs (needs pybigtools, numpy)."""
import json, os, subprocess, sys, tempfile
import numpy as np
import pybigtools as pb

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "bin", "csif_bigwig_bin.py")
SIZES = {"chr1": 12345, "chr2": 5000, "chrM": 3000, "chrUn_x": 90000, "chr3": 399}  # chr2 is an exact multiple of 1000
rng = np.random.default_rng(1)


def make_truth(sizes):
    """sparse per-base signal: piecewise-constant intervals with gaps (absent from the bigWig)"""
    truth = {}
    for c, L in sizes.items():
        arr = np.zeros(L)
        pos = 0
        while pos < L:
            gap = int(rng.integers(0, 400)); w = int(rng.integers(1, 300))
            pos += gap
            if pos >= L: break
            end = min(L, pos + w)
            arr[pos:end] = float(rng.integers(1, 9)) * 0.37
            pos = end
        truth[c] = arr
    return truth


def write_bw(path, truth, rename=lambda x: x, scale=1.0, sizes=None):
    sizes = sizes or SIZES
    def it():
        for c, arr in truth.items():
            nz = np.flatnonzero(np.diff(np.concatenate(([0.0], arr, [0.0]))) != 0)
            for a, b in zip(nz[:-1], nz[1:]):
                if arr[a] != 0: yield (rename(c), int(a), int(b), float(arr[a]) * scale)
    chroms = {rename(c): L for c, L in sizes.items()}
    # relative path: pybigtools reads "C:\\..." as a URL when writing on Windows
    w = pb.open(os.path.relpath(path), "w"); w.write(chroms, it())


def expected(truth, depth, with_chr=True, skip=()):
    # chromosomes kept by GenerateBins: no random/Un/alt/hap, length >= 2000
    order = [c for c, L in SIZES.items() if "Un" not in c and L >= 2000]
    masses = []
    for c in order:
        L = SIZES[c]; arr = truth[c] if c not in skip else np.zeros(L)
        n_full = L // 1000
        m = [arr[i * 1000:(i + 1) * 1000].sum() for i in range(n_full)]
        if L % 1000: m.append(arr[n_full * 1000:].sum())
        else: m.append(0.0)  # GenerateBins trailing empty bin
        masses.append(np.array(m))
    mass = np.concatenate(masses)
    return mass / mass.sum() * depth


def run(tmp, bw, prefix, *extra, expect_fail=False):
    cs = os.path.join(tmp, "t.chrom.sizes")
    open(cs, "w").write("".join("%s\t%d\n" % kv for kv in SIZES.items()))
    r = subprocess.run([sys.executable, TOOL, "--bigwig", bw, "--id", "S1", "--genome", cs, "--prefix", os.path.join(tmp, prefix), *extra],
                       capture_output=True, text=True)
    if expect_fail:
        assert r.returncode != 0, "expected failure, got success: " + r.stdout
        return r.stderr
    assert r.returncode == 0, r.stderr
    d = np.genfromtxt(os.path.join(tmp, prefix + "_rawCounts.txt"), skip_header=1, dtype=str, delimiter="\t")
    return d[:, 0], d[:, 1].astype(float), json.load(open(os.path.join(tmp, prefix + "_bigwig_qc.json")))


def main():
    with tempfile.TemporaryDirectory(dir=".") as tmp:
        tmp = os.path.abspath(tmp)
        truth = make_truth(SIZES)
        p1 = os.path.join(tmp, "a.bw"); write_bw(p1, truth)

        # 1. values equal an independent numpy computation; totals equal the requested depth
        lab, cnt, qc = run(tmp, p1, "a", "--pseudo-depth", "5000000")
        exp = expected(truth, 5e6)
        assert len(cnt) == len(exp), (len(cnt), len(exp))
        assert np.allclose(cnt, exp, rtol=2e-5, atol=1e-4), np.abs(cnt - exp).max()
        assert abs(cnt.sum() - 5e6) / 5e6 < 1e-4
        assert lab[0] == "chr1:1-1000" and lab[12] == "chr1:12001-12345" and "chr2:5001-5000" in lab  # GenerateBins labels and trailing empty bin
        assert not any(l.startswith("chrUn") for l in lab)
        print("ok 1: bin sums match numpy; labels as GenerateBins")

        # 2. a per-sample constant (RPM / spike-in factor) cancels
        p2 = os.path.join(tmp, "b.bw"); write_bw(p2, truth, scale=3.7)
        lab2, cnt2, _ = run(tmp, p2, "b", "--pseudo-depth", "5000000")
        assert np.allclose(cnt, cnt2, rtol=1e-4, atol=1e-4)
        print("ok 2: constant rescaling of the bigWig gives identical counts")

        # 3. chromosome naming: 1/2/MT vs chr1/chr2/chrM; labels follow the bigWig, counts are equal
        ren = lambda c: {"chrM": "MT"}.get(c, c.replace("chr", ""))
        p3 = os.path.join(tmp, "c.bw"); write_bw(p3, truth, rename=ren)
        lab3, cnt3, _ = run(tmp, p3, "c", "--pseudo-depth", "5000000")
        assert np.allclose(cnt, cnt3, rtol=1e-4, atol=1e-4) and lab3[0] == "1:1-1000" and any(l.startswith("MT:") for l in lab3)
        print("ok 3: chr / no-chr and M / MT handled")

        # 4. --reads sets the library size
        _, cnt4, qc4 = run(tmp, p1, "d", "--reads", "31000000")
        assert abs(cnt4.sum() - 31e6) / 31e6 < 1e-4 and qc4["pseudo_depth_source"] == "reads"
        print("ok 4: --reads sets total")

        # 5. negative values are rejected (input-subtracted / log-ratio tracks)
        neg = {c: a.copy() for c, a in truth.items()}; neg["chr1"][100:400] = -1.0
        p5 = os.path.join(tmp, "e.bw"); write_bw(p5, neg)
        assert "negative" in run(tmp, p5, "e", expect_fail=True)
        print("ok 5: negative signal rejected")

        # 6. wrong assembly (different chromosome length) is rejected unless overridden
        sz2 = dict(SIZES); sz2["chr1"] = 12000
        t6 = {c: a[:sz2[c]] for c, a in truth.items()}
        p6 = os.path.join(tmp, "f.bw"); write_bw(p6, t6, sizes=sz2)
        assert "wrong assembly" in run(tmp, p6, "f", expect_fail=True)
        run(tmp, p6, "g", "--allow-length-mismatch")
        print("ok 6: length mismatch rejected, override works")

        # 7. chromosome absent from the bigWig counts as zero signal, with a warning in the QC file
        t7 = {c: a for c, a in truth.items() if c != "chrM"}; sz7 = {c: L for c, L in SIZES.items() if c != "chrM"}
        p7 = os.path.join(tmp, "h.bw"); write_bw(p7, t7, sizes=sz7)
        _, cnt7, qc7 = run(tmp, p7, "h", "--pseudo-depth", "5000000")
        assert qc7["chroms_missing_in_bigwig"] == ["chrM"]
        m7 = expected({**truth, "chrM": np.zeros(SIZES["chrM"])}, 5e6)
        assert np.allclose(cnt7, m7, rtol=2e-5, atol=1e-4)
        print("ok 7: missing chromosome = zero signal, reported")

        # 8. all-zero bigWig is rejected
        z = {c: np.zeros(L) for c, L in SIZES.items()}; z["chr1"][5:10] = 0.0
        def it():
            yield ("chr1", 0, 10, 0.0)
        pz = os.path.join(tmp, "z.bw"); w = pb.open(os.path.relpath(pz), "w"); w.write(SIZES, it())
        assert "no signal" in run(tmp, pz, "z", expect_fail=True)
        print("ok 8: empty track rejected")
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
