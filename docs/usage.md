# ChIPseqSpikeInFree-bigwig-nf: Usage

Only what differs from the BAM pipeline (ChIPseqSpikeInFree-nf) is described here. Running on the
St. Jude HPC (Nextflow >= 24.10 in a conda environment, `-profile stjude`, shared Singularity cache)
is the same; see that pipeline's `docs/usage.md`.

## Samplesheet

`--input` takes a CSV (or TSV) with a header row:

| Column     | Required | Description                                                                    |
| ---------- | -------- | ------------------------------------------------------------------------------ |
| `sample`   | yes      | Unique sample name (letters, digits, `.`, `_`, `-`). Used as the sample ID.    |
| `bigwig`   | yes      | bigWig coverage track: local path, or `http(s)://` / `ftp://` URL. Ends in `.bw` or `.bigWig`. |
| `antibody` | yes      | ChIP target, e.g. `H3K27me3`. Samples are compared within each antibody.      |
| `group`    | yes      | Condition, e.g. `WT`, `K27M`.                                                  |
| `color`    | no       | R colour name or hex code for plots; used only when every row has one.        |
| `reads`    | no       | Mapped reads in the sample (whole number). Sets the pseudo read total of that sample. |

* A URL is **not downloaded**: the bigWig is read with HTTP range requests. Servers must send
  `Accept-Ranges: bytes` (GEO's `ftp.ncbi.nlm.nih.gov` does; `ftp://` links are rewritten to `https://`).
  Whole-genome 1 kb bins still need nearly all base-level data, so expect to transfer about as much as
  the file size (in the benchmark: 90-300 MB per 120-470 MB file, 20-75 s). Coarse bins with
  `--approx_summary` can use zoom levels and transfer less, at some cost in exactness.
* A local file is staged into the task's work directory.
* All rows are analysed together in one run. Scaling factors from different runs are not comparable,
  so analyse samples together that should be compared, and separate systems that should not be
  (the paper ran each cell system and mark separately).

## What the bigWig must contain

A non-negative coverage track from duplicate-removed, uniquely mapped reads. A constant per-sample
factor (RPM, RPGC, spike-in scaling) cancels and is allowed. Input-subtracted, ratio and log tracks are
not valid; negative values stop the run, other transformations cannot be detected. Check the
repository's description of how the tracks were made.

The pipeline checks: no negative values; bigWig chromosome lengths equal the genome's (wrong assembly
stops the run); chromosomes absent from the bigWig count as zero signal and are listed in the QC file.

## How the counts are made

`CSIF_BIN_BIGWIG` builds the bins of `GenerateBins()` (1 kb by default, verified identical for hg19,
hg38, mm9 and mm10) and takes the **sum of the per-base bigWig values** in each bin; bases missing from a
sparse bigWig count as 0. The masses are rescaled so that they add up to the sample's `reads`
(or `--pseudo_depth`, default 20,000,000) and written as fractional pseudo-counts.

Because `ParseReadCounts()` divides every sample by its own total, the rescaling cancels in the parsed
matrix (checked: random per-sample factors changed it by at most 1.1e-15). The read total therefore
only matters for the method's 1e6-read cut-off (a sample whose `reads` is below 1e6 is dropped, as it
would be for a BAM). Use `reads` when you know the library sizes.

**Scaling factors can still differ between runs.** The last step (`CalculateSF()`) is numerically
unstable with the default `--min_first_turn auto`: in the benchmark, noise of 1e-14 relative size changed
the factor of 28 of 33 samples (up to x1.78), and re-reading the parsed matrix from a file (as
`CSIF_SCALE` does) changed 13 of 33. Setting `--min_first_turn` to a fixed value (0.5, the value the
package comments call empirically good) removed this in the benchmark (0 of 33 changed). This applies to
BAM input too. See `benchmark/README.md`.

## Parameters

| Parameter          | Default | Notes                                                                 |
| ------------------ | ------- | --------------------------------------------------------------------- |
| `--genome`         | `hg38`  | `hg19`, `hg38`, `mm9`, `mm10` (chrom sizes bundled in `assets/chrom_sizes`) |
| `--chrom_sizes`    | none    | Custom chrom.sizes file; overrides `--genome`                         |
| `--bin_size`       | `1000`  | 200-10000 bp                                                          |
| `--pseudo_depth`   | `2e7`   | Pseudo read total when a row has no `reads`; must be >= 2e6           |
| `--approx_summary` | `false` | Allow zoom-level summaries (faster remote reads, approximate)         |
| `--min_first_turn`, `--max_last_turn`, `--cutoff_qc`, `--xmax` | as the BAM pipeline | Calibration; changing them re-runs only `CSIF_SCALE` |

## Container

`--chipseqspikeinfree_container` must contain R, ChIPseqSpikeInFree, Python 3 and `pybigtools`
(`containers/chipseqspikeinfree-bigwig/Dockerfile`). It has not been built yet.

## Without Nextflow

`bin/csif_bigwig_bin.py` is a stand-alone script (Python >= 3.9, `numpy`, `pybigtools`, and `requests`
for URLs). Per sample:

```bash
python3 bin/csif_bigwig_bin.py --bigwig <path-or-url> --id <sample>.bw --genome mm9 --prefix <sample> [--reads N]
```

writes `<sample>_rawCounts.txt` and `<sample>_bigwig_qc.json`. This is how the benchmark in
`benchmark/` was run, followed by the package's own `ParseReadCounts()`, `CalculateSF()`,
`PlotDistr()` and `BoxplotSF()`.
