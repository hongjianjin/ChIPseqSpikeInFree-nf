# ChIPseqSpikeInFree-bigwig-nf

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)

**Status: development (v0.1.0dev).** The Python binning tool and the R parse/scale steps are tested
(see [benchmark/](benchmark/)). The Nextflow wiring has **not been run** (it was written on a Windows
host without Java, Nextflow or Docker) and the container image has **not been built**. Run the stub
test and build the image before trusting a real run.

## Introduction

A variant of [ChIPseqSpikeInFree-nf](https://github.com/hongjianjin/ChIPseqSpikeInFree-nf) that takes
**bigWig coverage tracks instead of BAM files**. It is for data where only bigWigs are available, such
as the processed tracks of public repositories (GEO, ChIP-Atlas), and it can read bigWigs over HTTP
without downloading them.

ChIPseqSpikeInFree estimates scaling factors from how reads are spread over genome-wide bins: it
turns each sample into the cumulative proportion of reads below a series of bin-CPM cutoffs and
reads the scaling factor from turning points of that curve. `ParseReadCounts()` divides every sample
by its own bin total first, so only the *relative* signal of the bins matters, and a bigWig can
supply it.

1. `CSIF_BIN_BIGWIG`: one job per bigWig. Builds the same bins as `GenerateBins()`, sums the bigWig
   signal over each bin, rescales to pseudo read counts and writes `<sample>_rawCounts.txt`
   (`bin/csif_bigwig_bin.py`). Rejects negative-valued tracks and wrong assemblies.
2. `CSIF_MERGE_COUNTS`: merges the per-sample tables into one bin x sample table (numeric columns).
3. `CSIF_PARSE`: unchanged from the BAM pipeline (`run_ChIPseqSpikeInFree.R`, standard mode).
4. `CSIF_SCALE`: unchanged (replotting mode); changing calibration options re-runs only this step.

## What a bigWig must be

* a **non-negative coverage track** of the ChIP sample: read pileup from duplicate-removed, uniquely
  mapped reads (what the BAM pipeline asks for);
* any **constant per-sample scaling** is fine (RPM, RPGC, a spike-in factor). It cancels in the
  CPM step: in a test, scaling the counts by random per-sample factors changed the parsed matrix by
  at most 1.1e-15 (the final scaling factor can still move, see "Known limitation");
* **not** input-subtracted, ratio or log-transformed tracks, and not smoothed or peak-only tracks.
  Negative values are rejected; other transformations cannot be detected, so check how the track was made;
* on the genome given with `--genome` (chromosome lengths are checked; a mismatch stops the run).

[benchmark/](benchmark/) compares the bigWig route with the published BAM route and spike-in factors
for 33 public samples. It does **not** separate the effects of read extension, smoothing or
filtering in the bigWigs, and the BAM-route values are the published ones, not re-run.

## Known limitation: numerically unstable scaling factors

In the benchmark, `CalculateSF()` with the default `--min_first_turn auto` gave a different scaling
factor for 28 of 33 samples after a relative perturbation of 1e-14 of the parsed matrix (largest
change x1.78). The cause is the automatic first-turning-point step, whose second density is nearly
flat, so its maximum is decided by round-off. The same instability applies to BAM input. Consequences:

* the same data can give different factors in one-process runs and in this pipeline, whose
  `CSIF_SCALE` step re-reads the parsed matrix from a file (13 of 33 differed in the benchmark);
* a single sample's factor is only reliable to about +/-0.4 on the log2 scale (range seen: x1.0-1.78);
* fixing the first turning point removed the instability in the benchmark (`--min_first_turn 0.5`:
  0 of 33 samples changed). This is an exploratory finding on one dataset; it changes the estimator
  relative to the published values.

Group-level fold changes (mean factor of K27M over WT) were not clearly more stable: across 50 noise
replicates they varied by x1.13-1.37 (largest over smallest) per group, versus a median x1.0-1.5 for single
samples (`benchmark/dataset1_bigwig_group_fold_change_under_noise.csv`).

## Usage

Samplesheet, one row per bigWig (local path or http(s)/ftp URL):

```csv
sample,bigwig,antibody,group,color,reads
K27M_rep1,/path/to/bigwigs/H3K27me3-K27M.rep1.bw,H3K27me3,K27M,green,25000000
K27M_rep2,/path/to/bigwigs/H3K27me3-K27M.rep2.bw,H3K27me3,K27M,green,31000000
WT_rep1,https://example.org/tracks/H3K27me3-WT.rep1.bw,H3K27me3,WT,grey,28000000
WT_rep2,https://example.org/tracks/H3K27me3-WT.rep2.bw,H3K27me3,WT,grey,
```

`color` and `reads` are optional. `reads` is the sample's number of mapped reads; leave it empty to use
`--pseudo_depth`. All rows are analysed together in one run, so put every sample that has to be
comparable in the same samplesheet, and run each cell system / antibody combination you do not want
compared with each other as separate runs (scaling factors from different runs are not comparable).

```bash
nextflow run . -profile stjude --input samplesheet.csv --outdir results --cohort H3K27me3 --genome hg38
```

Wiring test (no data needed): `nextflow run . -profile test -stub --outdir results_test`.

### Before the first run

1. Make the tool executable: `chmod +x bin/csif_bigwig_bin.py` (or `git update-index --chmod=+x bin/csif_bigwig_bin.py`
   before pushing). Files written on Windows have no executable bit and Nextflow runs the script from `bin/`.
2. Build and push the container (`containers/chipseqspikeinfree-bigwig/Dockerfile`), or point
   `--chipseqspikeinfree_container` at your own image.
3. Run the stub test above.

### Without Nextflow

```bash
python3 bin/csif_bigwig_bin.py --bigwig sample.bw --id sample.bw --genome mm9 --prefix sample   # per sample
```

This step is tested. The remaining steps (merge the per-sample tables on the `bin` column into
`<prefix>_rawCounts.txt`, write the metadata file, then run
`Rscript run_ChIPseqSpikeInFree.R -L bam.lst -G mm9 -m sample_meta.txt -o <prefix>`, which loads the existing
rawCounts file) were **not run as a command here**: the benchmark called the package's
`ParseReadCounts()`, `CalculateSF()`, `PlotDistr()` and `BoxplotSF()` directly from R.

See [docs/usage.md](docs/usage.md) and [docs/output.md](docs/output.md).

## Container

`--chipseqspikeinfree_container` (default `ghcr.io/hongjianjin/chipseqspikeinfree-bigwig:0.1.0`,
**not published**). Build it from [containers/chipseqspikeinfree-bigwig/Dockerfile](containers/chipseqspikeinfree-bigwig/Dockerfile),
which extends the published ChIPseqSpikeInFree 1.2.6 image with Python and pybigtools.

## Tests

* `python3 tests/test_bigwig_bin.py`: 8 checks of the binning tool on synthetic bigWigs (needs
  `pybigtools`, `numpy`; passes with pybigtools 0.3.0, numpy 2.4.6).
* The bins were compared with the package's `GenerateBins()` for hg19, hg38, mm9 and mm10, with and
  without `chr`: identical (see benchmark/).

## Credits and license

Built on ChIPseqSpikeInFree-nf and the ChIPseqSpikeInFree R package by Hongjian Jin (St. Jude
Children's Research Hospital). Pipeline code: Apache License 2.0 ([LICENSE](LICENSE)); the R package
is GPL-3. The chrom.sizes files in `assets/chrom_sizes/` are copied from that package.

## Citations

If you use this pipeline, please cite:

> Jin H, Kasper LH, Larson JD, et al. ChIPseqSpikeInFree: a ChIP-seq normalization approach to
> reveal global changes in histone modifications without spike-in. *Bioinformatics*. 2020;36(4):1270-1272.
> [doi:10.1093/bioinformatics/btz720](https://doi.org/10.1093/bioinformatics/btz720)

and [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/). The pipeline structure follows the
[nf-core](https://nf-co.re) framework ([Ewels et al. 2020](https://doi.org/10.1038/s41587-020-0439-x)).
