# ChIPseqSpikeInFree-bigwig-nf: Changelog

## v0.1.0dev (unreleased)

First version, derived from ChIPseqSpikeInFree-nf v1.2.8. Implements component 1 of the ChIP-Atlas
global-shift audit plan: the ChIPseqSpikeInFree estimator from bigWig coverage.

### Pipeline

- `CSIF_BIN_BIGWIG` replaces `CSIF_COUNT_BAM`: one job per bigWig (local file or URL read with HTTP
  range requests). `bin/csif_bigwig_bin.py` builds the bins of `GenerateBins()` (checked identical for
  hg19, hg38, mm9, mm10, with and without `chr`), sums the bigWig signal per bin, rescales to pseudo
  read counts and writes `<sample>_rawCounts.txt` plus a QC JSON. It rejects negative-valued tracks and
  chromosome-length mismatches.
- `CSIF_MERGE_COUNTS` reads numeric instead of integer columns.
- `CSIF_PARSE` and `CSIF_SCALE` are unchanged.
- Samplesheet: `bigwig` replaces `bam`; optional `reads` column. Params: `--pseudo_depth`,
  `--approx_summary` replace `--paired_end`.
- Container: `containers/chipseqspikeinfree-bigwig/Dockerfile` (extends the ChIPseqSpikeInFree image
  with Python and pybigtools).

### Not done / not verified

- The Nextflow workflow has not been run (no Java, Nextflow or Docker on the development host); the
  `-stub` test and a real run are still to do, as is building the container.
- bigWig tracks other than the GEO tracks of GSE108344 have not been tested (see `benchmark/`).
