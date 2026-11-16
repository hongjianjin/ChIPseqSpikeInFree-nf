# ChIPseqSpikeInFree-nf: Changelog

## v1.2.7dev - [unreleased]

### Pipeline

- Read counting runs in parallel, one job per BAM (`CSIF_COUNT_BAM`, using ChIPseqSpikeInFree's
  `CountRawReads()`), instead of all BAMs in one job. `CSIF_MERGE_COUNTS` merges the per-BAM
  counts and `CSIF_PARSE` (formerly part of `CSIF_COUNT`) parses the merged table with the
  script's standard mode, so results match the single-job run.
- BAMs with < 1e6 counted reads are dropped with a warning in `<cohort>_merge.log`, as before.
- The cohort-wide "BAM file names must be unique" restriction is gone; each BAM is counted on its own.
- Outputs: `counts/<cohort>_count.log` is replaced by `_merge.log` and `_parse.log`.

## v1.2.6 - 2026-11-16

First release. Versioned to match the ChIPseqSpikeInFree release it runs (1.2.6).

### Pipeline

- `CSIF_COUNT`: read counting in genome-wide bins and count parsing for the whole cohort
  (standard mode of `run_ChIPseqSpikeInFree.R`); each BAM is linked as `<sample>.bam` so
  results are labelled by sample name; `--genome` (hg19, hg38, mm9, mm10) or `--chrom_sizes`,
  `--bin_size`, `--paired_end`
- `CSIF_SCALE`: scaling factors, QC calls, distribution and box plots from the parsed matrix
  (replotting mode); `--min_first_turn`, `--max_last_turn`, `--cutoff_qc`, `--xmax` re-run only
  this step with `-resume`
- nf-schema samplesheet (`sample,bam,antibody,group,color`) and parameter validation
- nf-core institutional profiles, including `-profile stjude` (LSF + Singularity)
- `-profile test -stub` wiring test on placeholder BAMs, runnable from any folder

### Container

- `ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6`: Bioconductor 3.20, Rsamtools,
  GenomicAlignments, ChIPseqSpikeInFree 1.2.6; built by GitHub Actions
- Host `R_LIBS` / `R_LIBS_USER` are kept out of the container (image and pipeline `env`)

### License

- Pipeline code: Apache-2.0. The ChIPseqSpikeInFree R package is GPL-3.
