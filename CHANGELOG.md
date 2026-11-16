# ChIPseqSpikeInFree-nf: Changelog

## v0.1.0dev - [unreleased]

Initial pipeline structure following nf-core conventions:

- `CSIF_COUNT`: read counting in genome-wide bins and count parsing for the whole cohort
  (standard mode of `run_ChIPseqSpikeInFree.R`)
- `CSIF_SCALE`: scaling factors, QC calls, distribution and box plots from the parsed matrix
  (replotting mode), so calibration options can be changed without recounting
- nf-schema samplesheet (`sample,bam,antibody,group,color`) and parameter validation
- nf-core institutional profiles, including `-profile stjude`
- Container `ghcr.io/hongjianjin/chipseqspikeinfree` built by GitHub Actions
- Apache-2.0 license
