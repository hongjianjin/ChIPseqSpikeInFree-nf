# ChIPseqSpikeInFree-nf

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**ChIPseqSpikeInFree-nf** is a Nextflow pipeline, laid out following nf-core conventions, that
computes spike-in free scaling factors for ChIP-seq data with
[ChIPseqSpikeInFree](https://github.com/stjude/ChIPseqSpikeInFree). The scaling factors reveal
global changes in histone modification occupancy (e.g. H3K27me3 loss in H3K27M cells) without
exogenous spike-in chromatin or peak calling.

1. `CSIF_COUNT_BAM`: counts reads in genome-wide bins with ChIPseqSpikeInFree's
   `CountRawReads()`, one job per BAM, so all BAMs are counted in parallel.
2. `CSIF_MERGE_COUNTS`: merges the per-BAM counts into one bin x sample table.
3. `CSIF_PARSE`: parses the merged counts into cumulative distributions
   (`run_ChIPseqSpikeInFree.R`, standard mode, which reuses the merged counts).
4. `CSIF_SCALE`: computes scaling factors and QC calls and draws the distribution and box plots
   from the parsed matrix (replotting mode). Changing the calibration options re-runs only this
   step, not the read counting.

## Usage

Prepare a samplesheet with one row per BAM:

```csv
sample,bam,antibody,group,color
K27M_rep1,/path/to/bams/H3K27me3-K27M.rep1.rmdup.bam,H3K27me3,K27M,green
K27M_rep2,/path/to/bams/H3K27me3-K27M.rep2.rmdup.bam,H3K27me3,K27M,green
WT_rep1,/path/to/bams/H3K27me3-WT.rep1.rmdup.bam,H3K27me3,WT,grey
WT_rep2,/path/to/bams/H3K27me3-WT.rep2.rmdup.bam,H3K27me3,WT,grey
```

`color` is optional (used in plots when every row has one). Use duplicate-removed BAMs with
low-quality reads filtered.

Run on the St. Jude HPC with the [nf-core `stjude` institutional profile](https://nf-co.re/configs/stjude/)
(LSF executor + Singularity):

```bash
nextflow run hongjianjin/ChIPseqSpikeInFree-nf \
    -r v1.2.6 \
    -profile stjude \
    --input samplesheet.csv \
    --outdir results \
    --cohort H3K27me3 \
    --genome hg38
```

Nextflow 24.10 or newer (with Java 17+) is required. If the cluster's `nextflow` modules are
older, use a conda environment, which bundles its own Java:

```bash
conda create -n nextflow -c conda-forge -c bioconda nextflow=24.10.5
conda activate nextflow
```

`-r v1.2.6` pins this release; use `-r main -latest` for the development version.

`-profile stjude` is pulled from [nf-core/configs](https://github.com/nf-core/configs) at runtime.
On nodes without internet access, clone nf-core/configs and pass `--custom_config_base /path/to/configs`.

Other profiles: `docker`, `singularity`, `apptainer`, `test`, `debug`.

See [docs/usage.md](docs/usage.md) for all options and [docs/output.md](docs/output.md) for the
output layout. `nextflow run hongjianjin/ChIPseqSpikeInFree-nf --help` lists every parameter.

## Container

Both processes run in one container with R and ChIPseqSpikeInFree, set with
`--chipseqspikeinfree_container` (default `ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6`).
The image is built from [containers/chipseqspikeinfree/Dockerfile](containers/chipseqspikeinfree/Dockerfile)
by the [Build ChIPseqSpikeInFree container](.github/workflows/build-container.yml) GitHub Actions
workflow, on pushes that change the Dockerfile or manually for another package tag:

```bash
gh workflow run build-container.yml -R hongjianjin/ChIPseqSpikeInFree-nf -f csif_ref=1.2.6
```

You can also point the parameter at a local `.sif` file.

## Testing

A wiring test on placeholder BAMs (downloaded from this repo, so it runs from any folder):

```bash
nextflow run hongjianjin/ChIPseqSpikeInFree-nf -r main -profile test -stub --outdir results_test
```

`-profile test` uses no container; add `stjude` (`-profile stjude,test`) to also check the
Singularity image and LSF submission.

## Credits

ChIPseqSpikeInFree-nf was written by Hongjian Jin (St. Jude Children's Research Hospital) on top
of the ChIPseqSpikeInFree R package.

## License

ChIPseqSpikeInFree-nf is released under the [Apache License 2.0](LICENSE). The ChIPseqSpikeInFree
R package it runs is GPL-3.

## Citations

If you use ChIPseqSpikeInFree-nf, please cite:

> Jin H, Kasper LH, Larson JD, et al. ChIPseqSpikeInFree: a ChIP-seq normalization approach to
> reveal global changes in histone modifications without spike-in. *Bioinformatics*. 2020;36(4):1270-1272.
> [doi:10.1093/bioinformatics/btz720](https://doi.org/10.1093/bioinformatics/btz720)

and [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/). The pipeline structure follows the
[nf-core](https://nf-co.re) framework ([Ewels et al. 2020](https://doi.org/10.1038/s41587-020-0439-x)).
