# ChIPseqSpikeInFree-nf: Usage

## Samplesheet

`--input` takes a CSV (or TSV) with a header row:

| Column     | Required | Description                                                                    |
| ---------- | -------- | ------------------------------------------------------------------------------ |
| `sample`   | yes      | Unique sample name (letters, digits, `.`, `_`, `-`). Used as the sample ID.    |
| `bam`      | yes      | Path to the BAM file (local, `s3://`, `https://`).                             |
| `antibody` | yes      | ChIP target, e.g. `H3K27me3`, or `INPUT` for input controls.                   |
| `group`    | yes      | Condition, e.g. `WT`, `K27M`.                                                  |
| `color`    | no       | R colour name or hex code for plots; used only when every row has one.        |

All rows are analysed together in one ChIPseqSpikeInFree run, and samples are compared within
each antibody. Each BAM is counted in its own job. Use
duplicate-removed BAMs with low-quality reads filtered; no BAM index is needed.

The pipeline writes the metadata file ChIPseqSpikeInFree expects (`ID`, `ANTIBODY`, `GROUP`,
`COLOR`) with `ID = <sample>.bam`; each BAM is linked under that name, so all outputs are labelled
by sample name. A copy is kept as `<outdir>/counts/<cohort>_sample_meta.txt`.

## Running on the St. Jude HPC

The HPC `nextflow/*` modules (up to 23.10) and default Java (1.8) are too old; the pipeline needs
Nextflow >= 24.10 with Java 17+. Create a conda environment once:

```bash
module load conda3/202402
conda create -y -n nextflow -c conda-forge -c bioconda nextflow=24.10.5
```

Then, in each session:

```bash
module load conda3/202402
conda activate nextflow
nextflow run hongjianjin/ChIPseqSpikeInFree-nf -r v1.2.7 -profile stjude --input samplesheet.csv --genome hg38 --outdir results
```

The `stjude` profile from nf-core/configs submits each process as an LSF job, picks the queue
from the requested resources, and runs tasks in Singularity. Launch Nextflow itself from a
long-running job, not on a login node, e.g.

```bash
bsub -q standard -n 1 -M 8000 -W 72:00 -J csif -oo csif.log \
    "source \$(conda info --base)/etc/profile.d/conda.sh; conda activate nextflow; \
     nextflow run hongjianjin/ChIPseqSpikeInFree-nf -r v1.2.7 -profile stjude --input samplesheet.csv --genome hg38 --outdir results -resume"
```

Set `NXF_SINGULARITY_CACHEDIR` to a shared folder so the container is converted only once.

## Parameters

### Read counting (`CSIF_COUNT_BAM` per BAM, `CSIF_MERGE_COUNTS`, `CSIF_PARSE`)

| Parameter       | Default | Script option  | Notes                                            |
| --------------- | ------- | -------------- | ------------------------------------------------ |
| `--genome`      | `hg38`  | `--genome`     | `hg19`, `hg38`, `mm9`, `mm10`                    |
| `--chrom_sizes` | none    | `--genome`     | Custom chrom.sizes file; overrides `--genome`    |
| `--bin_size`    | `1000`  | `--binSize`    | 200-10000 bp                                     |
| `--paired_end`  | `false` | `--paired`     | Count fragments of paired-end BAMs               |

Each BAM is counted in its own job (in parallel), the counts are merged, and `CSIF_PARSE`
runs the script's standard mode on the merged table (`--ncores` = the task's CPUs). Changing
any of these options re-runs the counting. BAMs with fewer than 1e6 counted reads are dropped
with a warning, as in ChIPseqSpikeInFree itself.

### Scaling factors (`CSIF_SCALE`, replotting mode)

| Parameter          | Default | Script option    | Notes                                         |
| ------------------ | ------- | ---------------- | --------------------------------------------- |
| `--min_first_turn` | auto    | `--minFirstTurn` | 0.2-0.8                                       |
| `--max_last_turn`  | `0.99`  | `--maxLastTurn`  | 0.90-0.99; try `0.97` for H3K9me3             |
| `--cutoff_qc`      | `1.2`   | `--cutoffQC`     | 0.90-1.20; smaller is looser; try `1` for H3K9me3 |
| `--xmax`           | auto    | `--xmax`         | x-axis maximum of the distribution plot       |

With `-resume`, changing only these options re-runs `CSIF_SCALE` from the parsed matrix; the
read counting is reused. Example for H3K9me3:

```bash
nextflow run hongjianjin/ChIPseqSpikeInFree-nf -r v1.2.7 -profile stjude --input samplesheet.csv \
    --genome hg38 --outdir results --max_last_turn 0.97 --cutoff_qc 1 -resume
```

## Resources

`CSIF_COUNT_BAM` requests 2 CPUs, 12 GB and 20 h per BAM (labels `process_low`, `process_long`);
`CSIF_MERGE_COUNTS` and `CSIF_PARSE` 2 CPUs, 32 GB and 8 h; `CSIF_SCALE` 1 CPU, 6 GB and 4 h.
Raise memory for small bins (e.g. `--bin_size 200`) or many samples. Requests double on a retry after an out-of-memory or time-limit
kill. Override them in a custom config passed with `-c`:

```groovy
process {
    withName: 'CSIF_COUNT_BAM' { memory = 24.GB }
    withName: 'CSIF_PARSE'     { cpus = 8; memory = 64.GB }
}
```

## Container

`--chipseqspikeinfree_container` accepts a registry image or a local Singularity image. See
`containers/chipseqspikeinfree/Dockerfile`.

## Troubleshooting

**`unable to load shared object '/home/<user>/R/.../*.so'` or `libRlapack.so: cannot open shared object file`**

R inside the container picked up a personal R library exported on the host (`R_LIBS`). The
pipeline and the image both override it. When running the image by hand, add `--cleanenv`:

```bash
singularity exec --cleanenv docker://ghcr.io/hongjianjin/chipseqspikeinfree:1.2.6 Rscript -e 'library(ChIPseqSpikeInFree)'
```

**`Plugin nf-schema@2.4.2 requires Nextflow version >=24.10.0`**

An old Nextflow is on the `PATH`; activate the conda environment above.

**`Process requirement exceeds available CPUs`**

`-profile test` runs tasks on the current node; real runs should use `-profile stjude` (LSF).
