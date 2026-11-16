# ChIPseqSpikeInFree-nf: Output

```
<outdir>/
├── counts/                                  # CSIF_COUNT
│   ├── <cohort>_sample_meta.txt             # metadata generated from the samplesheet (ID = <sample>.bam)
│   ├── <cohort>_rawCounts.txt               # reads per genome-wide bin per sample
│   ├── <cohort>_parsedMatrix.txt            # cumulative read distributions used for calibration
│   └── <cohort>_count.log
├── scaling_factors/                         # CSIF_SCALE
│   ├── <cohort>_SF.txt                      # scaling factors and QC calls per sample
│   ├── <cohort>_distribution.pdf            # cumulative distribution curves with turning points
│   ├── <cohort>_boxplot.pdf                 # scaling factors by group
│   └── <cohort>_scale.log
└── pipeline_info/
    ├── chipseqspikeinfree_software_versions.yml
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.html
```

## `<cohort>_SF.txt`

| Column     | Meaning                                                                          |
| ---------- | -------------------------------------------------------------------------------- |
| `ID`       | `<sample>.bam`                                                                   |
| `GROUP`, `ANTIBODY`, `COLOR` | from the samplesheet                                           |
| `QC`       | `pass`, or the reason a sample failed (e.g. complete loss, input or poor enrichment) |
| `TURNS`    | turning points of the cumulative distribution used for the calculation          |
| `SF`       | scaling factor (`NA` for failed samples)                                          |

To apply the scaling factors, multiply a sample's normalized signal (e.g. reads per million in
a bigWig) by its `SF`; see the ChIPseqSpikeInFree README for details.
