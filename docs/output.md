# ChIPseqSpikeInFree-bigwig-nf: Output

```
<outdir>/
├── bigwig_qc/                               # CSIF_BIN_BIGWIG
│   └── <sample>_bigwig_qc.json              # per-sample binning summary (see below)
├── counts/                                  # CSIF_MERGE_COUNTS + CSIF_PARSE
│   ├── <cohort>_rawCounts.txt               # pseudo read counts per genome-wide bin per sample
│   ├── <cohort>_merge.log                   # samples merged / dropped
│   ├── <cohort>_sample_meta.txt             # metadata generated from the samplesheet (ID = <sample>.bw)
│   ├── <cohort>_parsedMatrix.txt            # cumulative read distributions used for calibration
│   └── <cohort>_parse.log
├── scaling_factors/                         # CSIF_SCALE
│   ├── <cohort>_SF.txt                      # scaling factors and QC calls per sample
│   ├── <cohort>_distribution.pdf            # cumulative distribution curves with turning points
│   ├── <cohort>_boxplot.pdf                 # scaling factors by group
│   └── <cohort>_scale.log
└── pipeline_info/
```

`<cohort>_rawCounts.txt` holds **pseudo-counts**: bigWig signal mass per bin, rescaled to the sample's
`reads` (or `--pseudo_depth`). They are not read counts; only their relative values mean anything.

## `<cohort>_SF.txt`

As in the BAM pipeline, with `ID = <sample>.bw`: `GROUP`, `ANTIBODY`, `COLOR`, `QC` (`pass` or the
reason a sample failed), `TURNS` (turning points used), `SF` (scaling factor, `NA` for failed samples).
To apply a factor, multiply the sample's normalised signal by its `SF` (see the ChIPseqSpikeInFree
README).

## `bigwig_qc/<sample>_bigwig_qc.json`

| Field | Meaning |
| ----- | ------- |
| `n_bins`, `n_bins_with_signal`, `fraction_bins_with_signal` | bins built / bins with any signal |
| `total_signal_mass`, `mean_signal_per_bp` | sum of per-base values over all bins; per bin base |
| `pseudo_depth`, `pseudo_depth_source` | pseudo read total used and whether it came from `reads` |
| `fraction_bins_cpm_gt_150` | fraction of bins above 150 CPM; ChIPseqSpikeInFree zeroes such bins |
| `chroms_used`, `chroms_missing_in_bigwig`, `chroms_length_mismatch` | chromosome matching against the genome |
| `summary_statistics` | `exact` or zoom-level (`--approx_summary`) |
| `remote`, `http_requests`, `http_bytes_fetched` | whether the file was read over HTTP, requests made, bytes transferred |
