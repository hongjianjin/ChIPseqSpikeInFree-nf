#!/usr/bin/env python3
"""Compare scaling factors for the 33 samples of Dataset 1 (GSE108344) in btz720 Supplementary Table 1:
spike-in, ChIPseqSpikeInFree on BAM files (published), ChIPseqSpikeInFree on bigWig tracks (this pipeline).

Inputs (working layout used for the benchmark):
  dataset1_table1_bw.csv                           Table 1 values + GEO bigWig URLs
  bench/dataset1_bigwig_csif_SF.csv                bigWig route, single in-memory run (ParseReadCounts -> CalculateSF)
  bench/dataset1_bigwig_csif_SF_from_file.csv      bigWig route, parsedMatrix re-read from file (what CSIF_SCALE does)
  bench/dataset1_bigwig_csif_SF_jitter_long.csv    CalculateSF on the parsed matrix with 1e-14 / 1e-10 relative noise, 50 replicates
  bench/dataset1_bigwig_csif_SF_minFirstTurn0.5.csv  exploratory: minFirstTurn fixed at 0.5
  bench/counts/<GSM>_bigwig_qc.json                binning QC
Outputs go to out_bm/.
"""
import glob, json, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

OUT = 'out_bm'; os.makedirs(OUT, exist_ok=True)
t1 = pd.read_csv('dataset1_table1_bw.csv')
mem = pd.read_csv('bench/dataset1_bigwig_csif_SF.csv')[['gsm', 'SF', 'QC']].rename(columns={'SF': 'sf_bigwig_single_run', 'QC': 'qc_bigwig_single_run'})
fil = pd.read_csv('bench/dataset1_bigwig_csif_SF_from_file.csv')[['gsm', 'sf_file', 'QC_file']].rename(columns={'sf_file': 'sf_bigwig_pipeline', 'QC_file': 'qc_bigwig_pipeline'})
jit = pd.read_csv('bench/dataset1_bigwig_csif_SF_jitter_long.csv')
fx = pd.read_csv('bench/dataset1_bigwig_csif_SF_minFirstTurn0.5.csv')[['gsm', 'sf_fixed0.5']].rename(columns={'sf_fixed0.5': 'sf_bigwig_minFirstTurn0.5_exploratory'})
j14 = jit[jit.eps == 1e-14].groupby('gsm').sf.agg(sf_bigwig_noise_min='min', sf_bigwig_noise_median='median', sf_bigwig_noise_max='max', n_distinct_values='nunique').reset_index()

d = (t1.rename(columns={'sf_spike': 'sf_spikein', 'sf_csif_bam': 'sf_csif_bam_published', 'title': 'geo_sample_title'})
       .merge(mem, on='gsm').merge(fil, on='gsm').merge(j14, on='gsm').merge(fx, on='gsm'))
d['cell_system'] = d.group.str.split('.').str[0]
d['genotype'] = d.group.str.split('.').str[-1]
d['run'] = d.cell_system + ' ' + d.mark
for a, b, name in [('sf_bigwig_pipeline', 'sf_csif_bam_published', 'log2_bigwig_over_bam'), ('sf_bigwig_pipeline', 'sf_spikein', 'log2_bigwig_over_spikein'), ('sf_csif_bam_published', 'sf_spikein', 'log2_bam_over_spikein')]:
    d[name] = np.log2(d[a] / d[b])
cols = ['gsm', 'geo_sample_title', 'cell_system', 'mark', 'genotype', 'sf_spikein', 'sf_csif_bam_published',
        'sf_bigwig_pipeline', 'sf_bigwig_single_run', 'sf_bigwig_noise_min', 'sf_bigwig_noise_median', 'sf_bigwig_noise_max', 'n_distinct_values',
        'sf_bigwig_minFirstTurn0.5_exploratory', 'qc_bigwig_pipeline', 'qc_bigwig_single_run',
        'log2_bigwig_over_bam', 'log2_bigwig_over_spikein', 'log2_bam_over_spikein', 'https']
d = d.sort_values(['cell_system', 'mark', 'genotype', 'gsm'])
d[cols].rename(columns={'https': 'bigwig_url', 'n_distinct_values': 'n_distinct_sf_in_50_noise_replicates'}).round(3).to_csv(f'{OUT}/dataset1_scaling_factors_three_methods.csv', index=False)

# ---- summary statistics
def pr(a, b): return stats.pearsonr(np.log2(a), np.log2(b))[0]
def sr(a, b): return stats.spearmanr(a, b)[0]
rows = []
for label, col in [('bigWig, pipeline path (CalculateSF on re-read matrix)', 'sf_bigwig_pipeline'), ('bigWig, single in-memory run', 'sf_bigwig_single_run'), ('bigWig, minFirstTurn=0.5 (exploratory)', 'sf_bigwig_minFirstTurn0.5_exploratory')]:
    rows.append(dict(comparison=f'{label} vs published BAM route', n=len(d), pearson_r_log2=pr(d[col], d.sf_csif_bam_published), spearman_rho=sr(d[col], d.sf_csif_bam_published),
                     median_abs_log2_diff=np.median(np.abs(np.log2(d[col] / d.sf_csif_bam_published)))))
    rows.append(dict(comparison=f'{label} vs spike-in', n=len(d), pearson_r_log2=pr(d[col], d.sf_spikein), spearman_rho=sr(d[col], d.sf_spikein),
                     median_abs_log2_diff=np.median(np.abs(np.log2(d[col] / d.sf_spikein)))))
rows.append(dict(comparison='published BAM route vs spike-in', n=len(d), pearson_r_log2=pr(d.sf_csif_bam_published, d.sf_spikein), spearman_rho=sr(d.sf_csif_bam_published, d.sf_spikein),
                 median_abs_log2_diff=np.median(np.abs(np.log2(d.sf_csif_bam_published / d.sf_spikein)))))
pd.DataFrame(rows).round(3).to_csv(f'{OUT}/dataset1_summary_statistics.csv', index=False)

g_rows = []
for run, g in d.groupby('run'):
    k, w = g[g.genotype == 'K27M'], g[g.genotype == 'WT']
    g_rows.append(dict(group=run, n=len(g), rho_bigwig_vs_bam=sr(g.sf_bigwig_pipeline, g.sf_csif_bam_published), rho_bigwig_vs_spikein=sr(g.sf_bigwig_pipeline, g.sf_spikein),
                       rho_bam_vs_spikein=sr(g.sf_csif_bam_published, g.sf_spikein),
                       fold_K27M_over_WT_spikein=k.sf_spikein.mean() / w.sf_spikein.mean(), fold_K27M_over_WT_bam=k.sf_csif_bam_published.mean() / w.sf_csif_bam_published.mean(),
                       fold_K27M_over_WT_bigwig=k.sf_bigwig_pipeline.mean() / w.sf_bigwig_pipeline.mean(),
                       median_abs_log2_bigwig_vs_bam=np.median(np.abs(g.log2_bigwig_over_bam))))
pd.DataFrame(g_rows).round(3).to_csv(f'{OUT}/dataset1_per_group_summary.csv', index=False)

# ---- numeric jitter per sample (both noise levels)
jj = jit.groupby(['gsm', 'eps']).sf.agg(n_distinct='nunique', sf_min='min', sf_median='median', sf_max='max').reset_index()
jj['max_over_min'] = jj.sf_max / jj.sf_min
jj.round(3).to_csv(f'{OUT}/dataset1_bigwig_numeric_jitter_per_sample.csv', index=False)

# ---- Gate 1 criteria under each realization and under noise
def gate(sf, ref): 
    l = np.log2(sf); r = np.log2(ref); return stats.pearsonr(l, r)[0], np.median(np.abs(l - r))
ref = d.set_index('gsm').sf_csif_bam_published
g_rows = []
for name, col in [('single in-memory run (analysis fixed in advance)', 'sf_bigwig_single_run'), ('pipeline path (matrix re-read from file)', 'sf_bigwig_pipeline'), ('minFirstTurn=0.5 (exploratory)', 'sf_bigwig_minFirstTurn0.5_exploratory')]:
    r_, m_ = gate(d[col].values, d.sf_csif_bam_published.values)
    qc = (d['qc_bigwig_pipeline' if 'pipeline' in col else 'qc_bigwig_single_run'] != 'pass').sum() if col != 'sf_bigwig_minFirstTurn0.5_exploratory' else np.nan
    g_rows.append(dict(realization=name, P1_pearson_r=r_, P1_pass=r_ >= 0.90, P2_median_abs_log2=m_, P2_pass=m_ <= 0.15, P3_qc_failures=qc))
for eps in (1e-14, 1e-10):
    rs, ms = [], []
    for rep, g in jit[jit.eps == eps].groupby('rep'):
        x = g.set_index('gsm').sf.reindex(d.gsm).values
        r_, m_ = gate(x, d.sf_csif_bam_published.values); rs.append(r_); ms.append(m_)
    g_rows.append(dict(realization=f'50 noise replicates, relative noise {eps:g}: range', P1_pearson_r=f'{min(rs):.3f}-{max(rs):.3f}', P1_pass=f'{sum(r >= 0.9 for r in rs)}/50 pass',
                       P2_median_abs_log2=f'{min(ms):.3f}-{max(ms):.3f}', P2_pass=f'{sum(m <= 0.15 for m in ms)}/50 pass', P3_qc_failures=int((jit[jit.eps == eps].qc != 'pass').sum())))
pd.DataFrame(g_rows).to_csv(f'{OUT}/gate1_criteria_by_realization.csv', index=False)

# ---- binning QC table
q = [json.load(open(f)) for f in sorted(glob.glob('bench/counts/*_bigwig_qc.json'))]
qd = pd.DataFrame([{k: v for k, v in x.items() if k in ('id', 'n_bins', 'fraction_bins_with_signal', 'total_signal_mass', 'fraction_bins_cpm_gt_150', 'chroms_used', 'http_requests', 'http_bytes_fetched', 'seconds')} for x in q])
qd['id'] = qd.id.str.replace('.bw', '', regex=False)
qd.rename(columns={'id': 'gsm'}).round(5).to_csv(f'{OUT}/dataset1_bigwig_binning_qc.csv', index=False)

# ---- figure: sample-level scatter of the three scaling-factor sets
try:
    apply_figure_style()
except NameError:
    pass
runs = ['DIPG H3K27me3', 'H-NSC H3K27me3', 'DIPG H3K27Ac', 'DIPG H3K4me3']
colr = {'DIPG H3K27me3': '#0072B2', 'H-NSC H3K27me3': '#56B4E9', 'DIPG H3K27Ac': '#D55E00', 'DIPG H3K4me3': '#E69F00'}   # blue family: H3K27me3; orange family: other marks
fig, axs = plt.subplots(1, 3, figsize=(10.2, 3.5), gridspec_kw=dict(wspace=0.34))
lim = (0.8, 32); ticks = [1, 2, 4, 8, 16]
panels = [('sf_spikein', 'sf_csif_bam_published', 'Spike-in scaling factor', 'ChIPseqSpikeInFree, BAM (published)', 'BAM route vs spike-in', False),
          ('sf_spikein', 'sf_bigwig_pipeline', 'Spike-in scaling factor', 'ChIPseqSpikeInFree, bigWig (this pipeline)', 'bigWig route vs spike-in', True),
          ('sf_csif_bam_published', 'sf_bigwig_pipeline', 'ChIPseqSpikeInFree, BAM (published)', 'ChIPseqSpikeInFree, bigWig (this pipeline)', 'bigWig route vs BAM route', True)]
for k, (xc, yc, xl, yl, ttl, bars) in enumerate(panels):
    ax = axs[k]
    ax.plot(lim, lim, color='0.6', lw=0.8, zorder=1)
    for r in runs:
        g = d[d.run == r]
        if bars:
            lo, hi = g.sf_bigwig_noise_min, g.sf_bigwig_noise_max
            ax.vlines(g[xc], lo, hi, color=colr[r], lw=0.9, alpha=0.45, zorder=2)
        ax.scatter(g[xc], g[yc], s=20, color=colr[r], edgecolor='white', linewidth=0.4, zorder=3, label=r)
    ax.set_xscale('log', base=2); ax.set_yscale('log', base=2)
    for a_ in (ax.xaxis, ax.yaxis):
        a_.set_major_locator(FixedLocator(ticks)); a_.set_major_formatter(FixedFormatter([str(t) for t in ticks])); a_.set_minor_locator(FixedLocator([]))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(ttl, loc='left')
    ax.text(0.97, 0.05, f'r = {pr(d[xc], d[yc]):.2f} (log$_2$ SF)\nρ = {sr(d[xc], d[yc]):.2f}\nn = {len(d)} samples', transform=ax.transAxes, ha='right', va='bottom', fontsize=7)
    if k == 2:
        top = d.loc[d.sf_bigwig_pipeline.idxmax()]
        ax.annotate(top.gsm, (top[xc], top[yc]), xytext=(-46, -2), textcoords='offset points', fontsize=6, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc='lower center', ncol=4, frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.6, bbox_to_anchor=(0.5, -0.07))
fig.text(0.5, -0.135, 'Vertical bars in b and c: range of the bigWig-route factor in 50 recalculations with 1e-14 relative numerical noise', ha='center', fontsize=6.5, color='0.35')
try:
    for ax, L in zip(axs, 'abc'): panel_letter(ax, L)
except NameError:
    for ax, L in zip(axs, 'abc'): ax.text(-0.18, 1.06, L, transform=ax.transAxes, fontweight='bold', fontsize=11)
fig.savefig(f'{OUT}/fig_scaling_factors_three_methods.png', dpi=300, bbox_inches='tight')
print(pd.read_csv(f'{OUT}/dataset1_summary_statistics.csv').to_string(index=False))
print(pd.read_csv(f'{OUT}/gate1_criteria_by_realization.csv').to_string(index=False))
print(pd.read_csv(f'{OUT}/dataset1_per_group_summary.csv').to_string(index=False))
