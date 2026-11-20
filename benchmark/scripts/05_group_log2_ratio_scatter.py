#!/usr/bin/env python3
"""Group-mean log2 ratio scatter (like panel E of the ChIPseqSpikeInFree paper), for the benchmark results.
Ratio = group-mean scaling factor of the mutant (K27M, K36M) over the control (WT, K27M-KO or K27R).
a: all 10 comparisons of btz720 Supplementary Table 1, spike-in vs published ChIPseqSpikeInFree (BAM route).
b: the 4 comparisons of GSE108344 with the bigWig route (this pipeline) next to the published BAM route.
Inputs: out_bm/dataset1_scaling_factors_three_methods.csv, bench/dataset1_bigwig_csif_SF_jitter_long.csv, btz720_supp.txt (Table 1 as text)."""
import re, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
try: apply_figure_style()
except NameError: pass
OUT = 'out_bm'
# ---- all comparisons of Table 1 (datasets 1-4)
txt = open('btz720_supp.txt', encoding='utf-8').read(); tab = txt[txt.index('<<TABLE 1>>'):]
rows, ds = [], None
for line in tab.splitlines():
    m = re.match(r'Dataset (\d) \((.+?)\)', line)
    if m: ds = (int(m[1]), m[2]); continue
    m = re.match(r'(\S+) \| (\S+) \| ([\d.]+) \| ([\d.]+) \| (\S+)$', line)
    if m and ds: rows.append(dict(dataset=ds[0], ref=ds[1], sample=m[1], group=m[2], sf_spike=float(m[3]), sf_csif=float(m[4]), mark=m[5]))
T = pd.DataFrame(rows); T['geno'] = T.group.str.split('.').str[-1]; T['system'] = T.group.str.split('.').str[0]; T = T[T.dataset <= 4]
comps = []
for (ds_, sy, mk_), g in T.groupby(['dataset', 'system', 'mark']):
    a = g[g.geno.isin({'K27M', 'K36M'})]; b = g[g.geno.isin({'WT', 'K27M-KO', 'K27R'})]
    if len(a) and len(b):
        comps.append(dict(dataset=ds_, ref=g.ref.iloc[0], system=sy, mark=mk_, mutant=a.geno.iloc[0], control=b.geno.iloc[0], n_mutant=len(a), n_control=len(b),
                          log2_spike=np.log2(a.sf_spike.mean() / b.sf_spike.mean()), log2_csif_bam=np.log2(a.sf_csif.mean() / b.sf_csif.mean())))
C = pd.DataFrame(comps); C.round(3).to_csv(f'{OUT}/table1_all_comparisons_log2_ratios.csv', index=False)
# ---- GSE108344 groups, three methods + numeric-noise range
d = pd.read_csv(f'{OUT}/dataset1_scaling_factors_three_methods.csv'); jt = pd.read_csv('bench/dataset1_bigwig_csif_SF_jitter_long.csv')
geno = d.set_index('gsm')[['cell_system', 'mark', 'genotype']]
lr = lambda g, col: np.log2(g[g.genotype == 'K27M'][col].mean() / g[g.genotype == 'WT'][col].mean())
rows = []
for (cs, mk_), g in d.groupby(['cell_system', 'mark']):
    x = jt[jt.eps == 1e-14].merge(geno, left_on='gsm', right_index=True); x = x[(x.cell_system == cs) & (x.mark == mk_)]
    rep = x.groupby(['rep', 'genotype']).sf.mean().unstack(); rr = np.log2(rep['K27M'] / rep['WT'])
    rows.append(dict(cell_system=cs, mark=mk_, n_K27M=int((g.genotype == 'K27M').sum()), n_WT=int((g.genotype == 'WT').sum()),
                     log2_spikein=lr(g, 'sf_spikein'), log2_csif_bam_published=lr(g, 'sf_csif_bam_published'),
                     log2_csif_bigwig_pipeline=lr(g, 'sf_bigwig_pipeline'), log2_csif_bigwig_single_run=lr(g, 'sf_bigwig_single_run'),
                     log2_csif_bigwig_noise_min=rr.min(), log2_csif_bigwig_noise_max=rr.max(),
                     log2_csif_bigwig_minFirstTurn0p5_exploratory=lr(g, 'sf_bigwig_minFirstTurn0.5_exploratory')))
G = pd.DataFrame(rows).sort_values(['cell_system', 'mark']).reset_index(drop=True); G.round(3).to_csv(f'{OUT}/dataset1_group_log2_ratios_three_methods.csv', index=False)
rA = stats.pearsonr(C.log2_spike, C.log2_csif_bam); rB = stats.pearsonr(G.log2_spikein, G.log2_csif_bigwig_pipeline); rC = stats.pearsonr(G.log2_spikein, G.log2_csif_bam_published)
fp = lambda p: f'{p:.1e}' if p < 1e-3 else f'{p:.3f}'.rstrip('0').rstrip('.') if p < 0.1 else f'{p:.2f}'
# ---- figure
mk = {'H3K27Ac': 'D', 'H3K27me3': 'o', 'H3K4me3': 's', 'H3K36me2': '^', 'H3K36me3': 'v'}
dc = {1: '#0072B2', 2: '#E69F00', 3: '#009E73', 4: '#CC79A7'}
dn = {1: 'Larson 2019 (GSE108344)', 2: 'Pathania 2017', 3: 'Harutyunyan 2019', 4: 'Fang 2016'}
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.6, 3.7), gridspec_kw=dict(wspace=0.3))
for a_, lim in ((ax, (-1.2, 5.2)), (bx, (-1.2, 3.0))):
    a_.plot(lim, lim, ls='--', color='0.6', lw=0.9, zorder=1); a_.axhline(0, color='0.85', lw=0.6, zorder=0); a_.axvline(0, color='0.85', lw=0.6, zorder=0)
    a_.set_xlim(lim); a_.set_ylim(lim); a_.set_aspect('equal'); a_.set_xlabel('Spike-in, log$_2$ ratio'); a_.set_ylabel('ChIPseqSpikeInFree, log$_2$ ratio')
for r in C.itertuples():
    ax.scatter(r.log2_spike, r.log2_csif_bam, marker=mk[r.mark], s=38, color=dc[r.dataset], edgecolor='black', linewidth=0.6, zorder=3)
ax.set_title('Published BAM route, 10 comparisons', loc='left')
ax.text(0.04, 0.95, f'r = {rA[0]:.2f}, p = {fp(rA[1])}\nn = {len(C)} comparisons', transform=ax.transAxes, va='top', fontsize=7)
ax.legend(handles=[Line2D([], [], marker=m, ls='', color='white', markeredgecolor='black', markerfacecolor='0.8', markersize=5.5, label=k) for k, m in mk.items()],
          loc='lower right', frameon=False, fontsize=6.5, handletextpad=0.3, borderaxespad=0.2)
for r in G.itertuples():
    bx.vlines(r.log2_spikein, r.log2_csif_bigwig_noise_min, r.log2_csif_bigwig_noise_max, color=dc[1], lw=1.1, alpha=0.5, zorder=2)
    bx.scatter(r.log2_spikein, r.log2_csif_bigwig_pipeline, marker=mk[r.mark], s=40, color=dc[1], edgecolor='black', linewidth=0.6, zorder=4)
    bx.scatter(r.log2_spikein, r.log2_csif_bam_published, marker=mk[r.mark], s=40, facecolor='white', edgecolor=dc[1], linewidth=1.1, zorder=3)
off = {'DIPG H3K27Ac': (0, 13, 'center'), 'DIPG H3K27me3': (-9, 3, 'right'), 'DIPG H3K4me3': (8, -2, 'left'), 'H-NSC H3K27me3': (-9, -9, 'right')}
for r in G.itertuples():
    key = f'{r.cell_system} {r.mark}'; ox, oy, ha = off[key]
    bx.annotate(key, (r.log2_spikein, r.log2_csif_bigwig_pipeline), xytext=(ox, oy), textcoords='offset points', fontsize=6.5, ha=ha, va='center')
bx.set_title('bigWig route, 4 GSE108344 comparisons', loc='left')
bx.text(0.04, 0.95, f'bigWig route: r = {rB[0]:.2f}, p = {fp(rB[1])}\nBAM route: r = {rC[0]:.2f}, p = {fp(rC[1])}\nn = 4 comparisons', transform=bx.transAxes, va='top', fontsize=6.5)
bx.legend(handles=[Line2D([], [], marker='o', ls='', markerfacecolor=dc[1], markeredgecolor='black', markersize=5.5, label='bigWig route (this pipeline)'),
                   Line2D([], [], marker='o', ls='', markerfacecolor='white', markeredgecolor=dc[1], markeredgewidth=1.1, markersize=5.5, label='BAM route (published)')],
          loc='lower right', frameon=False, fontsize=6.5, handletextpad=0.3, borderaxespad=0.2)
ax.text(-0.2, 1.04, 'a', transform=ax.transAxes, fontweight='bold', fontsize=11); bx.text(-0.2, 1.04, 'b', transform=bx.transAxes, fontweight='bold', fontsize=11)
fig.legend(handles=[Line2D([], [], marker='o', ls='', markerfacecolor=dc[k], markeredgecolor='black', markersize=5.5, label=dn[k]) for k in dn],
           loc='lower center', ncol=4, frameon=False, fontsize=6.5, handletextpad=0.3, columnspacing=1.5, bbox_to_anchor=(0.5, -0.04))
fig.text(0.5, -0.09, 'Ratio of group-mean scaling factors, mutant (K27M, K36M) over control (WT, K27M-KO or K27R). Vertical bars in b: range over 50 recalculations with 1e-14 relative numerical noise.', ha='center', fontsize=6, color='0.35')
fig.savefig(f'{OUT}/fig_group_log2_ratio_scatter.png', dpi=300, bbox_inches='tight')
print(f'a: r={rA[0]:.3f} p={rA[1]:.3g} n={len(C)} | b: bigWig r={rB[0]:.3f} p={rB[1]:.3g}; BAM r={rC[0]:.3f} p={rC[1]:.3g}')
