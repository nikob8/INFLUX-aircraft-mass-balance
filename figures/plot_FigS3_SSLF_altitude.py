#!/usr/bin/env python3
"""Figure S3: dependence of the South Side Landfill (SSLF) estimates on transect
altitude (Section 3.1.3 of the paper; Reviewer 2, Line 512, and Reviewer 1).

Three panels, built from Table S1 (one row per downwind transect):
  (a) SSLF flux deviation vs. transect-height deviation, with the daily means
      removed as in Figure 8, for the days with at least two quantified plumes;
      linear fit with slope +/- 1 SE and its p-value
  (b) for the days on which the plume was quantified on some transects but
      not on others: within-day height deviation of the transects with and
      without a plume (mean +/- SE), with a rank-sum p-value
  (c) fraction of all transects with a quantified plume, binned by transect
      height over ABL height, with 95% Wilson intervals and a logistic fit
The script also prints (d): the citywide Downwind Edge emission rate on the
transects with and without a quantified plume of the same days (Section 3.1.3).

Input (same folder): Indy_Mass_Balance_Table_v5.csv
Output: FigS3_SSLF_altitude.png and the statistics printed.
Requires numpy, pandas, scipy, matplotlib.
"""
import numpy as np, pandas as pd
from scipy import stats
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t = pd.read_csv('Indy_Mass_Balance_Table_v5.csv')
t['det'] = t.Flux_Landfill_mols.notna()
t['h_dev'] = t.Height_AGL_m - t.groupby('Date').Height_AGL_m.transform('mean')    # metres from the mean height of the day
t['z_zi'] = t.Height_AGL_m / t.ABL_height_m

# ---- (a) flux vs height, daily means removed, days with >= 2 quantified plumes
q = t[t.det].copy()
nq = q.groupby('Date').Transect.transform('count')
a = q[nq >= 2].copy()
a['f_dev'] = a.Flux_Landfill_mols - a.groupby('Date').Flux_Landfill_mols.transform('mean')
a['hq_dev'] = a.Height_AGL_m - a.groupby('Date').Height_AGL_m.transform('mean')   # relative to the quantified transects of the day
fit = stats.linregress(a.hq_dev, a.f_dev)
print(f'(a) n = {len(a)} transects on {a.Date.nunique()} days; slope = {100*fit.slope:+.2f} +/- {100*fit.stderr:.2f} mol/s per 100 m, p = {fit.pvalue:.2f}')
a['zq_dev'] = a.z_zi - a.groupby('Date').z_zi.transform('mean')
fz = stats.linregress(a.zq_dev, a.f_dev)
print(f'    vs height/ABL: slope = {fz.slope:+.1f} +/- {fz.stderr:.1f} mol/s per unit z/zi, p = {fz.pvalue:.2f}')

# ---- (b) days with both quantified and non-quantified transects
g = t.groupby('Date').det.agg(['sum', 'count'])
mixed = g.index[(g['sum'] > 0) & (g['sum'] < g['count'])]
b = t[t.Date.isin(mixed)]
hd, hn = b.h_dev[b.det], b.h_dev[~b.det]
pb = stats.mannwhitneyu(hd, hn, alternative='two-sided').pvalue
print(f'(b) {len(mixed)} days; with plume mean {hd.mean():+.0f} m (n = {len(hd)}), without {hn.mean():+.0f} m (n = {len(hn)}); difference {hn.mean()-hd.mean():.0f} m; rank-sum p = {pb:.1e}')

# ---- (d) is the landfill CH4 missing on the transects without a quantified plume?
#      Citywide Downwind Edge emission rate on the same days, transects without minus with a
#      quantified plume (mean over days of the within-day difference, t-based 95% interval).
#      If the landfill plume were absent from those transects, the difference would be close
#      to minus the SSLF emission rate.
dd = np.array([s.Flux_DownwindEdge_mols[~s.det].mean() - s.Flux_DownwindEdge_mols[s.det].mean() for _, s in b.groupby('Date')])
f_with = np.mean([s.Flux_DownwindEdge_mols[s.det].mean() for _, s in b.groupby('Date')])
f_without = np.mean([s.Flux_DownwindEdge_mols[~s.det].mean() for _, s in b.groupby('Date')])
f_lf = b.groupby('Date').Flux_Landfill_mols.mean().mean()
se_d = dd.std(ddof=1) / np.sqrt(len(dd)); tq = stats.t.ppf(0.975, len(dd) - 1)
print(f'(d) citywide Downwind Edge flux on the same {len(dd)} days: {f_with:.1f} mol/s on transects with a quantified plume, '
      f'{f_without:.1f} without; mean within-day difference {dd.mean():+.1f} mol/s (95% CI {dd.mean() - tq*se_d:+.1f} to '
      f'{dd.mean() + tq*se_d:+.1f}); SSLF emission rate on these days {f_lf:.1f} mol/s')

# ---- (c) detection fraction vs z/zi
edges = [0, 0.25, 0.5, 0.75, 1.0, 1.5]
c = t.dropna(subset=['z_zi']).copy()
c['bin'] = pd.cut(c.z_zi, edges, right=False)
frac = c.groupby('bin', observed=True).det.agg(['mean', 'sum', 'count'])
def wilson(k, n, z=1.96):
    p = k/n; d = 1 + z*z/n; centre = (p + z*z/(2*n))/d; half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return centre - half, centre + half
lo, hi = zip(*[wilson(k, n) for k, n in zip(frac['sum'], frac['count'])])
# logistic regression det ~ z/zi (maximum likelihood, Wald test on the slope)
X = np.column_stack([np.ones(len(c)), c.z_zi.values]); y = c.det.values.astype(float)
nll = lambda beta: -np.sum(y*(X@beta) - np.log1p(np.exp(X@beta)))
res = minimize(nll, np.zeros(2), method='BFGS')
beta = res.x; pr = 1/(1 + np.exp(-(X@beta))); W = pr*(1-pr)
cov = np.linalg.inv((X*W[:, None]).T @ X); se = np.sqrt(np.diag(cov)); pz = 2*stats.norm.sf(abs(beta[1]/se[1]))
print('(c) detection fraction by z/zi bin:', ' '.join(f'{m:.2f} (n={n})' for m, n in zip(frac['mean'], frac['count'])),
      f'; logistic slope = {beta[1]:.2f} +/- {se[1]:.2f}, p = {pz:.2f}; above the ABL top: {int(frac["sum"].iloc[-1])} of {int(frac["count"].iloc[-1])}')
low = c[c.z_zi < 0.5].det.mean(); mid = c[(c.z_zi >= 0.75) & (c.z_zi < 1.0)].det.mean()
print(f'    lower half of the ABL: {100*low:.0f}%; 0.75 to 1.0: {100*mid:.0f}%')

# ---- figure
cB, cR = '#2f5597', '#c0504d'
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19.5, 6.3))
ax1.axhline(0, color='k', lw=0.8)
ax1.scatter(a.hq_dev, a.f_dev, s=50, c=cB, zorder=3)
xx = np.array([-650, 650]); ax1.plot(xx, fit.intercept + fit.slope*xx, 'r-', lw=2, label=f'Slope = ({100*fit.slope:.2f} ± {100*fit.stderr:.2f}) per 100 m, p = {fit.pvalue:.2f}')
ax1.set_xlim(-650, 650); ax1.set_xlabel('Meters from mean transect height', fontsize=14); ax1.set_ylabel('SSLF flux deviation (mol s$^{-1}$)', fontsize=14)
ax1.set_title(f'(a) SSLF flux vs. altitude (n = {len(a)} transects, {a.Date.nunique()} days)', fontsize=14); ax1.legend(loc='upper right', fontsize=12)

rng = np.random.default_rng(1)
for x0, v, col, lab in [(0, hd, cB, f'Plume detected\n(n = {len(hd)})'), (1, hn, cR, f'No plume\n(n = {len(hn)})')]:
    ax2.scatter(x0 + rng.uniform(-0.18, 0.18, len(v)), v, s=40, c=col, alpha=0.85, zorder=3)
    ax2.errorbar(x0, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)), fmt='_', color='k', ms=40, mew=4, capsize=6, zorder=4)
ax2.axhline(0, color='k', lw=0.8); ax2.set_xticks([0, 1]); ax2.set_xticklabels([f'Plume detected\n(n = {len(hd)})', f'No plume\n(n = {len(hn)})'], fontsize=13)
ax2.set_xlim(-0.6, 1.6); ax2.set_ylabel('Meters from mean transect height', fontsize=14)
ax2.set_title(f'(b) Same-day transect height ({len(mixed)} days), p = {pb:.0e}', fontsize=14)

centres = [(edges[i] + edges[i+1])/2 for i in range(len(edges)-1)]; widths = np.diff(edges)
ax3.bar(centres, frac['mean'], width=widths*0.96, color='#9db3d6', edgecolor='k', zorder=2)
ax3.errorbar(centres, frac['mean'], yerr=[frac['mean'] - np.array(lo), np.array(hi) - frac['mean']], fmt='none', ecolor='k', capsize=5, lw=1.6, zorder=3)
for xc, n in zip(centres, frac['count']): ax3.text(xc, 0.03, f'n={n}', ha='center', fontsize=11)
zz = np.linspace(0, 1.5, 100); ax3.plot(zz, 1/(1 + np.exp(-(beta[0] + beta[1]*zz))), 'r-', lw=2, label=f'Logistic fit, p = {pz:.2f}')
ax3.set_xlim(0, 1.5); ax3.set_ylim(0, 1); ax3.set_xlabel('Transect height / ABL height', fontsize=14); ax3.set_ylabel('Fraction of transects with detectable plume', fontsize=14)
ax3.set_title(f'(c) Plume detection vs. height in ABL (all {len(c)} transects)', fontsize=14); ax3.legend(loc='upper right', fontsize=12)
for ax in (ax1, ax2, ax3): ax.tick_params(labelsize=12)
fig.tight_layout()
fig.savefig('FigS3_SSLF_altitude.png', dpi=200)
print('saved FigS3_SSLF_altitude.png')
