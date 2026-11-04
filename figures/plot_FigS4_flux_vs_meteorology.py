#!/usr/bin/env python3
"""Figure S4: daily citywide emission estimates (Downwind Edge and Upwind Edge)
against the mean ABL height of the downwind transects, the mean wind speed and
the mean wind direction of the flight day (Section 3.1.1 of the paper).

Inputs (same folder):
    Indy_Mass_Balance_Table_v5.csv   per-transect results (Table S1)
    diagnostic_ratios.csv            diagnostic ratio of each flight day (Section 3.2)
Output: FigS4_flux_vs_meteorology.png and the correlation statistics printed.
Requires numpy, pandas, scipy, matplotlib.
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t = pd.read_csv('Indy_Mass_Balance_Table_v5.csv')
r = pd.read_csv('diagnostic_ratios.csv')
day = (t.groupby('Date')
         .agg(A=('Flux_DownwindEdge_mols', 'mean'), UE=('Flux_UpwindEdge_mols', 'mean'),
              zi=('ABL_height_m', 'mean'), ws=('Mean_wind_speed_ms', 'first'), wd=('Mean_wind_direction_deg', 'first'))
         .reset_index()
         .merge(r[['Date', 'Diagnostic_ratio']], on='Date'))
good = day.Diagnostic_ratio.values > 2


def pearson(x, y):
    ok = ~np.isnan(x) & ~np.isnan(y)
    rr, p = stats.pearsonr(x[ok], y[ok])
    return rr, p, ok.sum()


def circ_lin(theta_deg, y):
    """circular-linear correlation (Mardia and Jupp) with the chi-square p-value"""
    ok = ~np.isnan(theta_deg) & ~np.isnan(y)
    th, yy = np.deg2rad(theta_deg[ok]), y[ok]
    rxc = stats.pearsonr(yy, np.cos(th))[0]
    rxs = stats.pearsonr(yy, np.sin(th))[0]
    rcs = stats.pearsonr(np.cos(th), np.sin(th))[0]
    R = np.sqrt((rxc**2 + rxs**2 - 2*rxc*rxs*rcs) / (1 - rcs**2))
    return R, 1 - stats.chi2.cdf(ok.sum()*R**2, 2), ok.sum()


cG, cB = '#1f4e96', '#c0392b'
fig, axs = plt.subplots(2, 3, figsize=(15, 8.9))
rows = [('A', 'Downwind Edge emission rate (mol s$^{-1}$)'), ('UE', 'Upwind Edge emission rate (mol s$^{-1}$)')]
cols = [('zi', 'ABL height (m)'), ('ws', 'Mean wind speed (m s$^{-1}$)'), ('wd', 'Mean wind direction')]
letters = 'abcdef'
k = 0
for i, (fy, ylab) in enumerate(rows):
    y = day[fy].values
    for j, (fx, xlab) in enumerate(cols):
        ax = axs[i, j]; x = day[fx].values
        ok = ~np.isnan(y)
        ax.scatter(x[ok & good], y[ok & good], s=70, c=cG, edgecolor='k', linewidth=0.5, zorder=3, label='diagnostic ratio > 2')
        ax.scatter(x[ok & ~good], y[ok & ~good], s=70, facecolor='none', edgecolor=cB, linewidth=1.6, zorder=3, label='diagnostic ratio $\\leq$ 2')
        if fx == 'wd':
            Ra, pa, na = circ_lin(x, y); Rg, pg, ng = circ_lin(x[good], y[good])
            txt = f'all days: r = {Ra:.2f}, p = {pa:.2f}\nratio > 2: r = {Rg:.2f}, p = {pg:.2f}'
            ax.set_xlim(0, 360); ax.set_xticks([0, 90, 180, 270, 360]); ax.set_xticklabels(['N', 'E', 'S', 'W', 'N'])
            print(f'{fy} vs wind direction: circular-linear r = {Ra:.2f} (p = {pa:.2f}, n = {na}); ratio > 2: r = {Rg:.2f} (p = {pg:.2f}, n = {ng})')
        else:
            ra, pa, na = pearson(x, y); rg, pg, ng = pearson(x[good], y[good])
            txt = f'all days: r = {ra:+.2f}, p = {pa:.2f}\nratio > 2: r = {rg:+.2f}, p = {pg:.2f}'
            sl, ic = np.polyfit(x[ok], y[ok], 1)
            xx = np.linspace(np.nanmin(x), np.nanmax(x), 20); ax.plot(xx, ic + sl*xx, '-', color='0.3', lw=1.8, zorder=2)
            print(f'{fy} vs {fx}: r = {ra:+.2f} (p = {pa:.2f}, n = {na}); ratio > 2: r = {rg:+.2f} (p = {pg:.2f}, n = {ng})')
        ax.text(0.03, 0.96, txt, transform=ax.transAxes, va='top', fontsize=11, bbox=dict(fc='w', ec='0.4', boxstyle='square,pad=0.3'))
        ax.text(0.97, 0.96, f'({letters[k]})', transform=ax.transAxes, va='top', ha='right', fontsize=14, fontweight='bold'); k += 1
        ax.set_xlabel(xlab, fontsize=13); ax.grid(alpha=0.25); ax.tick_params(labelsize=12)
        ax.set_ylim(-60, 185)
        if j == 0: ax.set_ylabel(ylab, fontsize=13)
        if i == 0 and j == 0: ax.legend(loc='lower right', fontsize=11, framealpha=1)

# wind sectors and Kruskal-Wallis test (printed for the text)
sec = pd.cut(((day.wd + 45) % 360) // 90, bins=[-1, 0, 1, 2, 3], labels=['N', 'E', 'S', 'W'])
print(day.groupby(sec, observed=True).A.agg(['mean', 'count']).round(1))
print('Kruskal-Wallis across the sectors with more than one day:', stats.kruskal(*[day.A[sec == s] for s in ['N', 'E', 'S', 'W'] if (sec == s).sum() > 1]))
fig.tight_layout()
fig.savefig('FigS4_flux_vs_meteorology.png', dpi=250)
print('saved FigS4_flux_vs_meteorology.png')
