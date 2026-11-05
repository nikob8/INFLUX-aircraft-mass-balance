#!/usr/bin/env python3
"""Figure S2: daily citywide CH4 emission rates (Downwind Edge approach) by
meteorological season (Section 3.1.2 of the paper).

Input (same folder): Indy_Mass_Balance_Table_v5.csv (Table S1, one row per
    downwind transect). The daily value is the mean of the per-transect Downwind
    Edge emission rates of that day.
Output: FigS2_seasonal_flux.png; the seasonal counts and means, the overall mean
    and the Kruskal-Wallis p-value are printed.
Box = interquartile range, line = median, whiskers = full range, points = flight days.
Requires numpy, pandas, scipy, matplotlib.
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t = pd.read_csv('Indy_Mass_Balance_Table_v5.csv')
day = t.groupby('Date').Flux_DownwindEdge_mols.mean()
dates = pd.to_datetime(day.index)
season = np.full(len(day), 4)                      # 4 = SON
season[dates.month.isin([12, 1, 2])] = 1
season[dates.month.isin([3, 4, 5])] = 2
season[dates.month.isin([6, 7, 8])] = 3
names = ['Winter (DJF)', 'Spring (MAM)', 'Summer (JJA)', 'Fall (SON)']
groups = [day.values[season == s] for s in (1, 2, 3, 4)]
h, p_kw = stats.kruskal(*groups)
print('n per season:', [len(g) for g in groups], '| means:', [round(g.mean()) for g in groups],
      f'| overall mean {day.mean():.1f} (SD {day.std(ddof=1):.1f}) | Kruskal-Wallis H = {h:.2f}, p = {p_kw:.2f}')

rng = np.random.default_rng(1)
box_col, pt_col, w = '#9fb8d9', '#1f4e9c', 0.55
fig, ax = plt.subplots(figsize=(7.4, 5.4))
for s, v in enumerate(groups, start=1):
    q1, q2, q3 = np.percentile(v, [25, 50, 75])
    ax.add_patch(plt.Rectangle((s - w / 2, q1), w, q3 - q1, facecolor=box_col, edgecolor='k', lw=1, zorder=2))
    ax.plot([s - w / 2, s + w / 2], [q2, q2], 'k-', lw=1.8, zorder=3)
    ax.plot([s, s], [q3, v.max()], 'k-', lw=1, zorder=2); ax.plot([s, s], [v.min(), q1], 'k-', lw=1, zorder=2)
    ax.plot([s - w / 4, s + w / 4], [v.max(), v.max()], 'k-', lw=1); ax.plot([s - w / 4, s + w / 4], [v.min(), v.min()], 'k-', lw=1)
    ax.scatter(s + (rng.random(len(v)) - 0.5) * 0.24, v, s=26, c=pt_col, edgecolor='k', linewidth=0.3, zorder=4)
    ax.text(s, -52, f'n = {len(v)}, mean = {v.mean():.0f}', ha='center', va='top', fontsize=10)
hm = ax.axhline(day.mean(), color='r', ls='--', lw=1.3, label=f'All flights: mean = {day.mean():.0f} mol s$^{{-1}}$')
ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(names, fontsize=11)
ax.set_xlim(0.4, 4.6); ax.set_ylim(-70, 180)
ax.set_ylabel('Citywide CH$_4$ emission rate, Downwind Edge approach (mol s$^{-1}$)', fontsize=11)
ax.set_title(f'Daily emission rates by season ({len(day)} flight days)', fontsize=12)
ax.legend(loc='upper left', fontsize=10, framealpha=1)
ax.tick_params(labelsize=10)
fig.tight_layout()
fig.savefig('FigS2_seasonal_flux.png', dpi=300)
print('saved FigS2_seasonal_flux.png')
