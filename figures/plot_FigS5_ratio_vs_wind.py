#!/usr/bin/env python3
"""Figure S5: the diagnostic ratio of each flight day against the mean wind
direction and the mean wind speed (Section 3.2.3 of the paper).

Inputs (same folder):
    Indy_Mass_Balance_Table_v5.csv   per-transect results (Table S1): daily wind
    diagnostic_ratios.csv            diagnostic ratio of each flight day
Output: FigS5_ratio_vs_wind.png and the rank correlations printed.
Requires numpy, pandas, scipy, matplotlib.
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MINE_BEARING = 227.0        # bearing of the southern Illinois coal mines from Indianapolis (deg)
t = pd.read_csv('Indy_Mass_Balance_Table_v5.csv')
r = pd.read_csv('diagnostic_ratios.csv')
day = (t.groupby('Date').agg(ws=('Mean_wind_speed_ms', 'first'), wd=('Mean_wind_direction_deg', 'first')).reset_index()
         .merge(r[['Date', 'Diagnostic_ratio']], on='Date'))
ratio, ws, wd = day.Diagnostic_ratio.values, day.ws.values, day.wd.values
good = ratio > 2
ang = np.abs(((wd - MINE_BEARING) + 180) % 360 - 180)          # angular distance from the mine bearing
rho_a, p_a = stats.spearmanr(ang, ratio)
rho_a = 0.0 if abs(rho_a) < 0.005 else rho_a          # avoid printing -0.00
rho_s, p_s = stats.spearmanr(ws, ratio)
print(f'ratio vs angular distance from {MINE_BEARING:.0f} deg: Spearman rho = {rho_a:+.2f}, p = {p_a:.2f}')
print(f'ratio vs wind speed: Spearman rho = {rho_s:+.2f}, p = {p_s:.2f}')
in20 = ang <= 20
print(f'days within 20 deg of the bearing: {in20.sum()} ({(in20 & ~good).sum()} with ratio <= 2)')

fig, (a, b) = plt.subplots(1, 2, figsize=(16, 6.4), gridspec_kw=dict(width_ratios=[1.15, 1]))
a.axvspan(MINE_BEARING - 20, MINE_BEARING + 20, color='#e8c7a8', alpha=0.55, lw=0, label=f'within 20° of the coal-mine bearing ({MINE_BEARING:.0f}°)')
sc = a.scatter(wd, ratio, c=ws, s=120, cmap='viridis', vmin=2, vmax=12, edgecolor='k', linewidth=0.6, zorder=3)
a.axhline(2, ls='--', color='k', lw=1.4)
a.set_xlim(0, 360); a.set_xticks(np.arange(0, 361, 45)); a.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW', 'N'])
a.set_ylim(0, 8.6); a.set_xlabel('Mean wind direction (direction the wind blows from)', fontsize=14); a.set_ylabel('Diagnostic ratio', fontsize=14)
a.text(0.02, 0.97, f'ratio vs. angular distance from\nthe coal-mine bearing: ρ = {rho_a:.2f}, p = {p_a:.2f}', transform=a.transAxes, va='top', fontsize=12,
       bbox=dict(fc='w', ec='0.4', boxstyle='square,pad=0.3'))
a.legend(loc='lower left', fontsize=12, framealpha=1)
a.text(0.97, 0.04, '(a)', transform=a.transAxes, ha='right', fontsize=15, fontweight='bold')
cb = fig.colorbar(sc, ax=a, pad=0.02); cb.set_label('Mean wind speed (m s$^{-1}$)', fontsize=13)
a.grid(alpha=0.25); a.tick_params(labelsize=12)

b.scatter(ws[good], ratio[good], s=120, c='#1f4e96', edgecolor='k', linewidth=0.6, zorder=3, label=f'ratio > 2 ({good.sum()} days)')
b.scatter(ws[~good], ratio[~good], s=120, facecolor='none', edgecolor='#c0392b', linewidth=1.8, zorder=3, label=f'ratio $\\leq$ 2 ({(~good).sum()} days)')
b.axhline(2, ls='--', color='k', lw=1.4)
b.set_ylim(0, 8.6); b.set_xlabel('Mean wind speed (m s$^{-1}$)', fontsize=14); b.set_ylabel('Diagnostic ratio', fontsize=14)
b.text(0.02, 0.97, f'Spearman ρ = {rho_s:.2f}, p = {p_s:.2f}', transform=b.transAxes, va='top', fontsize=12, bbox=dict(fc='w', ec='0.4', boxstyle='square,pad=0.3'))
b.legend(loc='upper right', fontsize=12, framealpha=1)
b.text(0.97, 0.04, '(b)', transform=b.transAxes, ha='right', fontsize=15, fontweight='bold')
b.grid(alpha=0.25); b.tick_params(labelsize=12)
fig.tight_layout()
fig.savefig('FigS5_ratio_vs_wind.png', dpi=250)
print('saved FigS5_ratio_vs_wind.png')
