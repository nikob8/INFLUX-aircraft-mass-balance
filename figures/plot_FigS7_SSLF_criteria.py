#!/usr/bin/env python3
"""Figure S7: objective characterization of the SSLF plume selection (Text S3).

Input (same folder): SSLF_plume_criteria.csv, written by analyses/SSLF_plume_criteria.m,
    one row per downwind transect with the best-fit Gaussian spike near the projected
    position of the SSLF (amplitude A_ppb above the polynomial background, signal-to-noise
    ratio SNR = A / rms residual of the fit, full width at half maximum fwhm_km, measured
    peak above the edge values P_ppb), the SSLF emission rate with the fitted edges
    (Flux_alg_mols, edges at mu +/- 2.5 sigma) and with the edges chosen by the analysts
    (Flux_manual_mols), and the classification of the analysts (det_manual).
Output: FigS7_SSLF_criteria.png; the agreement statistics, the daily means and the
    threshold sensitivity table (Table S2) are printed.
The detection is recomputed here from the fitted quantities with the thresholds below
(those of Text S3), so the script does not depend on the det_alg column of the file.
Requires numpy, pandas, matplotlib.
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PMIN, SNRMIN, FWMAX = 10.0, 7.0, 10.0          # amplitude (ppb), signal-to-noise ratio, width (km)
RELAXED = (5.0, 4.0, 15.0)                       # permissive set quoted in Text S3

d = pd.read_csv('SSLF_plume_criteria.csv')
man = d.det_manual.astype(bool).values
A, snr, fw, P = d.A_ppb.values, d.SNR.values, d.fwhm_km.values, d.P_ppb.values
fa, fm = d.Flux_alg_mols.values, d.Flux_manual_mols.values

def select(pmin, snrmin, fwmax):
    return (A >= pmin) & (snr >= snrmin) & (fw <= fwmax) & (P >= pmin)

def daily(mask, flux):
    g = pd.Series(flux[mask]).groupby(d.Date.values[mask]).mean()
    return g.mean(), g.median(), 2 * g.std(ddof=1) / np.sqrt(len(g)), len(g)

def report(pmin, snrmin, fwmax, label):
    det = select(pmin, snrmin, fwmax)
    both, neither = (det & man).sum(), (~det & ~man).sum()
    mo, oo = (man & ~det).sum(), (det & ~man).sum()
    m, md, se, n = daily(det, fa)
    print(f'{label}: A >= {pmin:g} ppb, SNR >= {snrmin:g}, FWHM <= {fwmax:g} km: {det.sum()} detected; '
          f'agreement {both + neither} of {len(d)} ({100 * (both + neither) / len(d):.0f}%): both {both}, neither {neither}, '
          f'analysts only {mo}, automated only {oo}; daily mean {m:.1f} +/- {se:.1f} mol/s ({n} days, median {md:.1f})')
    return det

det = report(PMIN, SNRMIN, FWMAX, 'Text S3 thresholds')
report(*RELAXED, 'relaxed thresholds')
m, md, se, n = daily(man, fm)
print(f'selection of the analysts: {man.sum()} plumes; daily mean {m:.1f} +/- {se:.1f} mol/s ({n} days, median {md:.1f})')
both = det & man
r = np.corrcoef(fa[both], fm[both])[0, 1]
ratio = np.median(fa[both] / fm[both])
print(f'transects selected by both: n = {both.sum()}, r = {r:.2f}, median ratio fitted/manual {ratio:.2f}, '
      f'median difference {np.median(fa[both] - fm[both]):+.1f} mol/s')
for name, mask in (('quantified by the analysts', man), ('not quantified', ~man)):
    print(f'{name}: A median {np.median(A[mask]):.1f} ppb (10th pct {np.percentile(A[mask], 10):.1f}), '
          f'SNR median {np.median(snr[mask]):.1f} (10th pct {np.percentile(snr[mask], 10):.1f}), '
          f'FWHM median {np.median(fw[mask]):.1f} km (90th pct {np.percentile(fw[mask], 90):.1f})')

# fitted spike on the transects without a quantified plume, relative to the quantified transects of the same day
from scipy import stats
rat = []
for day, s_ in d.groupby('Date'):
    m_ = s_.det_manual.astype(bool)
    if m_.any() and (~m_).any():
        rat.append(s_.A_ppb[~m_].mean() / s_.A_ppb[m_].mean())
rat = np.array(rat)
print(f'same-day amplitude ratio, not quantified / quantified: median {np.median(rat):.2f} over {len(rat)} days, '
      f'smaller on {np.sum(rat < 1)} days, Wilcoxon p = {stats.wilcoxon(np.log(rat)).pvalue:.3f}')

print('\nTable S2: n detected / agreement / daily mean +/- 2 SE (days), amplitude threshold 10 ppb')
for s in (4, 5, 7, 8):
    row = []
    for f in (8, 10, 15):
        dd = select(10.0, s, f)
        m, md, se, n = daily(dd, fa)
        row.append(f'{dd.sum():3d} / {100 * ((dd & man).sum() + (~dd & ~man).sum()) / len(d):2.0f}% / {m:4.1f} +/- {se:.1f} ({n})')
    print(f'SNR >= {s}: ' + ' | '.join(row))

# ----------------------------------------------------------------------------- figure
BLUE, RED = '#1f4e96', '#c0392b'
fig, (a, b) = plt.subplots(1, 2, figsize=(16, 6.4))
a.scatter(snr[man], A[man], s=110, c=BLUE, edgecolor='k', linewidth=0.6, zorder=3,
          label=f'plume quantified by the analysts ({man.sum()})')
a.scatter(snr[~man], A[~man], s=110, facecolor='none', edgecolor=RED, linewidth=1.8, zorder=3,
          label=f'no plume quantified ({(~man).sum()})')
wide = fw > FWMAX
a.scatter(snr[wide], A[wide], s=60, marker='x', c='k', linewidth=1.4, zorder=4,
          label=f'width above {FWMAX:g} km ({wide.sum()})')
a.axvline(SNRMIN, ls='--', color='k', lw=1.4, label=f'thresholds: {PMIN:g} ppb, {SNRMIN:g} times the residual')
a.axhline(PMIN, ls='--', color='k', lw=1.4)
a.set_xscale('log'); a.set_yscale('log')
a.set_xlim(1, 50); a.set_ylim(2, 1200)
a.set_xticks([1, 2, 5, 10, 20, 50]); a.set_xticklabels(['1', '2', '5', '10', '20', '50'])
a.set_yticks([2, 5, 10, 20, 50, 100, 200, 500, 1000]); a.set_yticklabels(['2', '5', '10', '20', '50', '100', '200', '500', '1000'])
a.set_xlabel('Amplitude of the fitted spike / rms residual of the fit', fontsize=14)
a.set_ylabel('Amplitude of the fitted spike (ppb)', fontsize=14)
a.legend(loc='upper left', fontsize=11.5, framealpha=1)
a.text(0.97, 0.04, '(a)', transform=a.transAxes, ha='right', fontsize=15, fontweight='bold')
a.grid(alpha=0.25, which='both'); a.tick_params(labelsize=12)

lim = max(fa[both].max(), fm[both].max()) * 1.08
b.plot([0, lim], [0, lim], 'k--', lw=1.4, zorder=2)
b.scatter(fm[both], fa[both], s=110, c=BLUE, edgecolor='k', linewidth=0.6, zorder=3)
b.set_xlim(0, lim); b.set_ylim(0, lim); b.set_aspect('equal')
b.set_xlabel('SSLF emission rate, edges chosen by the analysts (mol s$^{-1}$)', fontsize=14)
b.set_ylabel('SSLF emission rate, fitted edges (mol s$^{-1}$)', fontsize=14)
b.text(0.03, 0.97, f'transects selected by both: n = {both.sum()}\nr = {r:.2f}, median ratio {ratio:.2f}',
       transform=b.transAxes, va='top', fontsize=12, bbox=dict(fc='w', ec='0.4', boxstyle='square,pad=0.3'))
b.text(0.97, 0.04, '(b)', transform=b.transAxes, ha='right', fontsize=15, fontweight='bold')
b.grid(alpha=0.25); b.tick_params(labelsize=12)
fig.tight_layout()
fig.savefig('FigS7_SSLF_criteria.png', dpi=250)
print('saved FigS7_SSLF_criteria.png')
