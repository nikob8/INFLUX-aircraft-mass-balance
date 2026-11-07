# INFLUX aircraft mass balance analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Code and per-transect results for

> Balashov, N., Barkley, Z. R., Davis, K. J., and Shepson, P. (2026). Evaluation of aircraft-based methane emissions estimates and estimation methodologies with six years of data from Indianapolis. *Journal of Geophysical Research: Atmospheres* (in revision).

The scripts compute Marion County (Indianapolis) CH4 emission rates from 37 aircraft mass balance flights of the Indianapolis Flux Experiment (INFLUX, 2011 to 2016) with four definitions of the background (Downwind Edge, Upwind, Adjusted ABL Upwind, Upwind Edge), the South Side Landfill emission rates, the sensitivity analyses of the Supporting Information, and all figures of the paper.

If you use this code or these data, please cite the paper and the Zenodo archive of this repository (the DOI badge above; replace XXXXXXX with the record number Zenodo assigns to the release).

## Contents

| Folder | What is in it |
|---|---|
| `data/` | `EPA_SSLF_GHGRP.csv` (EPA Greenhouse Gas Reporting Program estimates for the South Side Landfill, Figure 11); `Indy_Mass_Balance_Table_v5.csv` and `.xlsx` (Table S1: one row per downwind transect, emission rates of the four methods and of the landfill, backgrounds, ABL height, wind, cos θ); `diagnostic_ratios.csv` (diagnostic ratio of each flight day with the simulated wind and the quantities it is built from, Section 3.2); `Landfill_transects.csv` (landfill emission rate, plume width, enhancement and the plume end points of every transect with a quantified plume); `Trajectory_check_legs.csv` (per-transect results of the trajectory analysis, Text S2); `SSLF_plume_criteria.csv` (best-fit spike near the projected landfill position on every downwind transect, its classification by the automated selection and by the analysts, and the landfill emission rates with both edge definitions, Text S3) |
| `mass_balance/` | `MassBalanceMethodA.m` to `MassBalanceMethodD.m` (the four background methods, one script each; A = Downwind Edge, B = Upwind, C = Adjusted ABL Upwind, D = Upwind Edge); `Run_MethodsCD.m` (batch driver for C and D without figures); `Collect_Solutions.m` (gathers the per-day results into `All_Solutions.mat` and `All_Solutions.csv`); `MassBalanceLandfill.m` (mass balance of the South Side Landfill plume on each downwind transect, Section 2.4; writes the `Landfill/<yyyymmdd>_Loop<k>.mat` files; the plume end points were chosen by eye and are stored, so the script can rerun without interaction); `Collect_Landfill.m` (gathers the landfill results, including the stored plume end points, into `Landfill_transects.csv` and `Landfill_daily.csv`) |
| `analyses/` | `Wind_sensitivity_R2_5.m` (transect wind inside the integral of Equation 1, Section 2.3); `Trajectory_check_R2_7.m` (back trajectories with the measured and the WRF-Chem winds, Text S2 and Figure S6); `Check_PBLH_20121004.m` (WRF-Chem boundary layer height on 10/04/2012, Section 2.3.3); `Get_Flight_Wind_Directions.m` (transect-by-transect wind statistics for Table S1); `Diagnostic_ratio.m` (the WRF-Chem diagnostic ratio of each flight day, Section 3.2, computed from the released WRF subsets; writes `diagnostic_ratios.csv` and `Diag_ratio.mat`); `SSLF_plume_criteria.m` (objective characterization of the landfill plume selection, Text S3: Gaussian spike fit on every downwind transect, agreement with the analysts, sensitivity of the landfill estimate to the thresholds; writes `SSLF_plume_criteria.csv`, a summary figure and two inspection figures with the fit of every transect) |
| `wrf/` | `WRF_configuration.txt` (the complete WRF-Chem configuration of the simulations, recovered from the global attributes of the model output: grid, physics, dynamics, time step, eta levels, start dates of every run); `WRF_configuration_from_wrfout.m` (writes that file); `Extract_WRF_subset.m` (writes the subset of the model output used in the paper, one `.mat` file per flight day, see "WRF-Chem simulations" below) |
| `figures/` | one script per figure, see the table below; `alongTrackKm.m` is a helper used by `plot_Figure2_Edge.m` and `Flight_Atlas.m` |

### Figure scripts

| Figure | Script | Inputs |
|---|---|---|
| 1 | `plot_Figure1_seasonal.m` | `data/Indy_Mass_Balance_Table_v5.csv` |
| 2 | `plot_Figure2_Edge.m` (+ `alongTrackKm.m`) | aircraft data of 10/04/2012; flux labels from `data/Indy_Mass_Balance_Table_v5.csv` |
| 3, 5 | `plot_Figure3_5_Upwind.m` | aircraft data, Solution files of 10/04/2012 |
| 4 | `plot_Figure4_profiles.m` | aircraft data, `ElevationData.nc` |
| 7 | `Figure7_v4.m` | WRF-Chem output, colocated model and aircraft data of 03/01/2011 |
| 8 | `plot_Figure8_altitude.m` | `All_Solutions.mat` |
| 9 | `plot_Figure9_SSLF.m` | landfill results (`Landfill/<date>_Loop<k>.mat`, from `MassBalanceLandfill.m`), `All_Solutions.mat` |
| 10, S1 | `plot_FigS1_Fig10_trends.m` | `TempsAndTimes.mat`, `All_Solutions.mat`, landfill results |
| 11 | `plot_Figure11_SSLF_comparison.m` | `data/Landfill_daily.csv`, `data/EPA_SSLF_GHGRP.csv` (EPA GHGRP Equations HH-6 and HH-8 for the landfill, 2011 to 2023) |
| 12 | `Figure12_v2.m` | WRF-Chem maps and time series of 11/19/2014 |
| 13 | `plot_Figure13_diagnostic.m` | `Diag_ratio.mat` (diagnostic ratios, from `Diagnostic_ratio.m`), `All_Solutions.mat` |
| 14 | `plot_Figure14_bootstrap.m` | `Diag_ratio.mat`, `All_Solutions.mat` |
| S2 | `plot_FigS2_seasonal_flux.py` | `data/Indy_Mass_Balance_Table_v5.csv` |
| S3 | `plot_FigS3_SSLF_altitude.py` | `data/Indy_Mass_Balance_Table_v5.csv` |
| S4 | `plot_FigS4_flux_vs_meteorology.py` | `data/Indy_Mass_Balance_Table_v5.csv`, `data/diagnostic_ratios.csv` |
| S5 | `plot_FigS5_ratio_vs_wind.py` | same two csv files |
| S6 | `Trajectory_check_R2_7.m` | aircraft data, Solution files, WRF-Chem output |
| S7 | `plot_FigS7_SSLF_criteria.py` | `data/SSLF_plume_criteria.csv` |
| Document S2 (flight atlas) | `Flight_Atlas.m` | aircraft data, Solution files |

Figures S2, S3, S4, S5 and S7 (Python) and Figure 11 (MATLAB, no toolboxes) can be reproduced from the files in `data/` alone (copy the csv files into `figures/` and run, for example, `python3 plot_FigS4_flux_vs_meteorology.py`; `plot_FigS3_SSLF_altitude.py` also prints the comparison of the citywide emission rates on the transects with and without a quantified landfill plume of Section 3.1.3, and `plot_FigS7_SSLF_criteria.py` the agreement statistics and Table S2 of the Supporting Information). All other scripts need the input data described below.

## Input data

The scripts read the following files, which are not part of this repository:

* **Aircraft data**: one ASCII file per flight, `INFLUX_<yyyymmdd>`, with 1 Hz columns time (s UTC), latitude, longitude, altitude (m), CO2, CH4 (ppm), H2O, CO, pressure (hPa), temperature (K), wind speed (m/s), wind direction (deg). The INFLUX aircraft observations are archived as stated in the Open Research section of the paper.
* **Transect definitions**: `FlightList.mat` (list of flight days) and `LegAndProfileList/<yyyymmdd>.mat` with the index ranges of the upwind transect (`UpwindLeg`), the downwind transects (`LegList`) and the vertical profiles (`ProfileList`) in the aircraft file, and, for a few days, manually assigned boundary layer heights (`ManualABLHeights`) or a flagged CH4 series (`CH4Flagged`).
* **Terrain**: `ElevationData.nc` (used to convert altitude to height above ground).
* **Marion County corners**: hard-coded in the scripts (BottomRight, BottomLeft, TopLeft, TopRight).
* **Landfill results**: `Landfill/<yyyymmdd>_Loop<k>.mat`, one file per downwind transect, written by the landfill mass balance script (`Flag`, `Flux`, `Distance`, `Enhancement`, `Peak`; see `Collect_Landfill.m`).
* **WRF-Chem output** (`wrfout_d02_*` files, 3 km domain) for the diagnostic ratio, Figures 7 and 12, the boundary layer check and the trajectory analysis. The subset of the output that these analyses use is released separately, see the next section.
* **Third-party MATLAB functions**: `InterX.m` (curve intersections, MATLAB File Exchange 22441) is required by the Method scripts; `PlotINFLUX.m` and the Indianapolis shapefiles draw the roads on the maps and are optional.

Paths to these files are set in the settings block at the top of each script. Edit them before running.

## WRF-Chem simulations

The model fields of the paper come from WRF-Chem version 3.6.1 run as a tracer transport model (no chemistry): a 9 km outer domain over the Midwest and a 3 km nest (220 x 271 grid points, 51 levels) centred on Indianapolis, initial and boundary conditions from the North American Regional Reanalysis (NARR), Thompson microphysics, RRTM and Dudhia radiation, the MYNN 2.5 boundary layer scheme with its surface layer scheme, the Noah land surface model and no cumulus parameterization on either domain. Eight CH4 tracers carry the emissions of the regional source categories (oil and gas, coal, agriculture, landfills, other) and of Marion County (three tracers), each starting from a 300 ppb initial and boundary value. Each simulation started at 00 UTC on the day before the flight (24 h of spin-up) and ran for 48 h, except for four longer runs that cover consecutive flight days (11/09 to 11/22/2014, 11/30 to 12/04/2014, 03/03 to 03/08/2016, 03/16 to 03/19/2016). `wrf/WRF_configuration.txt` lists every namelist setting as recorded in the output files, for both domains, together with the start date and length of each simulation; the only setting that differs between runs is the time step (10 to 40 s on the 9 km domain, 3.3 to 13.3 s on the 3 km nest; 30 s and 10 s for most runs).

The full model output is too large to archive. `wrf/Extract_WRF_subset.m` writes, for each flight day, the fields the analyses use (earth-relative wind, potential temperature, height, pressure and boundary layer height of the lowest 25 levels on the 3 km grid within 160 km of the county centre, hourly from 8 h before the first transect to 1 h after the last, and the CH4 enhancement of each of the eight tracers, in ppb, during the flight hours and 16 to 20 UTC; the whole 3 km grid on 11/19/2014, the day of Figure 12) to `WRF_subset_<yyyymmdd>.mat`. The 37 files (about 100 MB each, 513 MB for 11/19/2014, 4.2 GB in total) are archived as a separate Zenodo data set (see the Open Research section of the paper); the variables are described inside each file (`info`).

## How to reproduce the emission rates

1. Run `MassBalanceMethodA.m` and `MassBalanceMethodB.m` (each loops over the 37 days and writes `<yyyymmdd>/Solution-MethodA.mat`, `-MethodB.mat` with the per-transect emission rates), then `Run_MethodsCD.m` for Methods C and D.
2. Run `Collect_Solutions.m` from the folder that holds the day folders. It writes `All_Solutions.mat` (structures `SolA` to `SolD`) and `All_Solutions.csv`, one row per transect with the emission rate of each method, the transect height, the boundary layer height, the mean wind and cos θ.
3. The daily means of `All_Solutions.csv` are the values of Tables 2 and 3 of the paper; Table S1 (`data/Indy_Mass_Balance_Table_v5.csv`) adds the landfill estimates and the backgrounds.
4. For the landfill, run `MassBalanceLandfill.m` (one `Landfill/<yyyymmdd>_Loop<k>.mat` per downwind transect; with `REUSE_BOUNDS = true` it takes the plume end points from `data/Landfill_transects.csv` instead of asking for them) and then `Collect_Landfill.m`, which applies cos θ of each transect and writes the per-transect and daily tables. `analyses/SSLF_plume_criteria.m` then characterizes the plume selection on every transect (Text S3).
5. The diagnostic ratio of each flight day (Section 3.2) is computed by `analyses/Diagnostic_ratio.m` from the WRF-Chem tracers along idealized downwind transects at about 300 m above ground, 17 to 20 UTC; it reads the released `WRF_subset_<yyyymmdd>.mat` files, so anyone can rerun it. `data/diagnostic_ratios.csv` holds the result for all 37 days.

The equations implemented by the scripts are Equations 1 to 4 of the paper. Each downwind transect gives one emission rate; the daily estimate is the mean of the transects of that day.

## Requirements

* MATLAB R2020b or later with the Statistics and Machine Learning Toolbox (`fitlm`) and the Mapping Toolbox (`worldmap`, `distance`, `shaperead`; used by the map figures and by the Method scripts for plotting).
* Python 3.9 or later with `numpy`, `pandas`, `scipy` and `matplotlib` for Figures S2, S3, S4, S5 and S7.

## License

MIT License (see `LICENSE`). The data in `data/` are released under the same terms; please cite the paper when using them.

## Contact

Nikolay Balashov (corresponding author; contact details in the paper)
