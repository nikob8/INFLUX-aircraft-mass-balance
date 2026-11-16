# INFLUX aircraft mass balance analysis

[![DOI](https://zenodo.org/badge/1445764261.svg)](https://doi.org/10.5281/zenodo.24184833)

Code and per-transect results for

> Balashov, N. V., Barkley, Z. R., Davis, K. J., and Shepson, P. B. (2026). Evaluation of aircraft-based methane emissions estimates and estimation methodologies with six years of data from Indianapolis. *Journal of Geophysical Research: Atmospheres* (in revision).

The scripts compute Marion County (Indianapolis) CH4 emission rates from 37 aircraft mass balance flights of the Indianapolis Flux Experiment (INFLUX, 2011 to 2016) with four definitions of the background (Downwind Edge, Upwind, Adjusted ABL Upwind, Upwind Edge), the South Side Landfill emission rates, the sensitivity analyses of the Supporting Information, and all figures of the paper.

If you use this code or these data, please cite the paper and the Zenodo archive of this repository. Each release is archived at Zenodo with its own DOI, and the badge above always points to the latest version. Version 1.1 is the version cited in the paper. It differs from version 1.0 (https://doi.org/10.5281/zenodo.24184834) as follows: the folders `input/` and `results/` were added, with the auxiliary files the scripts read and the published results of the mass balance scripts, and the scripts find these folders by themselves; five files were renamed (`Indy_Mass_Balance_Table.csv` and `.xlsx` in `data/`, `plot_Figure12_case_study.m` in `figures/`, `Trajectory_check.m` and `Wind_sensitivity.m` in `analyses/`); the spreadsheet has a page layout for printing; the documentation and the comments in the scripts were revised. The computations are unchanged. `CITATION.cff` holds the citation metadata.

## Contents

| Folder | What is in it |
|---|---|
| `data/` | `EPA_SSLF_GHGRP.csv` (EPA Greenhouse Gas Reporting Program estimates for the South Side Landfill, Figure 11); `Barkley2026_SSLF_annual.csv` (annual landfill emission rates of Barkley et al. (2026), their Figure 5: Tower 10, 2017 to 2024, mean and full range of their sensitivity analysis in kg CH4/h, provided by Z. Barkley; Figure 11, please cite that paper); `Indy_Mass_Balance_Table.csv` and `.xlsx` (Table S1: one row per downwind transect, emission rates of the four methods and of the landfill, backgrounds, ABL height, wind, cos θ; the second worksheet of the `.xlsx` lists the daily estimates of the five flights of 2011 next to those of Cambaliza et al. (2015), Section 3.1.2); `diagnostic_ratios.csv` (diagnostic ratio of each flight day with the simulated wind and the quantities it is built from, Section 3.2); `Landfill_transects.csv` (landfill emission rate, plume width, enhancement and the plume end points of every transect with a quantified plume); `Landfill_daily.csv` (daily landfill emission rates, the mean of the quantified transects of each of the 31 days); `Trajectory_check_legs.csv` (per-transect results of the trajectory analysis, Text S2); `SSLF_plume_criteria.csv` (best-fit spike near the projected landfill position on every downwind transect, its classification by the automated selection and by the analysts, and the landfill emission rates with both edge definitions, Text S3); `emissions/` (the CH4 emissions of Figure 6 as text files with latitude, longitude and emission of each grid cell in kg CH4 per hour: the five regional sectors on the 9 km and 3 km WRF-Chem grids, from the EPA gridded inventory of Maasakkers et al. (2016), and the 1 km Marion County emissions of Lamb et al. (2016) for natural gas, landfill and wastewater, as received; the simulations and Figure 6 use the landfill file multiplied by 1.4 and the wastewater file divided by 3.5; please cite those papers) |
| `input/` | Files the scripts read that are not results of the paper. `MassBalance/`: `FlightList.mat` (the 37 flight days), `LegAndProfileList/<yyyymmdd>.mat` (transects and profiles of each day) and `TempsAndTimes.mat` (takeoff temperature and day number of each flight day, Figures 10 and S1). `helpers/`: `ElevationData.nc` (terrain) and `PlotINFLUX.m` (draws the interstates on the maps). `shapefiles/`: interstates of Indianapolis and outline of Marion County. `WRF_20141119/`: WRF-Chem CH4 enhancement and wind on the 3 km grid at 16, 17 and 18 UTC on 11/19/2014 with the aircraft observations of each hour, and the model along the flight track (Figure 12). `colocated/`: WRF-Chem sampled along the flight track of 03/01/2011 (Figure 7). `aircraft/` and `WRF_subsets/` are the places for the aircraft files and the WRF-Chem subsets, which are archived elsewhere (see "Input data") |
| `results/` | Results of the mass balance scripts as used in the paper: `All_Solutions.mat` (per-transect results of the four methods for all flight days, structures `SolA` to `SolD`), the day folder `20121004/` with `Solution-MethodA.mat` to `Solution-MethodD.mat` (the example day of Figures 2, 3 and 5), `Landfill/<yyyymmdd>_Loop<k>.mat` (landfill mass balance of each downwind transect, written by `MassBalanceLandfill.m`) and `Diag_ratio.mat` (diagnostic ratios, written by `analyses/Diagnostic_ratio.m`). The scripts in `mass_balance/` write into this folder |
| `mass_balance/` | `MassBalanceMethodA.m` to `MassBalanceMethodD.m` (the four background methods, one script each; A = Downwind Edge, B = Upwind, C = Adjusted ABL Upwind, D = Upwind Edge); `Run_MethodsCD.m` (batch driver for C and D without figures); `Collect_Solutions.m` (gathers the per-day results into `results/All_Solutions.mat` and `All_Solutions.csv`); `MassBalanceLandfill.m` (mass balance of the South Side Landfill plume on each downwind transect, Section 2.4; writes the `results/Landfill/<yyyymmdd>_Loop<k>.mat` files; the plume end points were chosen by eye and are stored, so the script can rerun without interaction); `Collect_Landfill.m` (gathers the landfill results, including the stored plume end points, into `Landfill_transects.csv` and `Landfill_daily.csv`) |
| `analyses/` | `Wind_sensitivity.m` (transect wind inside the integral of Equation 1, Section 2.3); `Trajectory_check.m` (back trajectories with the measured and the WRF-Chem winds, Text S2 and Figure S6); `Check_PBLH_20121004.m` (WRF-Chem boundary layer height on 10/04/2012, Section 2.3.3 and Text S4); `Get_Flight_Wind_Directions.m` (transect-by-transect wind statistics for Table S1); `Diagnostic_ratio.m` (the WRF-Chem diagnostic ratio of each flight day, Section 3.2, computed from the released WRF subsets; writes `diagnostic_ratios.csv` and `Diag_ratio.mat`); `SSLF_plume_criteria.m` (objective characterization of the landfill plume selection, Text S3: Gaussian spike fit on every downwind transect, agreement with the analysts, sensitivity of the landfill estimate to the thresholds; writes `SSLF_plume_criteria.csv`, a summary figure and two inspection figures with the fit of every transect); `Table_S1_statistics.m` (statistics quoted in the paper that follow directly from Table S1: time of day of the downwind transects in UTC and local time, Section 2.2; order in which the transect altitudes were flown, Section 3.1.1; landfill plume and landfill emission rates by wind direction, Section 3.1.3 and Text S3) |
| `wrf/` | `WRF_configuration.txt` (the complete WRF-Chem configuration of the simulations, recovered from the global attributes of the model output: grid, physics, dynamics, time step, eta levels, start dates of every run); `WRF_configuration_from_wrfout.m` (writes that file); `Extract_WRF_subset.m` (writes the subset of the model output used in the paper, one `.mat` file per flight day, see "WRF-Chem simulations" below); `WRF_subsets_README.txt` (description of those files and of their variables, with reading examples for MATLAB and Python; the same text accompanies the Zenodo data set) |
| `figures/` | one script per figure, see the table below; `alongTrackKm.m` is a helper used by `plot_Figure2_Edge.m` and `Flight_Atlas.m` |

### Table S1 and Data Sets S1 to S4 of the paper

| In the paper | File |
|---|---|
| Table S1 | `data/Indy_Mass_Balance_Table.xlsx` |
| Data Set S1 | `data/Indy_Mass_Balance_Table.csv` |
| Data Set S2 | `data/diagnostic_ratios.csv` |
| Data Set S3 | `data/Trajectory_check_legs.csv` |
| Data Set S4 | `data/SSLF_plume_criteria.csv` |

### Figure scripts

| Figure | Script | Inputs |
|---|---|---|
| 1 | `plot_Figure1_seasonal.m` | `data/Indy_Mass_Balance_Table.csv` |
| 2 | `plot_Figure2_Edge.m` (+ `alongTrackKm.m`) | aircraft data of 10/04/2012; flux labels from `data/Indy_Mass_Balance_Table.csv` |
| 3, 5 | `plot_Figure3_5_Upwind.m` | aircraft data of 10/04/2012, `results/20121004/Solution-MethodB.mat` to `-MethodD.mat` |
| 4 | `plot_Figure4_profiles.m` | aircraft data of 10/04/2012, `input/helpers/ElevationData.nc` |
| 6 | `plot_Figure6_emissions.m` | `data/emissions/` (emission files of the 9 km and 3 km domains and of Marion County), `input/shapefiles/` |
| 7 | `plot_Figure7_WRF.m` | WRF-Chem output and aircraft data of 03/01/2011, `input/colocated/Colocated_20110301.mat` (the model sampled along the flight track), `input/shapefiles/` |
| 8 | `plot_Figure8_altitude.m` | `results/All_Solutions.mat` |
| 9 | `plot_Figure9_SSLF.m` | `results/Landfill/` (landfill results of each transect), `results/All_Solutions.mat` |
| 10, S1 | `plot_FigS1_Fig10_trends.m` | `input/MassBalance/TempsAndTimes.mat`, `results/All_Solutions.mat`, `results/Landfill/` (Figure 10) |
| 11 | `plot_Figure11_SSLF_comparison.m` | `data/Landfill_daily.csv`, `data/EPA_SSLF_GHGRP.csv` (EPA GHGRP Equations HH-6 and HH-8 for the landfill, 2011 to 2023), `data/Barkley2026_SSLF_annual.csv` (annual tower estimates, 2017 to 2024; the figure shows 2017 to 2019, `YEAR_MAX` in the script) |
| 12 | `plot_Figure12_case_study.m` | `input/WRF_20141119/` (WRF-Chem maps and time series of 11/19/2014), `input/shapefiles/` |
| 13 | `plot_Figure13_diagnostic.m` | `results/Diag_ratio.mat` (diagnostic ratios, from `Diagnostic_ratio.m`), `results/All_Solutions.mat` |
| 14 | `plot_Figure14_bootstrap.m` | `results/Diag_ratio.mat`, `results/All_Solutions.mat` |
| S2 | `plot_FigS2_seasonal_flux.py` | `data/Indy_Mass_Balance_Table.csv` |
| S3 | `plot_FigS3_SSLF_altitude.py` | `data/Indy_Mass_Balance_Table.csv` |
| S4 | `plot_FigS4_flux_vs_meteorology.py` | `data/Indy_Mass_Balance_Table.csv`, `data/diagnostic_ratios.csv` |
| S5 | `plot_FigS5_ratio_vs_wind.py` | same two csv files |
| S6 | `Trajectory_check.m` | aircraft data, Solution files of all flight days, WRF-Chem output |
| S7 | `plot_FigS7_SSLF_criteria.py` | `data/SSLF_plume_criteria.csv` |
| Document S1 (flight atlas) | `Flight_Atlas.m` | aircraft data, `results/All_Solutions.mat` |

The following scripts need only the files of the repository: Figures 1, 6, 8, 9, 10, 11, 12, 13, 14 and S1 (MATLAB), Figures S2, S3, S4, S5 and S7 (Python, for example `python3 plot_FigS4_flux_vs_meteorology.py`), `mass_balance/Collect_Landfill.m` and `analyses/Table_S1_statistics.m` (the latter in base MATLAB or GNU Octave). `plot_FigS3_SSLF_altitude.py` also prints the comparison of the citywide emission rates on the transects with and without a quantified landfill plume of Section 3.1.3, and `plot_FigS7_SSLF_criteria.py` the agreement statistics and Table S2 of the Supporting Information. All other scripts need the input data described below.

## Input data

Run each script from its own folder. The scripts address the folders `data/`, `input/` and `results/` relative to that folder, so there are no paths to edit, with one exception (the WRF-Chem output, below). Tables and figures are written into the current folder; the scripts in `mass_balance/` write into `results/`.

Three inputs are not in the repository. Each has its place in `input/`:

* **Aircraft data**, in `input/aircraft/`: one ASCII file per flight, `INFLUX_<yyyymmdd>`, with 1 Hz columns time (s UTC), latitude, longitude, altitude (m), CO2, CH4 (ppm), H2O, CO, pressure (hPa), temperature (K), wind speed (m/s), wind direction (deg). The INFLUX aircraft observations are archived as stated in the Open Research section of the paper.
* **`InterX.m`**, in `input/helpers/`: curve intersections by NS, MATLAB Central File Exchange 22441 (https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections). The mass balance scripts and the scripts that repeat their geometry call it.
* **WRF-Chem fields**: the subset of the model output that the paper uses is a separate Zenodo data set (next section); `analyses/Diagnostic_ratio.m` reads its files from `input/WRF_subsets/`. Five scripts read the complete model output (`wrfout_d02_*` files), which is too large to archive: `figures/plot_Figure7_WRF.m`, `analyses/Trajectory_check.m`, `analyses/Check_PBLH_20121004.m`, `wrf/Extract_WRF_subset.m` and `wrf/WRF_configuration_from_wrfout.m`. In these, the folder `WRF_ROOT` in the settings block has to be set.

The repository holds the other inputs:

* **Transect definitions**: `input/MassBalance/FlightList.mat` (list of the 37 flight days) and `input/MassBalance/LegAndProfileList/<yyyymmdd>.mat` with the index ranges of the upwind transect (`UpwindLeg`, 32 days), the downwind transects (`LegList`) and the vertical profiles (`ProfileList`) in the aircraft file, and, where needed, manually assigned boundary layer heights (`ManualABLHeights`, 18 days) or a flagged CH4 series (`CH4Flagged`, one day).
* **Terrain**: `input/helpers/ElevationData.nc`, global elevation on a 0.25 degree grid (written by T. Mitchell, JISAO, University of Washington, from the TerrainBase 5-minute data set), used to convert altitude to height above ground.
* **Maps**: `input/helpers/PlotINFLUX.m` (Z. Barkley) draws the interstates from `input/shapefiles/Indianapolis_Interstate_Mileposts.shp`, a layer of the open GIS data of the City of Indianapolis and Marion County; `Marion_County_Boundary.shp` is the outline of the county. The Marion County corners used in the calculations are hard-coded in the scripts (BottomRight, BottomLeft, TopLeft, TopRight).
* **Results of the mass balance**: `results/All_Solutions.mat` (all days, all four methods) and `results/20121004/Solution-Method?.mat`. The day folders of the other flight days are written by the Method scripts; `analyses/SSLF_plume_criteria.m` and `analyses/Trajectory_check.m` read them.
* **Landfill results**: `results/Landfill/<yyyymmdd>_Loop<k>.mat` (146 files, written by `mass_balance/MassBalanceLandfill.m`), one file per downwind transect k: `Flag` (larger than 0: plume quantified; -1: no landfill plume identified) and, for the quantified plumes, `Flux`, `Distance`, `Enhancement`, `Peak`, `Background`, the CH4 series of the plume (`CH4LF`) and its end points (`Start`, `End`); see `Collect_Landfill.m`. Figures 9 and 10 read them. `data/Landfill_transects.csv` lists the transects with a quantified plume.
* **Model fields of the two case studies**: `input/colocated/Colocated_20110301.mat` (Figure 7) and `input/WRF_20141119/` (Figure 12).

## WRF-Chem simulations

The model fields of the paper come from WRF-Chem version 3.6.1 run as a tracer transport model (no chemistry): a 9 km outer domain over the Midwest and a 3 km nest (220 x 271 grid points, 51 levels) centred on Indianapolis, initial and boundary conditions from the North American Regional Reanalysis (NARR), Thompson microphysics, RRTM and Dudhia radiation, the MYNN 2.5 boundary layer scheme with its surface layer scheme, the Noah land surface model and the Kain–Fritsch cumulus scheme on the 9 km domain only. Eight CH4 tracers carry the emissions of the regional source categories (1 oil and gas, 2 coal, 3 agriculture, 4 landfills, 5 other) and of Marion County (6 natural gas distribution, 7 landfill, 8 wastewater treatment). The Marion County emissions follow Lamb et al. (2016), with the landfill emission multiplied by 1.4 and the wastewater emission divided by 3.5 to bring them closer to the estimates of that study: 56, 28 and 2 mol/s in total. Each tracer starts from, and is held at the domain boundaries to, a constant value (300 in the raw model output, in which the emissions are scaled by 10^6), so the analyses use the enhancement in ppb, (tracer - 300)/10^6. Each simulation started at 00 UTC on the day before the flight (24 h of spin-up) and ran for 48 h, except for four longer runs that cover consecutive flight days (11/09 to 11/22/2014, 11/30 to 12/04/2014, 03/03 to 03/08/2016, 03/16 to 03/19/2016). `wrf/WRF_configuration.txt` lists every namelist setting as recorded in the output files, for both domains, together with the start date and length of each simulation; the only setting that differs between runs is the time step (10 to 40 s on the 9 km domain, 3.3 to 13.3 s on the 3 km nest; 30 s and 10 s for most runs).

The full model output is too large to archive. `wrf/Extract_WRF_subset.m` writes, for each flight day, the fields the analyses use (earth-relative wind, potential temperature, height, pressure and boundary layer height of the lowest 25 levels on the 3 km grid within 160 km of the county centre, hourly from 8 h before the first transect to 1 h after the last, and the CH4 enhancement of each of the eight tracers, in ppb, during the flight hours and 16 to 20 UTC; the whole 3 km grid on 11/19/2014, the day of Figure 12) to `WRF_subset_<yyyymmdd>.mat`. The 37 files (about 100 MB each, about 500 MB for 11/19/2014, 4.0 GB in total) are archived as a separate Zenodo data set (https://doi.org/10.5281/zenodo.24184894); the variables are described inside each file (`info`). `analyses/Diagnostic_ratio.m` reads the files from `input/WRF_subsets/`.

## How to reproduce the emission rates

1. Put the aircraft files into `input/aircraft/` and `InterX.m` into `input/helpers/`. In `mass_balance/`, run `MassBalanceMethodA.m` and `MassBalanceMethodB.m` (each loops over the 37 days and writes `results/<yyyymmdd>/Solution-MethodA.mat`, `-MethodB.mat` with the per-transect emission rates), then `Run_MethodsCD.m` for Methods C and D.
2. Run `Collect_Solutions.m`. It gathers the day folders into `results/All_Solutions.mat` (structures `SolA` to `SolD`) and `results/All_Solutions.csv`, one row per transect with the emission rate of each method, the transect height, the boundary layer height, the mean wind and cos θ. It stops if fewer than the 37 day folders are present, so that the complete file of the repository is not replaced by a partial one.
3. The daily means of `All_Solutions.csv` are the values of Tables 2 and 3 of the paper; Table S1 (`data/Indy_Mass_Balance_Table.csv`) adds the landfill estimates and the backgrounds.
4. For the landfill, run `MassBalanceLandfill.m` (one `results/Landfill/<yyyymmdd>_Loop<k>.mat` per downwind transect; it takes the plume end points from `data/Landfill_transects.csv`, or asks for them if `REUSE_BOUNDS = false`) and then `Collect_Landfill.m`, which applies cos θ of each transect and writes the per-transect and daily tables. `analyses/SSLF_plume_criteria.m` then characterizes the plume selection on every transect (Text S3).
5. The diagnostic ratio of each flight day (Section 3.2) is computed by `analyses/Diagnostic_ratio.m` from the WRF-Chem tracers along idealized downwind transects at about 300 m above ground, 17 to 20 UTC; it reads the released `WRF_subset_<yyyymmdd>.mat` files from `input/WRF_subsets/`, so anyone can rerun it. `data/diagnostic_ratios.csv` holds the result for all 37 days.

The equations implemented by the scripts are Equations 1 to 4 of the paper. Each downwind transect gives one emission rate; the daily estimate is the mean of the transects of that day.

## Requirements

* MATLAB R2020b or later with the Statistics and Machine Learning Toolbox (`fitlm`) and the Mapping Toolbox (`worldmap`, `distance`, `shaperead`; used by the map figures and by the Method scripts for plotting).
* Python 3.9 or later with `numpy`, `pandas`, `scipy` and `matplotlib` for Figures S2, S3, S4, S5 and S7.

## License

MIT License (see `LICENSE`). The data in `data/` are released under the same terms; please cite the paper when using them. `input/helpers/ElevationData.nc` and the shapefiles in `input/shapefiles/` come from the sources named under "Input data"; they are included for convenience and remain subject to the terms of their sources.

## Contact

Nikolay Balashov (corresponding author; contact details in the paper)
