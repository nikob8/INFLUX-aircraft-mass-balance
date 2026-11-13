WRF-Chem fields for the INFLUX aircraft mass balance flights, 2011 to 2016
==========================================================================

Companion data set of
  Balashov, N. V., Barkley, Z. R., Davis, K. J., and Shepson, P. B. (2026). Evaluation of
  aircraft-based methane emissions estimates and estimation methodologies with six
  years of data from Indianapolis. Journal of Geophysical Research: Atmospheres.
Analysis code: https://github.com/nikob8/INFLUX-aircraft-mass-balance
(version 1.0 archived at https://doi.org/10.5281/zenodo.24184834)

Contents
--------
WRF_subset_<yyyymmdd>.mat   one file per flight day (37 files), MATLAB v7.3 (HDF5)
                            format, readable with MATLAB (load), Python (h5py) or
                            any HDF5 library
WRF_configuration.txt       the complete model configuration (every namelist
                            setting of both domains as recorded in the output
                            files, the eta levels, the start date and length of
                            each simulation)

The model
---------
WRF-Chem version 3.6.1 run as a tracer transport model (no chemistry): 9 km outer
domain over the Midwest and a 3 km nest (220 x 271 grid points, 51 levels) centred
on Indianapolis; initial and boundary conditions from the North American Regional
Reanalysis; Thompson microphysics, RRTM longwave and Dudhia shortwave radiation,
MYNN 2.5 boundary layer scheme with the MYNN surface layer, Noah land surface
model, no cumulus parameterization on the 3 km domain. Each simulation started at
00 UTC on the day before the flight (24 h of spin-up). Eight CH4 tracers carry the
emissions of the regional source categories of the EPA gridded inventory
(Maasakkers et al., 2016): 1 oil and gas, 2 coal, 3 agriculture, 4 landfills,
5 other, and of Marion County (Lamb et al., 2016): 6 natural gas distribution,
7 landfill, 8 wastewater treatment (the landfill emission multiplied by 1.4 and the
wastewater emission divided by 3.5, toward the estimates of that study: 56, 28 and
2 mol/s in total). Every tracer
starts from, and is held at the domain boundaries to, a constant value (300 in the
raw model output, in which the emissions are scaled by 10^6); the TRACER variable
below therefore holds the enhancement in ppb, (tracer - 300)/10^6.

What each file contains
-----------------------
The fields are those used by the analyses of the paper: the 3 km grid within
160 km of the Marion County centre (the whole 3 km grid on 11/19/2014), the
lowest 25 model levels (about 3 km above ground), every hour from 8 h before the
first transect of the flight to 1 h after the last one (always including 16 to
20 UTC). Arrays are single precision, dimensions (x, y, level, hour) with
x = west to east and y = south to north on the model grid. The diagnostic ratio
of the paper (Section 3.2) can be recomputed from these files with
Diagnostic_ratio.m of the code repository.

grid.XLAT, grid.XLONG   latitude and longitude of the mass points (degrees)
grid.HGT                terrain height (m above sea level)
grid.COSALPHA, SINALPHA map rotation already applied to U and V
grid.i_range, j_range   indices of the sub-domain in the full 3 km grid
grid.ZNU, grid.ZNW      eta values of the mass and staggered levels kept
grid.z_agl_mean         mean height above ground of each level (m)
grid.global_attributes  all global attributes of the output file (namelist)
grid.run_directory      name of the simulation the fields come from
grid.flight_hours       UTC hours of the first and last aircraft transect
hours                   UTC hours of the stored fields
hasTracer               true for the hours with tracers: from one hour before
                        the first transect to one hour after the last one, and
                        16 to 20 UTC (the window of the diagnostic ratio)
tracerHours             the UTC hours for which the tracers are stored
U, V                    earth-relative wind on mass points (m/s)
TH                      potential temperature (K)
Z                       height of the mass levels above sea level (m); subtract
                        grid.HGT for height above ground
P                       pressure (Pa)
PBLH                    boundary layer height of the MYNN scheme (m), (x, y, hour)
TRACER                  CH4 enhancement of each tracer (ppb), dimensions
                        (x, y, level, tracer hour, tracer); regional background
                        enhancement = sum of tracers 1 to 5, city enhancement =
                        sum of tracers 6 to 8 (the 300 offset and the 1e6
                        emission scaling of the raw WRF-Chem tracers are
                        already removed)
info                    the same descriptions as text, inside the file

Example (MATLAB)
----------------
  S = load('WRF_subset_20121004.mat');
  city = sum(S.TRACER(:,:,:,:,6:8), 5);              % city CH4 enhancement (ppb)
  h = find(S.tracerHours == 17);                     % 17 UTC
  pcolor(S.grid.XLONG, S.grid.XLAT, city(:,:,8,h)); shading flat   % level 8, about 300 m

Example (Python)
----------------
  import h5py, numpy as np
  f = h5py.File('WRF_subset_20121004.mat', 'r')
  U = np.array(f['U'])            # note: h5py returns the dimensions reversed
                                  # (hour, level, y, x); TRACER is (tracer, hour, level, y, x)

License: Creative Commons Attribution 4.0 (CC BY 4.0). Please cite the paper and
this data set when using these files.

References
  Lamb, B. K., et al. (2016). Direct and indirect measurements and modeling of
    methane emissions in Indianapolis, Indiana. Environmental Science and
    Technology, 50, 8910 to 8917.
  Maasakkers, J. D., et al. (2016). Gridded national inventory of U.S. methane
    emissions. Environmental Science and Technology, 50, 13123 to 13133.
