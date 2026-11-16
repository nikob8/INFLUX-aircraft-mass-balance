Aircraft data (not part of this repository)
===========================================

The scripts expect one plain ASCII file per flight day in this folder, named
INFLUX_<yyyymmdd> (no extension, no header), with 1 Hz records and the columns

   1 time (s UTC)       5 CO2 (ppm)          9 pressure (hPa)
   2 latitude           6 CH4 (ppm)         10 temperature (K)
   3 longitude          7 H2O               11 wind speed (m/s)
   4 altitude (m)       8 CO                12 wind direction (deg)

The INFLUX aircraft observations are archived as stated in the Open Research
section of the paper.  The 37 flight days are listed in
input/MassBalance/FlightList.mat and in data/Indy_Mass_Balance_Table.csv.
