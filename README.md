JFV-strat
=========

Addition of a realistic stratospheric Newtonian cooling term in Held-Suarez-like General Circulation Models.


hs_forcing.f90
--------------

FORTRAN file which replaces the standard atmos_param/hs_forcing/hs_forcing.f90 file in the FMS file tree.
It adds the possibility to read in a Te and a tau profile from an input netCDF file, called INPUT/temp.nc and INPUT/tau.nc respectively.


Te_analytic.m
-------------

MATLAB script to compute analytic Te profiles and create the input file temp.nc for the FMS model.


tau_analytic.m
--------------

MATLAB script to compute analytic tau profiles and create the input file tau.nc for the FMS model.
