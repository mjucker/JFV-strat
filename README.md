JFV-strat
=========

Realistic stratospheric Newtonian cooling setup in Held-Suarez-like General Circulation Models.


Te_analytic.m
-------------

MATLAB script to compute analytic Te profiles and create the input file temp.nc for the FMS model.


tau_analytic.m
--------------

MATLAB script to compute analytic tau profiles and create the input file tau.nc for the FMS model.


writetempin.m
-------------

MATLAB script to create the needed netCDF file for the GCM. It computes the traditional Held-Suarez setup below a given pressure p_hsin, and uses an input profile Tin above p_bdin. For pressure surfaces between p_hsin and p_bdin, it interpolates linearly between the two profiles. The input Tin can come from Te_analytic, or any other source (e.g. temp_monthly_L10_full.nc).

writetauin.m
------------

MATLAB script to create the needed netCDF file for the GCM. Same functioning as writetempin.m. tauin can come from tau_analytic.m, or any other source (e.g. tau_monthly_L10_full.nc)


temp/tau_monthly_L10_full.nc
----------------------------

Example files for non-analytic Newtonian cooling profiles. These files come from work described in
M Jucker, S Fueglistaler, and GK Vallis (2013): Maintenance of the Stratospheric Structure in an Idealized General Circulation Model, J. Atmos. Sci. 70, 3341, DOI: 10.1175/JAS-D-12-0305.1 


hs_forcing.f90
--------------

Specific to the dry dynamical core of the GCM AM2 of the Geophysical Fluid Dynamics Laboratory (GFDL). The FORTRAN file  replaces the standard atmos_param/hs_forcing/hs_forcing.f90 file in the FMS file tree.
It adds the possibility to read in a Te and a tau profile from an input netCDF file, called INPUT/temp.nc and INPUT/tau.nc respectively (and easily produced with the MATLAB scripts provided in this package).

To make it work, add temp.nc and/or tau.nc to the folder INPUT in the run directory, and set equilibrium_t_option and/or equilibrium_tau_option to 'from_file' in the namelist hs_forcing_nml.
