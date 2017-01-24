JFV-strat [![DOI](https://zenodo.org/badge/19286204.svg)](https://zenodo.org/badge/latestdoi/19286204)
=========

Realistic stratospheric Newtonian cooling setup in Held-Suarez-like General Circulation Models, as described in 

[M Jucker, S Fueglistaler, and GK Vallis (2013): Maintenance of the Stratospheric Structure in an Idealized General Circulation Model, J. Atmos. Sci. 70, 3341, DOI: 10.1175/JAS-D-12-0305.1](http://dx.doi.org/10.1175/JAS-D-12-0305.1),

and the analytical approximation described in 

[M Jucker, S Fueglistaler, and GK Vallis (2014): Stratospheric sudden warmings in an idealized GCM, J. Geophys. Res., DOI: 10.1002/2014JD022170] (http://dx.doi.org/10.1002/2014JD022170)

This model is based on the Geophyiscal Fluid Dynamics Laboratory's (GFDL) AM2 code. For online documentation, see
[GFDL's AM2 documentation website](http://data1.gfdl.noaa.gov/~arl/pubrel/m/am2/doc/), and for references, see

[The GFDL Global Atmospheric Model Development Team, 2004: The new GFDL global atmosphere and land model AM2-LM2: Evaluation with prescribed SST simulations. Journal of Climate, 17(24), 4641-4673](http://dx.doi.org/10.1175/JCLI-3223.1)

[Delworth, T. L., et al., 2006: GFDL's CM2 Global Coupled Climate Models. Part I: Formulation and simulation characteristics. Journal of Climate, 19(5), 643-674](http://dx.doi.org/10.1175/JCLI3629.1).

[Lin, S-J., 2004: A "vertically Lagrangian" finite-volume dynamical core for global models. Monthly Weather Review, 132(10), 2293-2307](http://dx.doi.org/10.1175/1520-0493(2004)132%3C2293:AVLFDC%3E2.0.CO;2).

In particular, this code uses the spectral HSt42 experiment of the Riga release (December 2010). If you already have 
that experiment from the Flexible Modeling System (FMS), you probably only need to replace your hs_forcing.F90 with
fms_riga_hs_jucker/src/atmos_param/hs_forcing/hs_forcing.f90 and you should be good to go.


Te_analytic.m
-------------

MATLAB script to compute analytic Te profiles and create the input file temp.nc for the FMS model.


tau_analytic.m
--------------

MATLAB script to compute analytic tau profiles and create the input file tau.nc for the FMS model.

heldsuarez.m / heldsuarez_tau.m
-------------------------------

MATLAB scripts to compute Held and Suarez (1994) Te and tau profiles. These are used for the troposphere.

polvanikushner.m
----------------

MATLAB script to compute Polvani and Kushner (2002) stratosphere. Useful to compare to more basic model setup.

writetempin.m
-------------

MATLAB script to create the needed netCDF file for the GCM. It computes the traditional Held-Suarez setup below a given pressure p_hsin, and uses an input profile Tin above p_bdin. For pressure surfaces between p_hsin and p_bdin, it interpolates linearly between the two profiles. The input Tin can come from Te_analytic, or any other source (e.g. temp_monthly_L10_full.nc).

writetauin.m
------------

MATLAB script to create the needed netCDF file for the GCM. Same functioning as writetempin.m. tauin can come from tau_analytic.m, or any other source (e.g. tau_monthly_L10_full.nc)


temp / tau_monthly_L10_full.nc
----------------------------

Example files for non-analytic Newtonian cooling profiles. These files come from work described in
M Jucker, S Fueglistaler, and GK Vallis (2013): Maintenance of the Stratospheric Structure in an Idealized General Circulation Model, J. Atmos. Sci. 70, 3341, DOI: 10.1175/JAS-D-12-0305.1 


fms_riga_hs_jucker/src/atmos_param/hs_forcing/hs_forcing.f90
--------------

Specific to the dry dynamical core of the GCM AM2 of the Geophysical Fluid Dynamics Laboratory (GFDL). The FORTRAN file  replaces the standard atmos_param/hs_forcing/hs_forcing.f90 file in the FMS file tree.

There are two ways to get the stratospheric forcing working:

1) It adds the possibility to read in a Te and a tau profile from an input netCDF file, called INPUT/temp.nc and INPUT/tau.nc respectively (and easily produced with the MATLAB scripts provided in this package).

To make it work, add temp.nc and/or tau.nc to the folder INPUT in the run directory, and set equilibrium_t_option and/or equilibrium_tau_option to 'from_file' in the namelist hs_forcing_nml.

2) The analytical profile defined in the [2014 paper](http://onlinelibrary.wiley.com/doi/10.1002/2014JD022170/abstract) is hard coded into hs_forcing. Set equilibrium_t_option and/or equilibrium_tau_option to 'JFV', and set the parameters as described in the paper and the code. 
