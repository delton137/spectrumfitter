The module spectrum_fitter.py provides an obect-oriented framework for fitting dielectric spectra.

The fitting is performed using with the constraint of the f-sum rule heavily enforced.

The longitudinal dielectric function can also be fit with minor modification of the code. Plots of the longitudinal fit and data (1 - 1/eps(omega)) are provided an additional way to check the fit. 

An example of how to use spectrum_fitter.py is shown in "fit_dielectric_spectrum_transverse_test.py" using the refractive index data compiled in the file "Siegelstein.RI"

There are a few additional python scripts included in this git: 
-- epsL_epsT.py - plots dielectric functions, etc
-- EpsLT.py - plots dielectric functions from index of refraction data at a variety of temperatures and attempts to locate peak maxima and park them with a line
-- kw.py - code for calculating k-dependent susceptibilities from k-dependent polarization-polarization correlation function data

2015 Daniel C. Elton 