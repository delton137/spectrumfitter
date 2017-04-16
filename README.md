[![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()

**spectrumfitter** 

The ` obect-oriented framework for fitting dielectric spectra.

The fitting is performed with the constraint of the f-sum rule.

The longitudinal dielectric function can also be fit with minor modification of the code. Plots of the longitudinal fit and data (1 - 1/eps(omega)) are provided as an additional way to check the fit. 

An example of how to use spectrum_fitter.py is given in "fit_dielectric_spectrum_transverse_test.py" using the refractive index data compiled by Sigelstein, which can be found in the file "Siegelstein.RI"

There are a few additional python scripts included in this git: 
-- epsL_epsT.py - plots dielectric functions, etc
-- EpsLT.py - plots dielectric functions from index of refraction data at a variety of temperatures and attempts to locate peak maxima and park them with a line
-- kw.py - code for calculating k-dependent susceptibilities from k-dependent polarization-polarization correlation function data


