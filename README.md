**spectrumfitter** 

The module spectrum_fitter.py provides an obect-oriented framework for fitting dielectric spectra.

The fitting is performed with the constraint of the f-sum rule.

The longitudinal dielectric function can also be fit with minor modification of the code. Plots of the longitudinal fit and data (1 - 1/eps(omega)) are provided as an additional way to check the fit. 

An example of how to use spectrum_fitter.py is given in "fit_dielectric_spectrum_transverse_test.py" using the refractive index data compiled by Sigelstein, which can be found in the file "Siegelstein.RI"

There are a few additional python scripts included in this git: 
-- epsL_epsT.py - plots dielectric functions, etc
-- EpsLT.py - plots dielectric functions from index of refraction data at a variety of temperatures and attempts to locate peak maxima and park them with a line
-- kw.py - code for calculating k-dependent susceptibilities from k-dependent polarization-polarization correlation function data



The MIT License (MIT)

Copyright (c) 2015-2016 Daniel C. Elton 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

