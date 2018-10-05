[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://pypi.python.org/pypi/spectrumfitter)
[![PyPI version](https://badge.fury.io/py/spectrumfitter.svg)](https://badge.fury.io/py/spectrumfitter)

# spectrumfitter

The **spectrumfitter** package is an obect-oriented framework for fitting dielectric spectra.

The following lineshapes are supported (bold text refers to the class name): 
* **Debye** - Debye dielectric relxation 
* **DHO** - Damped harmonic oscillator 
* **Gaussian()** - Gaussian function for imaginary part and corresponding real part. 
* **Constant()** - constant function for real part, 0 imaginary part. 
* **ColeCole()** - Cole-Cole dielectric relaxation 
* **VanVleck()** - Van Vleck & Weisskopf lineshape [*Rev. Mod. Phys.*, **17**:227 236, (1945)]
* **BrendelDHO** - Brendel DHO [Ref: *J. Appl. Phys*. **71**, 1 (1992) ]
* **PowerLawDebye()** - Debye lineshape with truncated power law wing. [*J. Phys. Chem. B* **109** (12), 6031 (2005)] 
* **StrExp** - Stretched exponential relaxation lineshape (a non-analytic function).


Besides basic fitting it has several advanced features, such as: 
* Use the *f*-sum rule as a constraint 
* Use the generalized LST (gLST) relation as a constraint 
* Fit and plot the longitudinal dielectric function

Example codes are shown in the `examples` directory. 



