# asimon: Adaptive Two-Stage Designs for Cancer Trials by Yong Lin & Weichung J. Shih

This package calculates parameters for adaptive two-stage design (modified Simon's design [1]) proposed by Yong Lin & Weichung J. Shih [2].

Search algorithm and C code are straightforward with little optimisations, but still work for simple designs (one stage design n ~ 55 or less).

For optimized version see branch [optimized](https://github.com/yoursdearboy/asimon/tree/optimized) and even better [cached](https://github.com/yoursdearboy/asimon/tree/cached). Feel free to improve. Thank you!

## Requirements

GSL. On OS X `brew install gsl`, on Ubuntu: `sudo apt-get install libgsl0 libgsl0-dev`;

Also, On OS X you would like to use `gcc` instead of `clang` to enable parallelization. Install it with `brew install gcc` and put next lines to `~/.R/Makevars` (you may need to change -9 to your gcc version number):
```
CC=gcc-9
CXX=g++-9
CXX1X=g++-9
CXX11=g++-9
CXX14=g++-9
CXX17=g++-9

SHLIB_OPENMP_CFLAGS= -fopenmp
SHLIB_OPENMP_CXXFLAGS= -fopenmp
SHLIB_OPENMP_FCFLAGS= -fopenmp
SHLIB_OPENMP_FFLAGS= -fopenmp
```

If you want to use clang anyway or other compiler, put this to `~/.R/Makevars`.
```
SHLIB_OPENMP_CFLAGS=
SHLIB_OPENMP_CXXFLAGS=
SHLIB_OPENMP_FCFLAGS=
SHLIB_OPENMP_FFLAGS=
```

## Installation

Using `devtools`:

```
install.packages('devtools') # if you don't have it
library(devtools)
devtools::install_github('yoursdearboy/asimon')
```

Or manually:

```
# clone git branch you want (or download archive and unzip it)
git clone --single-branch --branch optimized https://github.com/yoursdearboy/asimon.git
cd asimon
# install dependencies and package
R -e "install.packages('clinfun')"
R CMD INSTALL .
```

## References

[1] Simon, Richard. "Optimal two-stage designs for phase II clinical trials." Controlled clinical trials 10, no. 1 (1989): 1-10.

[2] Lin, Yong, and Weichung J. Shih. "Adaptive Two‐Stage Designs for Single‐Arm Phase IIA Cancer Clinical Trials." Biometrics 60, no. 2 (2004): 482-490.
