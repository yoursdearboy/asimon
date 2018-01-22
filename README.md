# asimon: Adaptive Two-Stage Designs for Cancer Trials by Yong Lin & Weichung J. Shih

This package calculates parameters for adaptive two-stage design (modified Simon's design [1]) proposed by Yong Lin & Weichung J. Shih [2].

Search algorithm and C code are straightforward with little optimisations, but still works after 2 to 10 minutes. That's my first attempt in writing R package with C, feel free to improve. Thank you!

## Requirements

GSL. On OS X `brew install gsl`, on Ubuntu: `sudo apt-get install libgsl0 libgsl0-dev`;

Also, On OS X you would like to use `gcc` instead of `clang` to enable parallelization. Install it with `brew install gcc` and put next lines to `~/.R/Makevars` (you may need to change version number):
```
CC=gcc-7
CXX=g++-7
CXX1X=g++-7
CXX11=g++-7
CXX14=g++-7
CXX17=g++-7

SHLIB_OPENMP_CFLAGS= -fopenmp
SHLIB_OPENMP_CXXFLAGS= -fopenmp
SHLIB_OPENMP_FCFLAGS= -fopenmp
SHLIB_OPENMP_FFLAGS= -fopenmp
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
# clone using git (or you can download archive and unzip it)
git clone https://github.com/yoursdearboy/asimon
cd asimon
# install
R CMD INSTALL .
```


## References

[1] Simon, Richard. "Optimal two-stage designs for phase II clinical trials." Controlled clinical trials 10, no. 1 (1989): 1-10.

[2] Lin, Yong, and Weichung J. Shih. "Adaptive Two‐Stage Designs for Single‐Arm Phase IIA Cancer Clinical Trials." Biometrics 60, no. 2 (2004): 482-490.
