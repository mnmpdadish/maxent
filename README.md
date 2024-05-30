# BASIC README

OmegaMaxEnt is a tool for the analytic continuation of Matsubara data.
This is a modification of the original code OmegaMaxEnt by available 
at https://github.com/dbergeron1/OmegaMaxEnt


# License:

It is published under the [GNU Public License, version 3][license]
[license]: http://www.gnu.org/licenses/gpl.html
Note that the new part of the code here comes from another code
published under MIT license.


# To compile:

check makefile first (armadillo include header directory ad linking)
have armadillo installed
have lapack and blas (or mkl) installed
have gnuplot installed (available from command line).

Then type:

```shell
$ make
```

# To run:

```shell
$ cd testMaxEnt
$ ../maxEnt fictiveGreenI.dat
```

then check A.dat with gnuplot. Kind of basic for now. For now, there is no parameter file, and it is read directly into maxEnt_data.cpp. 

