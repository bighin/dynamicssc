# Alignment dynamics of a molecule strongly-coupled to a many-body environment

This code allows one to study the alignment dynamics of a molecule strongly-coupled to a many-body environment.

The code allows for simulating a wide range of molecules, an arbitrary shape and strength of the alignment laser, as well as outputting a large number of time-dependent observables.

The theoretical model is a time-dependent version of the strong-coupling angulon variational Ansatz.

# Requirements

A modern C compiler (GCC or Clang) and the following libraries:

- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- NCurses

The following libraries/codes are included in the repository:

- [Cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature_(Multi-dimensional_integration)), for adaptive multidimensional integration.
- [inih](https://github.com/benhoyt/inih) for parsing .ini files.
- [libprogressbar](https://github.com/doches/progressbar) to display a nice progress bar.
- spline.c for spline interpolation

# Usage

Compile the code using the provided Makefile (`make`), eventually adapting it to point to the correct location of your C compiler and of your libraries. Prepare a .ini file with the physical details of the system you want to simulate, the `dsc.ini` contains all the possible options, along with comments. Finally run the code (`./dsc dsc.ini`).

