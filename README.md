<!-- ---------------------------------------------------------------------------
README.md
Joshua M. Rady
Woodwell Climate Research Center
2024

Note: This document is a GitHub ReadMe and used the GitHub markdown dialect.
---------------------------------------------------------------------------- -->

# Fireweed
This project contains the **Fireweed wildfire code library**.  The library provides implementations of fire behavior model equations from the wildfire modeling scientific literature.  This is part of the [Fireweed Project](https://whrc.github.io/FireweedDocs/) to update the fire implementation in DVM-DOS-TEM.  The purpose of the library is to allow flexible code testing and exploration.  Expect this code to change rapidly.  The primary documentation is embedded as comments in the code itself.  These comments will be migrated to Doxygen in the future.

The library currently implements the Rothermel fire spread model (Rothermel 1972) with the modifications of Albini (Albini 1976) and related equations from several subsequent papers.  R was used for the in initial implementation and verification of the Rothermel Albini fire spread model prior to ported to C++.  Both versions are implemented as close to the original model equations as possible.

## R Code
The R code is an informal library.  It is source-able but it not structured as a R package and there is no plan make it into one at this time.  Rather the priority was to keep the R code in a form that is easy to read and follow.  The code is written with minimal R dependancies.  For maintainability efforts have been made to keep the code structure as parrallel to the C++ code as possible.

## C++ Code
The C++ code is structured so it can be compiled as a shared library or linked into an application.  The code includes R interface functions that allow the main code entry points in the shared library to be called from R code via the .C() function.

### Shared Library Compilation

The code can be compiled as a shared library via g++, a part of [GCC](https://gcc.gnu.org).  We are successfully using GCC 13 on MacOS / Darwin, but other versions and platforms should work.  To compile do something like the following:

```
$ /usr/local/path/gcc/13.1.0/bin/g++-13 -c -o FireweedRAFireSpread.o FireweedRAFireSpread.cpp
$ /usr/local/path/gcc/13.1.0/bin/g++-13 -shared -o FireweedRAFireSpread.so FireweedRAFireSpread.o
```

## Rothermel and Albini Notation

Rothermel & Albini is abbreviated as RA in code file names.

The Rothermel and Albini fire spread model has a lot of variables.  The notation used for them in both the R and C++ is described [here](FireweedVariableNotation.md).
