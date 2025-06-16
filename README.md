<!-- ---------------------------------------------------------------------------
README.md
Joshua M. Rady
Woodwell Climate Research Center
2024

Note: This document is a GitHub ReadMe and uses the GitHub markdown dialect.
---------------------------------------------------------------------------- -->

# Fireweed
This project contains the **Fireweed Wildfire Code Library**.  The library provides implementations of wildfire behavior model equations from the wildfire modeling scientific literature.

## Overview and Purpose

The library currently implements equations to calculate the properties of surface fires.  At the center is the Rothermel fire spread model (Rothermel 1972) with the modifications of Albini (Albini 1976).  Additionally functions for managing and manipulating fire behavior fuel models and estimating fuel moisture are also provided (see [Library Modules](#library-modules) for more specifics).  This library is a work in progress.  Expect it to grow and change. 

The library was started as part of the [Fireweed Project](https://whrc.github.io/FireweedDocs/) to update the wildfire implementation in [DVM-DOS-TEM](https://github.com/uaf-arctic-eco-modeling/dvm-dos-tem/) (the _Dynamic Vegetation Model-Dynamic Organic Soil-Terrestrial Ecosystem Model_).  When we started work we found that public open source code for some commonly used mathematical models of Wildfire behavior was not readily available.  We felt that sharing our code might help others who want to learn about wildfire modeling, experiment with these models, or who want to couple these models to their ecosystem models (like we are doing).  With this in mind we have tried to make the library:

- **Flexible:** The code is high modular and has minimal dependancies.  It can be used for interactive experimentation on your laptop or can be compiled into a larger projects and deployed on any computing platform.

- **Easy to Understand:** We have [documented the code](#the-code) in great detail.  We showed our work.  We make the provenance of our calculations clear, citing all sources, and implemented our code as close to the original model equations as possible.  If you are new to wildfire modeling or are teaching it to others this should make it easier to follow and understand.  If you just want to run the code each function and it's inputs are clearly described.  _Use examples are in the works._


- **Reliable:** To the greatest extent possible we have verified that this code gives the same results as shown in the original papers documenting the respective models.

<!-- 
Easy to Use: This is a tricky one...

Open
 -->

## Library Contents <!-- Modules -->

### Surface Fire

The Rothermel and Albini fire spread rate model has been implemented for homogeneous and heterogenous fuels.  Related equations from several subsequent papers are also included.  The code can take use both the original US customary units and metric units.  In addition to the central spread rate calculation, the meny internal component calculations and be returned.  ...

**Add more on primary entry points...**

### Fuel Models

The Rothermel and Albini spread model use the concept of **fuel models** to characterize the surface fuel bed structure and properties.  The library includes a fuel model object model that simplifies the tasks of loading and manipulating fuel models.  In addition to making the spread rate equations easier to use fuel model objects provide unit handling and printing utilities.  Fuel models can be used to calculate live fuel curing for the dynamic fuel models added with "the 50" fuel models of Scott & Burgan 2005.

<!-- 
### Units
...
 -->

### Fuel Moisture

Equations for calculating dead fuel moisture are provided using the Rothermel 1983 / NWCG variant of the Fosberg model.  This method has the advantage of only requiring current day weather conditions but its estimates are made through simple lookup table based calculation.  These predictions will not be a reliable and calculations that consider past weather.  We plan to add a more robust method in the near future.

Dead fuel moisture can be calculated using an implementation of the GSI based system used in NFDRS 2016.

### Other

Utilities have been provided for some basic tasks.  These utilities are provided to prevent the need for external libraries and are fairly minimal.

## The Code

The Fireweed Wildfire Code Library is provided in two largely parallel versions, one in R, and one in C++.  R is particularly useful for experimentation since it can be used interactively.  We use it for one-off experiments, working out ideas, and prototyping.  C++ is significantly faster.  We use the C++ modules for integration with other scientific software.

### R Code
The R code is an informal library.  It is source-able but it not structured as a R package and there is no plan make it into one at this time.  Rather the priority was to keep the R code in a form that is easy to read and follow.  The code is written with minimal R dependancies (it is base R).  For maintainability efforts have been made to keep the code structure as parrallel to the C++ code as possible.

### C++ Code
The C++ code is structured so it can be compiled as a shared library or linked into an application.  The code includes R interface functions that allow some of the main code entry points in the shared library to be called from R code via the .C() function.

## Input Data
Input files containing the 53 standard fuel model parameters and the tables need for dead fuel moisture calculations are included.

#### Shared Library Compilation

The C++ code can be compiled as a shared library via g++, a part of [GCC](https://gcc.gnu.org).  We are successfully using GCC 13 on MacOS / Darwin, but other versions and platforms should work.  To compile do something like the following:

```
$ /usr/local/path/gcc/13.1.0/bin/g++-13 -c -o FireweedRAFireSpread.o FireweedRAFireSpread.cpp
$ /usr/local/path/gcc/13.1.0/bin/g++-13 -shared -o FireweedRAFireSpread.so FireweedRAFireSpread.o
```

## Documentation

The Fireweed Code Library is documented on this website and in the code itself.  Library level documentation can be found on this website.  More detailed documentation is embedded in the source code itself.  The C++ code functions and objects are documented with Doxygen.  (We plan to integrate this into the website tin the future.)    Source file headers contain additional details and references for each file.  Code within functions is further commented for those who what to understand how the code works.

Additionally the library will be documented in an future paper.

<!-- ## Rothermel and Albini Notation -->
### Notation

Rothermel & Albini is abbreviated as RA in code file names.

The Rothermel and Albini fire spread model has a lot of variables.  The notation used for them in both R and C++ is described [here](FireweedVariableNotation.md).  _Notation for some of the more recent components still needs to be added._

