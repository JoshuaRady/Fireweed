/***************************************************************************************************
FireweedRAFireSpread.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/29/2024
Reference: Proj. 11 Exp. 14

	This header file declares a set of functions implementing the Rothermel fire spread model
(Rothermel 1972) with the modifications of Albini (Albini 1976).  This is part header file is part
of the Fireweed fire code library.

***************************************************************************************************/
#ifndef FIREWEEDRAFIRESPREAD_H
#define FIREWEEDRAFIRESPREAD_H

#include <math.h>//Or cmath?????
#include <vector>
#include <numeric>
#include <iostream>

enum UnitsType {US, Metric};

void Warning(const char* message);

#endif
