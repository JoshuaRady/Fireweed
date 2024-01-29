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
#include <iostream>//std::cout

//Globals:------------------------------------------------------------------------------------------
enum UnitsType {US, Metric};

//Specify the units to use.  The default is United States customary units.
//This should not be set directly.  Use SetModelUnits().
//private 
UnitsType ModelUnits = US;


double PackingRatio(double rho_b, double rho_p);
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth);
double OptimumPackingRatiofunction(double SAV, UnitsType units = ModelUnits);

void Warning(const char* message);
void Stop(const char* message);

#endif
