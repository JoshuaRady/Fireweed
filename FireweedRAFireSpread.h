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

//Fuel class weights from Albini ...
//May be converted into a class.
struct FuelWeights {
	std::vector<double> f_ij;
	std::vector<double> f_i;//Always length 2!
	std::vector<double> g_ij;
};

//The value of these two categories is forced to match the matching array indexes so they may be
//used to access arrays of form X_i:
//enum FuelCategory {Dead = 0, Live = 1};//Having some issues with casting this.
const int Dead = 0;//Temporary?????
const int Live = 1;//Temporary?????


double BulkDensity(double w_o, double fuelBedDepth);
double MeanBulkDensity(std::vector<double> w_o_ij, double fuelBedDepth);
double PackingRatio(double rho_b, double rho_p);
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth);
double OptimumPackingRatiofunction(double SAV, UnitsType units = ModelUnits);


////Utilities:

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<int> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4);

//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, FuelCategory liveDeadCat);
double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat);
bool FloatCompare(double val1, double val2, double precision = 0.0001);

void Warning(const char* message);
void Stop(const char* message);

#endif
