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
	//f_i will always be length 2 so we could use an array or initialize the length with a constructor. 
	std::vector<double> f_i;//(2, 0);
	std::vector<double> g_ij;
};

//The value of these two categories is forced to match the matching array indexes so they may be
//used to access arrays of form X_i:
//enum FuelCategory {Dead = 0, Live = 1};//Having some issues with casting this.
const int Dead = 0;//Temporary?????
const int Live = 1;//Temporary?????

//Bulk Density:
double BulkDensity(double w_o, double fuelBedDepth);
double MeanBulkDensity(std::vector<double> w_o_ij, double fuelBedDepth);

//Packing Ratio:
double PackingRatio(double rho_b, double rho_p);
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth);
double OptimumPackingRatio(double SAV, UnitsType units = ModelUnits);

//Weighting Factors:
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           std::vector<double> rho_p_ij, std::vector<int> liveDead,
                           UnitsType units = ModelUnits);
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           double rho_p, std::vector<int> liveDead, UnitsType units = ModelUnits)

//Fuel Bed Surface-area-to-volume Ratio:
double FuelBedSAV(std::vector<double> SAV_ij, std::vector<double> f_ij, std::vector<double> f_i,
                  std::vector<int> liveDead);

//Net Fuel Load:
double NetFuelLoad_Homo(double w_o, double S_T);
std::vector <double> NetFuelLoad_Het(std::vector <double> w_o_ij, std::vector <double> S_T_ij,
                                     std::vector <double> g_ij, std::vector <int> liveDead);

//Damping Coefficients:
double MoistureDampingCoefficient_Homo(double M_f, double M_x);
std::vector <double> MoistureDampingCoefficient_Het(std::vector <double> M_f_ij,
                                                    std::vector <double> M_x_i,
                                                    std::vector <double> f_ij,
                                                    std::vector <double> liveDead);
double LiveFuelMoistureOfExtinction(std::vector <double> M_f_ij, double M_x_1,
                                    std::vector <double> w_o_ij, std::vector <double> SAV_ij,
                                    std::vector <int> liveDead, UnitsType units = ModelUnits);
double MineralDampingCoefficient_Homo(double S_e);
std::vector <double> MineralDampingCoefficient_Het(std::vector <double> S_e_ij,
                                                   std::vector <double> f_ij,
                                                   std::vector <int> liveDead);

//Slope and Wind Factors:
double SlopeFactor(double packingRatio, double slopeSteepness);
double WindFactor(double SAV, double packingRatio, double optPackingRatio, double U,
                  UnitsType units = ModelUnits);
double WindFactorC(double SAV, UnitsType units = ModelUnits);
double WindFactorB(SAV, UnitsType units = ModelUnits);
double WindFactorE(SAV, UnitsType units = ModelUnits);
double WindFactorA(double SAV, double packingRatio, double optPackingRatio,
                   UnitsType units = ModelUnits);

//Heat Source Components:
double OptimumReactionVelocity(double packingRatio, double SAV, UnitsType units = ModelUnits);
double LiveDeadHeatContent(std::vector <double> h_ij, std::vector <double> f_ij,
                           std::vector <double> liveDead)
double ReactionIntensity_Homo (double GammaPrime, double w_n, double h, double eta_M, double eta_s);
double ReactionIntensity_Het(double GammaPrime, std::vector <double> w_n_i,
                             std::vector <double> h_i, std::vector <double> eta_M_i,
                             std::vector <double> eta_s_i)
double PropagatingFluxRatio(double packingRatio, double SAV, UnitsType units = ModelUnits);

//Heat Sink Components:
double EffectiveHeatingNumber(double SAV, UnitsType units = ModelUnits);
double HeatOfPreignition(double M_f, UnitsType units = ModelUnits);
std::vector <double> HeatOfPreignition(std::vector <double> M_f_ij, UnitsType units = ModelUnits);

//Utilities:

double StdHeatContent(UnitsType units = ModelUnits);
double StdRho_p(UnitsType units = ModelUnits);

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<int> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4);

bool InRange(double value, double low, double high);
bool InRange(std::vector<double> value, double low, double high);

bool ValidProportion(double value);
bool ValidProportion(std::vector<double> value);

//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, FuelCategory liveDeadCat);
double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat);
bool FloatCompare(double val1, double val2, double precision = 0.0001);

void Warning(const char* message);
void Stop(const char* message);

#endif
