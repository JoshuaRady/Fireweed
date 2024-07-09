/***************************************************************************************************
FireweedRAFireSpread.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/29/2024
Reference: Proj. 11 Exp. 14

	This file is part of the Fireweed fire code library.  It declares a set of functions
implementing the Rothermel fire spread model (Rothermel 1972) with the modifications of Albini
(Albini 1976) of the Fireweed fire code library.

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

//Fuel class weights from Albini 1976:
//This may be converted into a class.
struct FuelWeights {
	std::vector<double> f_ij;
	//f_i will always be length 2 so we could use an array or initialize the length with a constructor. 
	std::vector<double> f_i;//(2, 0);
	std::vector<double> g_ij;
};

//The values of the live dead categories are forced to match the matching array indexes so they may
//be used to access arrays of the form X_i (values are language specific):
//enum FuelCategory {Dead = 0, Live = 1};//Having some issues with casting this.
const int Dead = 0;
const int Live = 1;

//A data structure that holds the intermediate calculations of the spread rate calculation:
//This currently assumes heterogenous fuels.
//The x_i variables could probably be safely made fixed arrays of length 2.  However if we were to
//expand this to work with homogeneous fuels that could be a problem
//A print() function could be useful.
struct SpreadCalcs {//Name????? SpreadComponents
	UnitsType units;//The unit type for the values.
	double R;//Rate of spread in ft/min | m/min.
	FuelWeights weights;//Weights.
	//Heat source components:
	double GammaPrime;//Optimum reaction velocity (min^-1).
	std::vector <double> w_n_i;//Net fuel load for live/dead fuel categories (lb/ft^2 | kg/m^2).
	std::vector <double> h_i;//Heat content for live/dead fuel (Btu/lb | kJ/kg).
	std::vector <double> eta_M_i;//Moisture damping coefficient for live/dead fuel categories (unitless).
	std::vector <double> eta_s_i;//Mineral damping coefficient for live/dead fuel categories (unitless).
	double I_R;//Reaction intensity (Btu/ft^2/min | kJ/m^2/min).
	double xi;//Propagating flux ratio (dimensionless ratio).
	double phi_s;//Slope factor (dimensionless).
	double phi_w;//Wind factor (dimensionless).
	double heatSource;
	//Heat sink components:
	double rho_b_bar;//Mean bulk density
	//double epsilon;//Effective heating number.  Only calculated for homogeneous fuels!
	std::vector <double> Q_ig_ij;//Head of preignition
	double heatSink;
	//Other components that can be informative:
	double cSAV;//Fuel bed characteristic SAV
	double meanPackingRatio;
	double optimumPR;//Optimum packing ratio
	//RelativePR = (meanPackingRatio / optPackingRatio),#Overkill?
};

//Bulk Density:
double BulkDensity(double w_o, double fuelBedDepth);
double MeanBulkDensity(std::vector<double> w_o_ij, double fuelBedDepth);

//Packing Ratio:
double PackingRatio(double rho_b, double rho_p);
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth);
double OptimumPackingRatio(double cSAV, UnitsType units = ModelUnits);

//Weighting Factors:
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           std::vector<double> rho_p_ij, std::vector<int> liveDead,
                           UnitsType units = ModelUnits);
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           double rho_p, std::vector<int> liveDead, UnitsType units = ModelUnits);

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
                                                    std::vector <int> liveDead);
double LiveFuelMoistureOfExtinction(std::vector <double> M_f_ij, double M_x_1,
                                    std::vector <double> w_o_ij, std::vector <double> SAV_ij,
                                    std::vector <int> liveDead, UnitsType units = ModelUnits);
double MineralDampingCoefficient_Homo(double S_e);
std::vector <double> MineralDampingCoefficient_Het(std::vector <double> S_e_ij,
                                                   std::vector <double> f_ij,
                                                   std::vector <int> liveDead);

//Slope and Wind Factors:
double SlopeFactor(double packingRatio, double slopeSteepness);
double WindFactor(double cSAV, double packingRatio, double optPackingRatio, double U,
                  UnitsType units = ModelUnits);
double WindFactorC(double cSAV, UnitsType units = ModelUnits);
double WindFactorB(double cSAV, UnitsType units = ModelUnits);
double WindFactorE(double cSAV, UnitsType units = ModelUnits);
double WindFactorA(double cSAV, double packingRatio, double optPackingRatio,
                   UnitsType units = ModelUnits);
double WindLimit(double U, double I_R, UnitsType units = ModelUnits);

//Heat Source Components:
double OptimumReactionVelocity(double packingRatio, double cSAV, UnitsType units = ModelUnits);
std::vector <double> LiveDeadHeatContent(std::vector <double> h_ij, std::vector <double> f_ij,
                                         std::vector <int> liveDead);
double ReactionIntensity_Homo(double GammaPrime, double w_n, double h, double eta_M, double eta_s);
double ReactionIntensity_Het(double GammaPrime, std::vector <double> w_n_i,
                             std::vector <double> h_i, std::vector <double> eta_M_i,
                             std::vector <double> eta_s_i);
double PropagatingFluxRatio(double packingRatio, double cSAV, UnitsType units = ModelUnits);

//Heat Sink Components:
double EffectiveHeatingNumber(double SAV, UnitsType units = ModelUnits);
double HeatOfPreignition(double M_f, UnitsType units = ModelUnits);
std::vector <double> HeatOfPreignition(std::vector <double> M_f_ij, UnitsType units = ModelUnits);

//
double StdHeatContent(UnitsType units = ModelUnits);
double StdRho_p(UnitsType units = ModelUnits);

//Spread Rate Calculations:
double SpreadRateRothermelAlbini_Homo(double SAV, double w_o, double fuelBedDepth, double M_x,
                                      double M_f, double U, double slopeSteepness,
                                      double heatContent = StdHeatContent(),
                                      double S_T = 0.0555, double S_e = 0.01,
                                      double rho_p = StdRho_p(),
                                      bool useWindLimit = true,
                                      UnitsType units = US,
                                      bool debug = false);
double SpreadRateRothermelAlbini_Het(std::vector <double> SAV_ij,
                                     std::vector <double> w_o_ij,
                                     double fuelBedDepth,
                                     double M_x_1,
                                     std::vector <double> M_f_ij,
                                     double U, double slopeSteepness,//Both could default to 0.
                                     std::vector <double> h_ij,
                                     std::vector <double> S_T_ij,
                                     std::vector <double> S_e_ij,
                                     std::vector <double> rho_p_ij,
                                     std::vector <int> liveDead = {Dead, Dead, Dead, Live, Live},//Standard fuel model 5 classes.
                                     bool useWindLimit = false,
                                     UnitsType units = US,
                                     bool debug = false);
double SpreadRateRothermelAlbini_Het(std::vector <double> SAV_ij,
                                     std::vector <double> w_o_ij,
                                     double fuelBedDepth,
                                     double M_x_1,
                                     std::vector <double> M_f_ij,
                                     double U, double slopeSteepness,
                                     double h = StdHeatContent(),
                                     double S_T = 0.0555,
                                     double S_e = 0.01,
                                     double rho_p = StdRho_p(),
                                     std::vector <int> liveDead = {Dead, Dead, Dead, Live, Live},//Standard fuel model 5 classes.
                                     bool useWindLimit = false,
                                     UnitsType units = US,
                                     bool debug = false);

SpreadCalcs SpreadCalcsRothermelAlbini_Het(std::vector <double> SAV_ij,
                                           std::vector <double> w_o_ij,
                                           double fuelBedDepth,
                                           double M_x_1,
                                           std::vector <double> M_f_ij,
                                           double U, double slopeSteepness,//Both could default to 0.
                                           std::vector <double> h_ij,
                                           std::vector <double> S_T_ij,
                                           std::vector <double> S_e_ij,
                                           std::vector <double> rho_p_ij,
                                           std::vector <int> liveDead = {Dead, Dead, Dead, Live, Live},//Standard fuel model 5 classes.
                                           bool useWindLimit = false,
                                           UnitsType units = US,
                                           bool debug = false);

//Utilities:

//double StdHeatContent(UnitsType units = ModelUnits);
//double StdRho_p(UnitsType units = ModelUnits);

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<int> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4, std::vector<double> arg5, std::vector<double> arg6,
                 std::vector<double> arg7);
bool SameLengths(std::vector<double> dbl1, std::vector<double> dbl2, std::vector<double> dbl3,
                 std::vector<double> dbl4, std::vector<double> dbl5, std::vector<double> dbl6,
                 std::vector<double> dbl7, std::vector<int> int8);

bool InRange(double value, double low, double high);
bool InRange(std::vector<double> value, double low, double high);

bool ValidProportion(double value);
bool ValidProportion(std::vector<double> value);

//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, FuelCategory liveDeadCat);
double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat);
bool FloatCompare(double val1, double val2, double precision = 0.0001);

//Logging:
void LogMsg(const char* message);
void LogMsg(const char* message, double value);
void LogMsg(const char* message, std::vector<double> value);
void LogMsg(const char* message, std::vector<int> value);
void Warning(const char* message);
void Stop(const char* message);

//Related Fire Property Equations:
double EffectiveWindSpeed(double U, double phi_w, double phi_s, double meanPackingRatio,
                          double optPackingRatio, double cSAV, UnitsType units = ModelUnits);
double ResidenceTime(double cSAV, UnitsType units = ModelUnits);
double HeatPerUnitArea(double I_R, double t_r);
double ByramsFirelineIntensity(double H_A, double R);
double ByramsFlameLength(double I_B, UnitsType units = ModelUnits);

#endif
