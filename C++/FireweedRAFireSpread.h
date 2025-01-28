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

//#include "FireweedUnits.h"
#include "FireweedFuelModels.h"

//Globals:------------------------------------------------------------------------------------------

/** Specify the default units to use.  The default is United States customary units.
 * This should not be set directly.
 * 
 * Note: This is only provided in the header to allow it to be used for default values in the
 * functions that follow.  This functionality could be achieved via function overloading or some
 * other method.  Currently ModelUnits is not used by any outside code and should be private.
 */
extern UnitsType ModelUnits;

/** @struct SpreadCalcs
 * @brief A data structure that holds fuel class weights from Albini 1976.
 *
 * @note This may be converted into a class?????
 */
struct FuelWeights {
	std::vector<double> f_ij;
	//f_i will always be length 2 so we could use an array or initialize the length with a constructor. 
	std::vector<double> f_i;//(2, 0);
	std::vector<double> g_ij;
};

/** See FireweedFuelModels.h for explanation of live dead categories */

/** @struct SpreadCalcs
 * @brief A data structure that holds the intermediate calculations of the spread rate calculation.
 *
 * The homogeneous and heterogenous forms of the spread calculations have some differences in their
 * intermediate calculations.  Members that have different interpretations based on the fuel model
 * type are indicated by trailing x subscripts (i.e. _x).  See the comments below.
 */
struct SpreadCalcs {//Name????? SpreadComponents
	UnitsType units;//The unit type for the values.
	bool homogeneous;//True if homogeneous fuel model output, false if heterogeneous fuel model output.
	
	double R;//Rate of spread in ft/min | m/min.
	FuelWeights weights;//Weights.  Only calculated for heterogeneous fuels!
	
	//Heat source components:
	double GammaPrime;//Optimum reaction velocity (min^-1).
	//Net fuel load for live/dead fuel categories (lb/ft^2 | kg/m^2): 
	std::vector <double> w_n_x;//homogeneous: w_n, heterogeneous: w_n_i
	//Heat content for live/dead fuel (Btu/lb | kJ/kg):
	std::vector <double> h_x;//homogeneous: h, heterogeneous: h_i
	//Moisture damping coefficient for live/dead fuel categories (unitless):
	std::vector <double> eta_M_x;//homogeneous: eta_M, heterogeneous: eta_M_i
	//Mineral damping coefficient for live/dead fuel categories (unitless):
	std::vector <double> eta_s_x;//homogeneous: eta_s, heterogeneous: eta_s_i
	double I_R;//Reaction intensity (Btu/ft^2/min | kJ/m^2/min).
	double xi;//Propagating flux ratio (dimensionless ratio).
	double phi_s;//Slope factor (dimensionless).
	double phi_w;//Wind factor (dimensionless).
	double heatSource;
	
	//Heat sink components:
	double rho_b_x;//(Mean) bulk density: homogeneous: rho_b, heterogeneous: rho_b_bar
	double epsilon;//Effective heating number.  Only calculated for homogeneous fuels! (-1 otherwise)
	std::vector <double> Q_ig_x;//Heat of preignition: homogeneous: Q_ig, heterogeneous: Q_ig_ij
	double heatSink;
	
	//Other components that can be informative:
	double cSAV;//Fuel bed characteristic SAV (= SAV for homogeneous)
	double packingRatio;//(Mean) packing ratio		xPackingRatio?????
	double optimumPR;//Optimum packing ratio
	
	std::ostream& Print(std::ostream& output) const;//Print the struct contents.
};

std::ostream& operator<<(std::ostream& output, const SpreadCalcs& fm);//Overloaded stream print operator.

//Functions:----------------------------------------------------------------------------------------

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

//These are utilities that need to be declared here to be available for the following functions:
double StdHeatContent(UnitsType units = ModelUnits);
double StdRho_p(UnitsType units = ModelUnits);

//Spread Rate Calculations:
double SpreadRateRothermelAlbini_Homo(double SAV, double w_o, double fuelBedDepth, double M_x,
                                      double M_f, double U, double slopeSteepness,
                                      double heatContent = StdHeatContent(),
                                      double S_T = 0.0555, double S_e = 0.01,
                                      double rho_p = StdRho_p(),
                                      bool useWindLimit = true,//True for all homo versions!!!!!
                                      UnitsType units = US,
                                      bool debug = false);
double SpreadRateRothermelAlbini_Homo(FuelModel fuelModel,
                                      double M_f, double U, double slopeSteepness,
                                      bool useWindLimit = true,
                                      bool debug = false);

SpreadCalcs SpreadCalcsRothermelAlbini_Homo(double SAV, double w_o, double fuelBedDepth, double M_x,
                                            double M_f, double U, double slopeSteepness,
                                            double heatContent = StdHeatContent(),
                                            double S_T = 0.0555, double S_e = 0.01,
                                            double rho_p = StdRho_p(),
                                            bool useWindLimit = true,
                                            UnitsType units = US,
                                            bool debug = false);
SpreadCalcs SpreadCalcsRothermelAlbini_Homo(FuelModel fuelModel,
                                            double M_f, double U, double slopeSteepness,
                                            bool useWindLimit = true,
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
double SpreadRateRothermelAlbini_Het(FuelModel fuelModel,
                                     double U, double slopeSteepness,
                                     std::vector <double> M_f_ij = {},
                                     bool useWindLimit = false,
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
SpreadCalcs SpreadCalcsRothermelAlbini_Het(FuelModel fuelModel,
                                           double U, double slopeSteepness,
                                           std::vector <double> M_f_ij = {},
                                           bool useWindLimit = false,
                                           bool debug = false);

//Utilities:
double SumByFuelCat(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat);

//Related Fire Property Equations:
double EffectiveWindSpeed(double U, double phi_w, double phi_s, double meanPackingRatio,
                          double optPackingRatio, double cSAV, UnitsType units = ModelUnits);
double ResidenceTime(double cSAV, UnitsType units = ModelUnits);
double HeatPerUnitArea(double I_R, double t_r);
double ByramsFirelineIntensity(double H_A, double R);
double ByramsFlameLength(double I_B, UnitsType units = ModelUnits);

#endif
