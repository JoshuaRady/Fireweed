/***************************************************************************************************
FireweedFuelModels.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  This header file declares a fire
behavior fuel object for use with the Rothermel Albini (Rothermel 1972, Albini 1976) fire spread
model and related models.

***************************************************************************************************/
#ifndef FIREWEEDFUELMODELS_H
#define FIREWEEDFUELMODELS_H

#include <algorithm>//For count().
#include <string>
#include <vector>
#include "FireweedUnits.h"

//Constants and Flags:------------------------------------------------------------------------------
enum FuelModelType {Static, Dynamic};

/* Fuel Category Constants
 *
 * The values of the live dead categories are forced to match the matching array indexes so they may
 * be used to access arrays of the form X_i (values are language specific).
 *
 * Note: These are temporarily duplicated from FireweedRAFireSpread.h.
*/
//enum FuelCategory {Dead = 0, Live = 1};//Having some issues with casting this.
const int Dead = 0;
const int Live = 1;

/* @class FuelModel
 * @brief The FuelModel class represents the fuel properties and quantities necessary for a location.
 * [More?????]
 *
 * For more information on fuel model data members see the documentation in FireweedRAFireSpread.cpp ...
 *
 * Number of fuel classes:
 *   The standard fuel models have five fuel classes but this is not strictly fixed.  Homogeneous
 * fuel models with one fuel type are described in the Rothermel 
 * Some of the standard fuel models, 1 and 3 from the original 13, are in fact homogenous.  In
 * published tables they are presented with empty columns.  Additionally application of curing to
 * dynamic fuel models introduces a sixth dead herbaceous fuel class.
 *   Out implementation allow for ana arbitrary number of fuel classes.  The liveDead flag and   is used
 * 
 *   The homogeneous fuels case can be represented as a FuelModel object with a single fuel class.
 * However, the homogenous fuels version of the spread rate calculation uses different,
 * unsubscripted, notation than the heterogenous version.  For convenience homogeneous notation
 * aliases are provided for the moisture of extinction and fuel particle properties.
 * Note: These aliases are turning out to be a problem because they can't be updated when the
 * heterogeneous members are resized.
 */
class FuelModel {
	public:

	//Fuel model identifiers:
	//Together the model number and code represent different model IDs.  For the 13 the code ~ number.
	int number;//The model's standard fuel model number.
	std::string code;//Alphanumeric code identifying the model.
	std::string name;//Descriptive name.
	
	//Model model properties:
	FuelModelType type;//Static vs. Dynamic
	UnitsType units;//The model units type.
	//Should add member to indicate units of fuel loading.
	bool cured;//For dynamic models, has curing been applied to the herbaceous fuels?
	int numClasses;//The number of fuel classes.  Can be inferred but...
	UnitCode w_o_Units;//Units used for w_o_ij.
	UnitCode M_x_Units;//Units used for M_x / M_x_1.

	//Fuel model parameters / data members:

	//Fuel array:
	std::vector <double> SAV_ij;//Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).

	std::vector <double> w_o_ij;//An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
	double delta;//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
	std::vector <int> liveDead;//The live / dead category of each fuel class.

	//Dead fuel moisture of extinction (fraction: water weight/dry fuel weight).
	double M_x_1;
	double& M_x = M_x_1;

	//The fuel moisture content (M_f, M_f_ij) is an environmental fuel property which is supplied
	//separately from the fuel model but could be added to the object.
	//It is also needed for curing/

	//Fuel particle properties:
	
	//Heat content of the fuel types (Btu/lb | kJ/kg):
	std::vector <double> h_ij;
	double& h = h_ij[0];

	//S_T_ij = Total mineral content for each fuel type (unitless fraction):
	std::vector <double> S_T_ij;
	double& S_T = S_T_ij[0];
	
	//Effective mineral content (unitless fraction (mineral mass – mass silica) / total dry mass):
	std::vector <double> S_e_ij;
	double& S_e = S_e_ij[0];

	//Fuel particle density (lb/ft^3 | kg/m^3):
	std::vector <double> rho_p_ij;
	double& rho_p = rho_p_ij[0];

	//Precalculated columns (could be removed):
	double cSAV;//Characteristic SAV of the fuel bed. [characteristicSAV in R!!!!!]
	double bulkDensity;
	double relativePackingRatio;

	//Constructors:
	FuelModel();
	FuelModel(const std::string& fuelModelFilePath, int modelNumber, bool spreadModelUnits = true);
	FuelModel(const std::string& fuelModelFilePath, std::string modelCode, bool spreadModelUnits = true);
	
	void ConvertUnits(UnitsType newUnits);

	std::ostream& Print(std::ostream& output) const;
	
	void SetFuelMoisture(std::vector <double> M_f_ij);
	void CalculateDynamicFuelCuring(std::vector <double> M_f_ij, bool warn = true);
	void CalculateDynamicFuelCuring(double curing, bool warn = true);

	private:
	//Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).  Initially empty:
	std::vector <double> M_f_ij;
	//The percent herbaceous fuel curing.  Will be 0 when cured = false:
	//Currently this only stored for reference and is not used for anything.
	double curing;
	
	void Initialize();//Init(), InitBlank()?
	void LoadFromDelimited(const std::string& fuelModelFilePath, int modelNumber,
	                       std::string modelCode, bool spreadModelUnits = true, char delimiter = ',');
	void DynamicFuelCuringCore(double cureFrac);
};

//External functions:

FuelModel GetFuelModelFromCSV(const std::string fuelModelFilePath, int modelNumber,
                              bool spreadModelUnits = true);
FuelModel GetFuelModelFromCSV(const std::string fuelModelFilePath, std::string modelCode,
                              bool spreadModelUnits = true);
std::ostream& operator<<(std::ostream& output, const FuelModel& fm);

std::vector<std::string> SplitDelim(const std::string& str, char delimiter);
std::vector<std::string> SplitDelim(const std::string& str, char delimiter, bool allowQuotes);
std::ostream& PrintVector(std::ostream& output, const std::vector <double>& vec, std::string separator = ", ");
std::string VectorToStr(const std::vector <double> vec, std::string separator = ", ");

#endif
