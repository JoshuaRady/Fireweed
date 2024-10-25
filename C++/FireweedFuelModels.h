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
 * @brief This class represents the fuel properties and quantities for a fuel bed.
 * 
 * FuelModel objects are used to hold and manipulate fire behavior fuel models for use with the
 * Rothermel Albini fire spread model.
 * 
 * Similar fuel models are used for other purposes such as fire danger.  This class could be
 * expanded in the future for those applications.
 *
 * @par Number of Fuel Types:
 *   The Rothermel Albini fire spread model includes calculations for homogenous (single fuel) and
 * heterogenous (multiple fuel) cases.  The standard fuel models have five fuel types and this
 * convention is commonly used for other related US models.  The number of fuel types is not fixed
 * and the number can vary in practice.
 *
 * @par
 *   Additionally, the number of fuel types in a fuel model may not represent the actual number of
 * fuel types.  For example some of the standard fuel models, 1 and 3 from the original 13, are in
 * fact homogenous.  In published tables they are presented with empty columns.  Other fuel models
 * leave some fuel columns empty as well.
 *
 * @par
 * Further, application of curing for dynamic fuel models introduces a sixth dead herbaceous fuel
 * class (assuming the five standard staring classes).
 *
 * @par
 * Our fuel model implementation allows for an arbitrary number of fuel types.
 *
 * @par Homogeneous Fuels:
 * The homogeneous fuels case can be represented as a FuelModel object with a single fuel class.
 * Such a fuel model can be fed into the homogeneous versions of the spread model equations.  A
 * complicating factor is the fact that the homogenous fuels version of the spread rate calculation
 * use slightly different, unsubscripted notation than the heterogenous versions.
 *
 * @par
 * Attempts were made to include both notations in the FuelModel object with some problems.  For
 * convenience homogeneous notation aliases are provided for the moisture of extinction and fuel
 * particle properties.  Aliases were used to prevent homogenous and heterogenous members from
 * getting out of sync.  These aliases turned out to be a problem because they can't be updated when
 * the heterogeneous members are resized.  This is a ubiquitous problem because currently the
 * members start empty and so have to be resized during initialization.
 *
 * @par
 * The homogeneous members remain pending revision but should not be used.  If we can't find a way
 * to make them work they could be replaced with getter functions.
 *
 * @par Fuel Moisture:
 * Since fuel moisture is highly dynamic it is not included in fuel model tables and traditionally
 * may not be treated as part of a fuel model itself.  We have made the choice to optionally allow
 * fuel moisture to be recorded in the FuelModel object.  There are several reasons for this.  First
 * is is convient to keep all the fuel properties together.  More importantly the fuel moisture
 * values are linked with number and types of fuels.  Adding them to the fuel model allows us to
 * check them and make sure they are valid and consistent the fuels.  Finally, and perhaps most
 * importantly, this was necessary to represent dynamic fuel model behavior since calculating dynamic
 * fuel model moisture changes the structure of the fuel model itself.
 *
 * @par
 * The fuel moisture representation is currently focused on heterogeneous fuel moisture.  The
 * homogeneous case is less essential since M_f is a single value that is static.  It may be added
 * in the future.
 *
 * @par Fuel Model Properties:
 * For more information on fuel model data members see the documentation in FireweedRAFireSpread.cpp.
 * Add more here!!!!! ...
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

	//Fuel Type Index Functions:
	int NumDeadClasses() const;
	int NumLiveClasses() const;
	int LiveHerbaceousIndex() const;
	int LiveWoodyIndex() const;
	int DeadHerbaceousIndex() const;
	bool LiveHerbaceousPresent() const;
	bool LiveWoodyPresent() const;

	//Fuel Moisture functions:
	void SetFuelMoisture(std::vector <double> M_f_ij);
	std::vector <double> GetM_f_ij();
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
int FuelClassIndex(std::vector<int> liveDead, int liveDeadCat, int sizeIndex);

#endif
