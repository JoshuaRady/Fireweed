/***************************************************************************************************
FireweedFuelModels.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 19

Description:--------------------------------------------------------------------------------------
	This file is part of the Fireweed fire code library.  It contains a C++ implementation of a
fire behavior fuel object for use with the Rothermel Albini (Rothermel 1972, Albini 1976) fire
spread model and related models.

	This file is documented using Doxygen formatted comments.
***************************************************************************************************/

/** Fuel model table file format
 *
 * The fuel model input fiel contains data for the standard fuel models in a CSV format.
 * [The latest version is D3.]
 * The file starts with three lines of human readable header information followed by a machine
 * readable header line containing condensed column names.  The column names must match those
 * expected but the code below has been written to not assume a specific column order.
 *
 * Fields:
 * ...
 * All field values represent decimal (double) data other than the Number (int), Code (string
 * without whitespace), Name (string with whitespace), and Type (enum) fields.
 * 
 * Some models are do not have one or both live fuel classes.  This is indicated in the file as SAV
 * values of NA.  These values are converted to 0s in the FuelModel representation.  A SAV of 0 is
 * meaningless and was chosen because it can, unlike NA, be used in spread rate calculations without
 * causing problems.
 * The weights for any missing SAV class should always be zero [and we should enforce that!!!!!].
 *
 * [Find other notes on missing value notation.]
 *
 * ToDo:
 * - Consider adding version information or a format specifier to the file format.
 */

#include <algorithm>//For find(), count().
#include <fstream>
#include <iostream>
#include <math.h>//For pow()
#include <sstream>
#include <vector>

#include "FireweedFuelModels.h"
#include "FireweedMessaging.h"
#include "FireweedStringUtils.h"
#include "FireweedUtils.h"

/** Default (empty) constructor
 *
 * By default create a generic empty / no fuel loading model.
 */
FuelModel::FuelModel()
{
	Initialize();
}

/** File constructor: Initialize the fuel model specified by number from the specified file.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelFilePath The path to the CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 *                         in the Rothermel & Albini spread model.
 */
FuelModel::FuelModel(const std::string& fuelModelFilePath, int modelNumber,
                     bool spreadModelUnits)
{
	Initialize();
	LoadFromDelimited(fuelModelFilePath, modelNumber, "", spreadModelUnits);
}

/** File constructor: Initialize the fuel model specified by code from the specified file.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelFilePath The path to the CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 *                         in the Rothermel & Albini spread model.
 */
FuelModel::FuelModel(const std::string& fuelModelFilePath, std::string modelCode,
                     bool spreadModelUnits)
{
	Initialize();
	LoadFromDelimited(fuelModelFilePath, -1, modelCode, spreadModelUnits);
}

/** Convert the fuel models units.
 *
 * @param newUnits The units to convert to.  If the same as the current units do nothing.
 */
void FuelModel::ConvertUnits(UnitsType newUnits)
{
	if (units == newUnits)
	{
		Warning("Fuel model is already in requested units.");
	}
	else
	{
		if (newUnits == US)
		{
			//Leave Number, Code, Name, and Type unchanged.
			//Cured and NumClasses don't change.
	
			delta = delta / mPerFt;//m -> ft
	
			//Leave liveDead unchanged.
			//M_x and M_x_1 are either unitless fractions or percentages and can be left unchanged.
	
			//h = h / kJPerBtu * kgPerLb;//kJ/kg -> Btu/lb
	
			//S_T, S_T_ij, S_e, and S_e_ij are unitless fractions and can be left unchanged.
	
			//rho_p = rho_p / lbPerFtCuToKgPerMCu;
	
			for (int i = 0; i < numClasses; i++)
			{
				SAV_ij[i] = SAV_ij[i] * cmPerFt;//1/cm -> 1/ft | cm^2/cm^3 -> 3ft^2/ft^3
				w_o_ij[i] = w_o_ij[i] / kgPerLb * pow(mPerFt, 2);//kg/m^2 -> lb/ft^2
				h_ij[i] = h_ij[i] / kJPerBtu * kgPerLb;//kJ/kg -> Btu/lb
				rho_p_ij[i] = rho_p_ij[i] / lbPerFtCuToKgPerMCu;
			}
	
			cSAV = cSAV * cmPerFt;//cm^2/cm^3 -> ft^2/ft^3
	
			bulkDensity = bulkDensity / lbPerFtCuToKgPerMCu;//kg/m^3 -> lb/ft^3
	
			//RelativePackingRatio is a dimensionless ratio.
			//M_f_ij is a fraction.
			//Curing is a percent.
	
			//The metric units are alway kg/m^2. We don't include an option to convert to ton/Ac:
			w_o_Units = lbPer_ft2;
			units = US;
		}
		else//if (newUnits == Metric)
		{
			//Leave Number, Code, Name, and Type unchanged.
			//Cured and NumClasses don't change.
	
			delta = delta * mPerFt;//ft -> m
	
			//Leave liveDead unchanged.
			//M_x and M_x_1 are either unitless fractions or percentages and can be left unchanged.
	
			//h = h * kJPerBtu / kgPerLb;//Btu/lb -> kJ/kg
	
			//S_T, S_T_ij, S_e, and S_e_ij are unitless fractions and can be left unchanged.
	
			//rho_p = rho_p * lbPerFtCuToKgPerMCu;
	
			for (int i = 0; i < numClasses; i++)
			{
				SAV_ij[i] = SAV_ij[i] / cmPerFt;//1/ft -> 1/cm | ft^2/ft^3 -> cm^2/cm^3
				
				if (w_o_Units == tonPerAc)//Move out of loop?????
				{
					w_o_ij[i] = w_o_ij[i] * tonsPerAcToLbPerSqFt;//Convert loadings to lb/ft^2.
				}
				//We could check for invalid w_o_Units here.
				w_o_ij[i] = w_o_ij[i] * kgPerLb / pow(mPerFt, 2);//lb/ft^2 -> kg/m^2
	
				h_ij[i] = h_ij[i] * kJPerBtu / kgPerLb;//Btu/lb -> kJ/kg
				rho_p_ij[i] = rho_p_ij[i] * lbPerFtCuToKgPerMCu;
			}
	
			cSAV = cSAV / cmPerFt;//ft^2/ft^3 -> cm^2/cm^3
	
			bulkDensity = bulkDensity * lbPerFtCuToKgPerMCu;//lb/ft^3 -> kg/m^3
	
			//RelativePackingRatio is a dimensionless ratio.
			//M_f_ij is a fraction.
			//Curing is a percent.
	
			w_o_Units = kgPer_m2;
			units = Metric;
		}
	}
}

/** Print the fuel model data to an output stream.
 *
 * @param output The output stream to print to.
 *
 * Note: The homogenous aliases are pointing to bad addresses or something because the heterogeneous
 * members start as empty.  This needs to be fixed.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& FuelModel::Print(std::ostream& output) const
{
	output << "Fuel model:" << std::endl;
	output << "Number: " << number << std::endl;
	output << "Code: " << code << std::endl;
	output << "Name: " << name << std::endl;
	
	if (type == Static)
	{
		output << "Type: Static" << std::endl;
	}
	else
	{
		output << "Type Dynamic" << std::endl;
	}
	
	if (units == US)
	{
		output << "Units: US" << std::endl;
	}
	else
	{
		output << "Units: Metric" << std::endl;
	}

	//Should add member to indicate units of fuel loading.
	if (cured)
	{
		output << "Cured: true" << std::endl;
		output << "Percent curing: " << curing << std::endl;
	}
	else
	{
		output << "Cured: false" << std::endl;
	}

	output << "numClasses: " << numClasses << std::endl;

	output << "SAV_ij: ";
	PrintVector(output, SAV_ij);

	//This isn't in member order but it make sense to put it with the variable itself:
	output << "w_o_Units: ";
	if (w_o_Units == lbPer_ft2)
	{
		output << "lb/ft^2" << std::endl;
	}
	else if (w_o_Units == tonPerAc)
	{
		output << "ton/Ac" << std::endl;
	}
	else if (w_o_Units == kgPer_m2)
	{
		output << "kg/m^2" << std::endl;
	}
	//else error !!!!!!

	output << "w_o_ij: ";
	PrintVector(output, w_o_ij);
	
	output << "Delta: " << delta << std::endl;

	output << "LiveDead: ";
	bool firstClass = true;
	for (int ld : liveDead)
	{
		//Don't print a hanging comma:
		if (!firstClass)
		{
			output << ", ";
		}
		else
		{
			firstClass = false;
		}

		if (ld == Dead)
		{
			output << "Dead";
		}
		else//(ld == Live)
		{
			output << "Live";
		}
	}
	output << std::endl;

	output << "M_x_Units: ";
	if (M_x_Units == Percent)
	{
		output << "Percent" << std::endl;
	}
	else if (M_x_Units == Fraction)
	{
		output << "Fraction" << std::endl;
	}
	//else error!!!!!

	output << "M_x / M_x_1: " << M_x_1 << std::endl;

	//Note: The homogeneous aliases are omitted for now. !!!!!
	//output << "h: " << h << std::endl;//Not really necessary.
	output << "h_ij: ";
	PrintVector(output, h_ij);

	//output << "S_T: " << S_T << std::endl;//Not really necessary.
	output << "S_T_ij: ";
	PrintVector(output, S_T_ij);
	
	//output << "S_e: " << S_e << std::endl;//Not really necessary.
	output << "S_e_ij: ";
	PrintVector(output, S_e_ij);

	//output << "rho_p: " << rho_p << std::endl;//Not really necessary.
	output << "rho_p_ij: ";
	PrintVector(output, rho_p_ij);

	//These are precalcualted values populated when we load from a file.  These are not recalculated
	//when the fuel model is updated so could be out of date.  They could be removed in the future.
	//For now we label them to be clear.
	output << "cSAV: " << cSAV << " (precalculated)" << std::endl;
	output << "bulkDensity: " << bulkDensity << " (precalculated)" << std::endl;
	output << "relativePackingRatio: " << relativePackingRatio << " (precalculated)" << std::endl;
	
	if (M_f_ij.empty())
	{
		output << "M_f_ij not set." << std::endl;
	}
	else
	{
		output << "M_f_ij: ";
		PrintVector(output, M_f_ij);
	}
	
	return output;
}

//--------------------------------------------------------------------------------------------------
/** @par Fuel Type Index Functions:
 * For standard fuel models the number and order of fuel types is fixed but there are complicating
 * factors in practice.  First, many of the fuel models only represent either one of the herbaceous
 * or woody fuel types.  Both index positions will be present but the underrepresented class should
 * not be used without making changes to the model.  Adding fuel loadings for these 'missing' size
 * classes should not prevent the Rothermel Albini spread model from running but the fuel will
 * effectively be ignored and no have impact on fire behavior output, though this is not well
 * tested.  Second, when dynamic fuel models are used the number of fuel types and indexes may
 * change.
 *
 * The public liveDead member can be used to infer the fuel categories and the FuelClassIndex()
 * function can be used to get indexes by category more conveniently.  Still additional assumptions
 * must be used to know exactly what each index represents.  Additionally we could change the way we
 * represent fuel types in the future.  For example we may add a live moss category in the future.
 *
 * To futureproof the code and to reduce the need for calling code to understand the way the fuel
 * types are organized the following accessor and information functions are provided.
 *
 * @note These are currently C++ only.
 */

/** 
 * We have the public numClasses member to give the number of fuel types (numTypes might be a
 * better name).  It is provided in part for convenience since the number of fuels can be inferred
 * from the length of many of the data members.  I can also be used as the standard to check against
 * to ensure other members are correct.  However, I have worried about having it be public since it
 * should be read only.  It would be safer to make it private, if slightly less convenient.
 */
// int FuelModel::NumFuelTypes() const
// {
// 	return numClasses;
// }

/** How many size classes are there in the dead fuel category?
 *
 */
int FuelModel::NumDeadClasses() const
{
	return std::count(liveDead.begin(), liveDead.end(), Dead);
}

/** How many size classes are there in the live fuel category?
 *
 */
int FuelModel::NumLiveClasses() const
{
	return std::count(liveDead.begin(), liveDead.end(), Live);
}

/** Return the index (k) of the live herbaceous fuel type.
 *
 */
int FuelModel::LiveHerbaceousIndex() const
{
	return FuelClassIndex(liveDead, Live, 0);//In the future the 0 index may not be guaranteed.
}

/** Return the index (k) of the live woody fuel type.
 *
 */
int FuelModel::LiveWoodyIndex() const//Or just WoodyIndex()?
{
	return FuelClassIndex(liveDead, Live, 1);//In the future the 1 index may not be guaranteed.
}

/** Return the index (k) of the dead herbaceous fuel type, if present.  If not present returns -1.
 * Check cured first to see if it is present.
 */
int FuelModel::DeadHerbaceousIndex() const
{
	if (cured)
	{
		//We put the dead herbaceous fuel in the second dead position.  See DynamicFuelCuringCore().
		//This could change in the future.
		return FuelClassIndex(liveDead, Dead, 1);
	}
	else
	{
		return -1;
	}
}

/** Is the live herbaceous fuel type active in this fuel model?
 *
 */
bool FuelModel::LiveHerbaceousPresent() const
{
	//We indicate that a live fuel is not present with a SAV of 0, which is unique to our implementation. 
	return (SAV_ij[LiveHerbaceousIndex()] != 0);
}

/** Is the live woody fuel type active in this fuel model?
 *
 */
bool FuelModel::LiveWoodyPresent() const
{
	//We indicate that a live fuel is not present with a SAV of 0, which is unique to our implementation. 
	return (SAV_ij[LiveWoodyIndex()] != 0);
}

/** Is the dead herbaceous fuel type present in this fuel model?
 * We could be include this for completeness but checking cured gives the same result.
 */
// bool FuelModel::DeadHerbaceousPresent() const
// {
// 	return cured;
// }

//Fuel Moisture Functions:--------------------------------------------------------------------------

/** Record the fuel moisture content values.
 *
 * @param M_f_ij Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
 *
 * This function checks M_f_ij and stores it if valid, without apply curing.  If curing has been
 * previously applied the value can not be overwritten and an error is generated.  If the value has
 * been previously set but curing has not been applied we allow the value to be overwritten.  The
 * utility of overwriting is uncertain so we currently warn when this happens.
 *
 * @note This function could overloaded to take a sigle value (M_f) for the homogenous case.
 */
void FuelModel::SetFuelMoisture(std::vector <double> M_f_ij)
{
	//Check length is appropriate:
	if (M_f_ij.size() != numClasses)
	{
		//Stop("M_f_ij not of proper length.");//Add lengths of inputs and numClasses?????
		Stop("M_f_ij has length " + std::to_string(M_f_ij.size()) + " while numClasses = " + std::to_string(numClasses));
	}

	//Confirm values are valid:
	for (int i = 0; i < M_f_ij.size(); i++)
	{
		if (M_f_ij[i] < 0.0)//Clearly invalid:
		{
			Stop("M_f_ij contains negative value(s): " + VectorToStr(M_f_ij));
		}
		else
		{
			//The maximum fuel moisture varies by fuel type but also by moisture model:
			//These limits could be researched further but making this a warning lowers the stakes.
			const double maxDeadFM = 1.5;//100% may be fine but go higher.
			const double maxHerbFM = 2.5;//GSI maxes out at 250%.
			const double maxWoodyFM = 2.0;//GSI maxes out at 200%.

			if (liveDead[i] == Dead && M_f_ij[i] > maxDeadFM)
			{
				Warning("Dead fuel moisture value seems high: " + VectorToStr(M_f_ij));
			}
			else if (i == LiveHerbaceousIndex() && M_f_ij[i] > maxHerbFM)
			{
				Warning("Herbaceous fuel moisture value seems high: " + VectorToStr(M_f_ij));
			}
			else if (i == LiveWoodyIndex() && M_f_ij[i] > maxWoodyFM)
			{
				Warning("Woody fuel moisture value seems high: " + VectorToStr(M_f_ij));
			}
		}
	}

	//Checking if curing has already been calculated for this fuel model:
	if (cured == true)
	{
		Stop("Fuel model has already had curing applied.");
	}
	
	if (!this->M_f_ij.empty())
	{
		Warning("M_f_ij is being overwritten.");
	}

	this->M_f_ij = M_f_ij;
}

/** The (heterogenous) fuel moisture.
 *
 * The vector returned may be empty if the fuel moisture has not been assigned.  The calling code
 * should check for this.
 *
 * @note This function is specific about returning M_f_ij, which is alway a vector.  Homogenous
 * moisture M_f is a scalar which could be added in the future.
 */
std::vector <double> FuelModel::GetM_f_ij()//Or just M_f_ij()?????
{
	return M_f_ij;
}

/** Calculate and apply the curing of herbaceous fuels based on the herbaceous fuel moisture (per Scott & Burgan 2005).
 *
 * @param M_f_ij Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
 *
 * For dynamic fuels curing moves some live herbaceous fuel to a new dead herbaceous fuel class.
 * As a result the number of fuel classes may increase with this call.
 *
 * Trying to apply curing to a static fuel model has no effect except to store the moisture values.
 * By default the code posts a warning when this occurs.
 */
void FuelModel::CalculateDynamicFuelCuring(std::vector <double> M_f_ij, bool warn)
{
	SetFuelMoisture(M_f_ij);//Check validity of M_f_ij and save.

	if (type == Dynamic)
	{
		//The live herbaceous is the first dead fuel.  The index should be 3 (base 0) for standard fuel models:
		//int liveHerbIndex = std::find(liveDead.begin(), liveDead.end(), Live);
		std::vector <int>::iterator liveHerbIt = std::find(liveDead.begin(), liveDead.end(), Live);
		int liveHerbIndex = liveHerbIt - liveDead.begin();
		
		//Curing is a function of live herbaceous fuel moisture:
		double M_f_21 = M_f_ij[liveHerbIndex];

		//Calculate the transfer of herbaceous fuel loading from live to dead:
		//Curing is 0 at 120% fuel moisture and 1 at 30%:
		//Equation from Andrews 2018 pg. 36:
		//T = -1.11(M_f)_21 + 1.33, 0<=T<=1.0
		//This implements this exactly but is not numerically exact at the ends of the curing range:
		//cureFrac = -1.11 * M_f_21 + 1.33
		//Note: We can't use T, since T = TRUE in R.
		//This is exact at 30 and 120%:
		double cureFrac = 1.0 - ((M_f_21 - 0.3) / 0.9);//1.2 - 0.3 = 0.9
		
		if (cureFrac < 0.0)
		{
			cureFrac = 0.0;
			//Could break out here as no curing needs to be applied.
		}
		else if (cureFrac > 1.0)
		{
			cureFrac = 1.0;
		}

		DynamicFuelCuringCore(cureFrac);
	}
	else if (warn)
	{
		 Warning("Fuel model is static. No curing applied.");
	}
}

/** Calculate and apply the curing of herbaceous fuels based on percent curing.
 *
 * @param curing The percent herbaceous fuel curing to apply.
 *
 * This function is provided as an alternative to specifying moisture content.  Curing is normally
 * calculated based on the live herbaceous fuel moisture but it can be useful to specify the curing
 * directly, especially for testing.
 */
void FuelModel::CalculateDynamicFuelCuring(double curing, bool warn)
{
	if (type == Dynamic)
	{
		//Curing is generally presented as a percentage and that is what we expect:
		if (curing < 0 || curing > 100)
		{
			Stop("We expect fuel curing as a percentage.");
		}
		
		DynamicFuelCuringCore(curing / 100);//Convert to a fraction.
	}
	else if (warn)
	{
		 Warning("Fuel model is static. No curing applied.");
	}
}

//Private Functions:--------------------------------------------------------------------------------

/** Initialize members to a consistant 'no fuel' initial state.
 *
 * This function sets the fuel model values to reasonable default values.  It assumes the model is
 * structured like a standard fuel model with five fuel types.  It sets fuel loading to zero and
 * parameters to values that produce a fire spread rate of zero.
 *
 * This initial fuel model state should be further modified.  For a fuel model to be useful it needs
 * to have its parameters set to meaningful values, either from a file or manually.
 */
void FuelModel::Initialize()
{
	//Alternatively we could default the IDs to one of the nonburnable fuel models:
	number = -1;//Set to impossible value.
	code = "NA";//Alphanumeric code identifying the model.
	name = "Empty fuel model";//Descriptive name.
	type = Static;
	units = US;//Default to the units of the original papers.
	cured = false;
	numClasses = 5;//Standard fuel model.

	w_o_Units = lbPer_ft2;//Model units.
	M_x_Units = Fraction;//Model units.

	SAV_ij.resize(5, 0);//All fuel classes are absent.
	w_o_ij.resize(5, 0);//Important.

	delta = 0;//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
	liveDead.assign({Dead, Dead, Dead, Live, Live});//Standard fuel model.
	M_x_1 = 50;//Tough to choose a default here.  50% for now.

	h_ij.resize(5, 8000);//Standard for all but one standard fuel model.
	S_T_ij.resize(5, 0.0555);//Standard for all standard fuel models.
	S_e_ij.resize(5, 0.010);//Standard for all standard fuel models.
	rho_p_ij.resize(5, 32);//Standard for all standard fuel models.

	cSAV = 0;
	bulkDensity = 0;
	relativePackingRatio = 0;
	
	//M_f_ij is left as empty.
	
	curing = 0;//Percent curing = none.  When cured = false there should be no reason to look at this but we set a valid value.
}

/** Load a fuel model from the specified file.
 *
 * @param fuelModelFilePath The delimited text file containing the table of fuel models.  (See file
                            format specification.)
 * @param modelNumber The standard fuel model number of the fuel model requested.  -1 if not used.
 * @param modelCode The unique alphanumeric code of the fuel model requested.  Blank if not used.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 *                         in the Rothermel & Albini spread model.
 * @param delimiter The delimiter character.  Defaults to ',' for CSV.
 * 
 * Only the modelNumber or the modelCode should be passed in.  This is enforced through the calling
 * code. 
 *
 * It is assumed that the the FuelModel has the default five fuel classes.  All parameters will be
 * overwritten if the fuel model is found.  If not found the FuelModel will be returned unchanged.
 * We may change this to have a failed search throw an error.  Initializing the object to the
 * default state at the start of this routine also an option but this implementation may have some
 * advantages.
 *
 * Note: This code currently assumes the units of the file are in United States customary units with
 * loadings in ton/acre and moisture of extinction in percent.
 */
void FuelModel::LoadFromDelimited(const std::string& fuelModelFilePath, int modelNumber,
                                  std::string modelCode, bool spreadModelUnits, char delimiter)
{
	std::string line;//To hold lines of the input file.
	bool found = false;	
	int theModelNumber;
	std::string theModelCode;

	//Open the file:
	std::ifstream fmCSV(fuelModelFilePath);
	if (!fmCSV)
	{
		Stop("Couldn't open fuel model file at " + fuelModelFilePath);
	}

	//Skip the first 3 lines, which are human readable column headers.
	for (int i = 3; i > 0; i--)
	{
		std::getline(fmCSV, line);
	}
	
	//Get the parsable header line:
	std::getline(fmCSV, line);
	
	//Extract the column names from the header:
	std::vector<std::string> colNames = SplitDelim(line, delimiter);
	
	//Search rows until a match is found:
	while(std::getline(fmCSV, line))
	{
		std::stringstream lineStr(line);//This line as a stream.
		std::string field;//A field or token extracted from the line.

		//Given the code devoted to making the column order agnostic below we ought to search for these:

		//The first field is the fuel model number:
		std::getline(lineStr, field, delimiter);
		theModelNumber = stoi(field);

		//The second field is the fuel model code:
		std::getline(lineStr, field, delimiter);
		theModelCode = field;

		//Choose the selection method:
		if (modelNumber != -1)//Search with the model number:
		{
			//We could pre-check that the model number is in the range of valid values.
			if (theModelNumber == modelNumber)
			{
				found = true;
				break;
			}
		}
		else//If an invalid model number was passed search using the model code:
		{
			if (theModelCode == modelCode)
			{
				found = true;
				break;
			}
		}
	}

	//Extract fields from the matching row:
	if (found == true)
	{
		//std::cout << "Match found!" << std::endl;//Temp debugging.
		//std::cout << "Match found with model number: " << theModelNumber << std::endl;//Temp debugging!!!!!

		std::vector<std::string> fields = SplitDelim(line, delimiter, true);

		number = theModelNumber;
		code = theModelCode;
		//The remaining fields, other than the name and type, are numeric so could be converted here?

		units = US;//Assumed to always be the case for now.  Should be determined or set.
		cured = false;

		//Load the field values into the appropriate data members:
		//This is a bit of extra processing that allows us to not worry about the field order.
		for (int j = 0; j < fields.size(); j++)
		{
			//std::cout << j << ": " << fields[j] << std::endl;//Temp debugging!!!!!

			if (colNames[j] == "Name")
			{
				name = fields[j];
			}
			else if (colNames[j] == "Type")
			{
				if (fields[j] == "Static")
				{
					type = Static;
				}
				else if (fields[j] == "Dynamic")
				{
					type = Dynamic;
				}
				else//Error handling:
				{
					Stop("Invalid value for fuel model type.");
				}
			}
			else if (colNames[j] == "SAV_11")
			{
				SAV_ij[0] = stof(fields[j]);
			}
			else if (colNames[j] == "SAV_12")
			{
				SAV_ij[1] = stof(fields[j]);
			}
			else if (colNames[j] == "SAV_13")
			{
				SAV_ij[2] = stof(fields[j]);
			}
			else if (colNames[j] == "SAV_21")
			{
				if (fields[j] == "NA")
				{
					SAV_ij[3] = 0.0;
				}
				else
				{
					SAV_ij[3] = stof(fields[j]);
				}
			}
			else if (colNames[j] == "SAV_22")
			{
				if (fields[j] == "NA")
				{
					SAV_ij[4] = 0.0;
				}
				else
				{
					SAV_ij[4] = stof(fields[j]);
				}
			}
			else if (colNames[j] == "w_o_11")
			{
				w_o_ij[0] = stof(fields[j]);
			}
			else if (colNames[j] == "w_o_12")
			{
				w_o_ij[1] = stof(fields[j]);
			}
			else if (colNames[j] == "w_o_13")
			{
				w_o_ij[2] = stof(fields[j]);
			}
			else if (colNames[j] == "w_o_21")
			{
				w_o_ij[3] = stof(fields[j]);
			}
			else if (colNames[j] == "w_o_22")
			{
				w_o_ij[4] = stof(fields[j]);
			}
			else if (colNames[j] == "delta")
			{
				delta = stof(fields[j]);
			}
			else if (colNames[j] == "M_x")
			{
				M_x_1 = stof(fields[j]);
			}
			else if (colNames[j] == "h")
			{
				std::fill(h_ij.begin(), h_ij.end(), stof(fields[j]));
			}
			else if (colNames[j] == "S_T")
			{
				std::fill(S_T_ij.begin(), S_T_ij.end(), stof(fields[j]));
			}
			else if (colNames[j] == "S_e")
			{
				std::fill(S_e_ij.begin(), S_e_ij.end(), stof(fields[j]));
			}
			else if (colNames[j] == "rho_p")
			{
				std::fill(rho_p_ij.begin(), rho_p_ij.end(), stof(fields[j]));
			}
			else if (colNames[j] == "CharacteristicSAV")
			{
				cSAV = stof(fields[j]);
			}
			else if (colNames[j] == "BulkDensity")
			{
				bulkDensity = stof(fields[j]);
			}
			else if (colNames[j] == "RelativePackingRatio")
			{
				relativePackingRatio = stof(fields[j]);
			}
			else//Unrecognized field.  Report it?
			{
				//!!!!!
			}
		}

		//Typically the units of some parameters are different in the published tables than in the
		//model equations:
		//It would be an improvement to detect the units used the file.  That information is currently
		//in the header information but not in a form that would be ideal to parse.
		if (spreadModelUnits)
		{
			for (int i = 0; i < numClasses; i++)
			{
				w_o_ij[i] = w_o_ij[i] * lbsPerTon / ft2PerAcre;//ton/acre to lb/ft^2
			}
			M_x_1 = M_x_1 / 100;//% to fraction

			//Record the units used:
			w_o_Units = lbPer_ft2;
			M_x_Units = Fraction;
		}
		else
		{
			w_o_Units = tonPerAc;
			M_x_Units = Percent;
		}
	}
	else
	{
		//Not finding the fuel model is an error.
		//There is currently no way for the calling code to know that a failure occured and there
		//is little the code could do with an empty fuel model so it is better to throw an error.
		if (modelNumber != -1)//The search was by the model number:
		{
			Stop("No matching fuel model was found for number " + std::to_string(modelNumber));
		}
		else//The search used the model code:
		{
			Stop("No matching fuel model was found for number " + modelCode);
		}
	}

	fmCSV.close();
}

/** Perform the core curing calculation and update the fuel model.
 *
 * @param cureFrac The fractional curing of herbaceous fuel to apply.
 */
void FuelModel::DynamicFuelCuringCore(double cureFrac)
{
	//Checking if this function has already been run for this fuel model:
	//This check could be moved to start of the calling functions.
	if (cured == true)
	{
		Stop("Fuel model has already had curing applied.");
	}

	//The live herbaceous is the first live fuel.  The index should be 3 (base 0) for standard fuel models:
	//int liveHerbIndex = std::find(liveDead.begin(), liveDead.end(), Live)
	std::vector <int>::iterator liveHerbIt = std::find(liveDead.begin(), liveDead.end(), Live);
	int liveHerbIndex = liveHerbIt - liveDead.begin();

	//Expand the number of fuel classes, inserting the cured herbaceous at second dead position,
	//(Dead, 1) in 0 based space:
	w_o_ij.insert(w_o_ij.begin() + 1, 0);//Initial loading = 0.
	//Curing doesn't change SAV.  Inherit from the live herbaceous:
	SAV_ij.insert(SAV_ij.begin() + 1, SAV_ij[liveHerbIndex]);

	//For the standard fuel models all values for these parameters should be the same but don't
	//assume that for robustness.  Inherit from the live class:
	h_ij.insert(h_ij.begin() + 1, h_ij[liveHerbIndex]);
	S_T_ij.insert(S_T_ij.begin() + 1, S_T_ij[liveHerbIndex]);
	S_e_ij.insert(S_e_ij.begin() + 1, S_e_ij[liveHerbIndex]);
	rho_p_ij.insert(rho_p_ij.begin() + 1, rho_p_ij[liveHerbIndex]);

	//Expand liveDead:
	liveDead.insert(liveDead.begin() + 1, Dead);

	//Update the moisture content vector if present:
	if (!M_f_ij.empty())
	{
		M_f_ij.insert(M_f_ij.begin() + 1, M_f_ij[0]);//Inherit from 1-hr dead moisture.
	}

	numClasses = numClasses + 1;
	liveHerbIndex = liveHerbIndex + 1;//Update after all data members are restructured.

	//Transfer the loading from live to dead:
	w_o_ij[1] = cureFrac * w_o_ij[liveHerbIndex];
	w_o_ij[liveHerbIndex] = w_o_ij[liveHerbIndex] - w_o_ij[1];

	cured = true;//Record that curing has been applied.
	curing = cureFrac * 100;//Record the curing percentage.
}

//External functions:-------------------------------------------------------------------------------

/** Find a fuel model by number in the specified file and return it as a FuelModel object.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelFilePath The path to the CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 *                         in the Rothermel & Albini spread model.
 */
FuelModel GetFuelModelFromCSV(const std::string fuelModelFilePath, int modelNumber,
                              bool spreadModelUnits)
{
	FuelModel fm(fuelModelFilePath, modelNumber, spreadModelUnits);

	return fm;
}

/** Find a fuel model by alphanumeric code in the specified file and return it as a FuelModel object.
 *
 * @param modelCode The unique alphanumeric code of the fuel model requested.
 * @param fuelModelFilePath The path to the CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 *                         in the Rothermel & Albini spread model.
 */
FuelModel GetFuelModelFromCSV(const std::string fuelModelFilePath, std::string modelCode, 
                              bool spreadModelUnits)
{
	FuelModel fm(fuelModelFilePath, modelCode, spreadModelUnits);

	return fm;
}

/* Overloaded stream print operator for FuelModel.
 *
 */
std::ostream& operator<<(std::ostream& output, const FuelModel& fm)
{
	fm.Print(output);
	return output;
}

/** Return the index (k) in a variable array of the form X_ij given the (live/dead, size) index pair.
 *
 * Representing fuel model X_ij varaibles as vectors has the disadvantage of making mapping to
 * individual classes awkward.  This function makes this simple but remains somewhat inelegant.
 *
 * @param liveDead An array indicating if each index in the input variables represents a dead or
 *                 live fuel category.
 * @param liveDeadCat The Live / Dead catagory value to get.
 * @param sizeIndex The index (0 based) of the size class to get.
 *
 * @note This function is currently C++ only.
 */
int FuelClassIndex(std::vector<int> liveDead, int liveDeadCat, int sizeIndex)//Or FMClassIndex, ClassIndex
{
	int numDead = std::count(liveDead.begin(), liveDead.end(), Dead);
	int numLive = liveDead.size() - numDead;//or std::count(liveDead.begin(), liveDead.end(), Live);

	if (liveDeadCat == Dead)
	{
		if (sizeIndex < 0 || sizeIndex >= numDead)
		{
			Stop("Invalid dead size index.");
		}
		else
		{
			return sizeIndex;
		}
	}
	else if (liveDeadCat == Live)
	{
		if (sizeIndex < 0 || sizeIndex >= numLive)
		{
			Stop("Invalid live size index.");
		}
		else
		{
			return numDead + sizeIndex;
		}
	}
	else
	{
		Stop("Invalid live / dead category.");
	}

	return -1;//We could make the size index errors above warnings and return.
}
