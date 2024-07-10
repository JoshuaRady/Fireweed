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

***************************************************************************************************/

#include "FireweedFuelModels.h"

/** Default (empty) constructor
 *
 * A fuel model needs to initialized with data.  A new fuel model will be empty of data.  The
 * default constructor labels it as such and puts it in a no fuel state that will produce no fire
 * behavior.
 */
void FuelModel()
{
	number = -1;//Set to impossible value.
	code = "NA";//Alphanumeric code identifying the model.
	name "Empty fuel model";//Descriptive name.
	//type = ;
	
	numClasses = 0;
	units = "US";//Default to the units of the original papers.
	cured = FALSE;
	
	//Setting everything to 0 should allow the spread rate calculations to complete with a spread rate of zero.
	
	cSAV = 0;
	bulkDensity = 0;
	relativePackingRatio = 0;
}

//External functions:
//I think it makes the most sense to keep these functions outside the class.

/** Find a fuel model by number in the specified file and return it as a FuelModel object.
 *
 * @param modelNumber The unique number of the fuel model.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param originalUnits If true then the fuel model table file is in the original United States customary units.
 * @param expand If true expand properties provided as single values to the length of fuel classes. (Should always be true?)
 * 
 * Incomplete!!!!!
 */
FuelModel GetFuelModelFromCSV(int modelNumber, std::string fuelModelTableFile,//fuelModelPath = 
                              bool originalUnits, bool expand)
{
	FuelModel fm;

	//...

	return fm;
}

/** Find a fuel model by alphanumeric code in the specified file and return it as a FuelModel object.
 *
 * @param modelNumber The unique number of the fuel model.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param originalUnits If true then the fuel model table file is in the original United States customary units.
 * @param expand If true expand properties provided as single values to the length of fuel classes. (Should always be true?)
 *
 * I think it makes the most sense to keep these functions outside the class.
 * 
 * Incomplete!!!!!
 */
FuelModel GetFuelModelFromCSV(std::string modelCode, std::string fuelModelTableFile,
                              bool originalUnits, bool expand)
{
	FuelModel fm;

	//...

	return fm;
}
