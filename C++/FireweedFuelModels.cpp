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

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "FireweedFuelModels.h"

/** Default (empty) constructor
 *
 * A fuel model needs to initialized with data.  A new fuel model will be empty of data.  The
 * default constructor labels it as such and puts it in a no fuel state that will produce no fire
 * behavior.
 */
void FuelModel::FuelModel()
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



//Private functions:--------------------------------------------------------------------------------

/** Load a fuel model from the specified file.
 *
 * @param fuelModelTableFile The CSV file containing the table of fuel models.			Path?????
 * @param modelNumber The standard fuel model number of the fuel model requested.  -1 if not used.
 * @param modelCode The unique alphanumeric code of the fuel model requested.  Blank if not used.
 * @param originalUnits If true then the fuel model table file is in the original United States
 * customary units with loading in ton/acre..
 * @param expand If true expand properties provided as single values to the length of fuel classes. (Should always be true?)
 * 
 * Incomplete!!!!!
 */
void FuelModel::LoadFromCSV(const std::string& fuelModelTableFile,//fuelModelPath = 
                            int modelNumber, std::string modelCode,
                            bool originalUnits)//, bool expand)//Neither yet used!
{
	char delimiter = ',';
	std::string str;
	//FuelModel fm;//Return value.


	//Initialize members to a consistant state?????
	this.Initialize();





	//Open the file:
	std::ifstream fmCSV(fuelModelTableFile);
	//Error handling?????

	//Skip the first 3 lines, which are human readable column headers.
	for (i = 3; i > 0; i--)
	{
		std::getline(fmCSV, str);
	}
	
	//Get the parsable header line:
	std::getline(fmCSV, str);
	
	//Extract the column names from the header:
	std::vector<std::string> colNames = SplitDelim(str, delimiter);
	
	//Search rows until a match is found:
	bool found = false;
	while(std::getline(fmCSV, str))
	{
		std::stringstream strStr(str);//This line as a stream.
		std::string field;//Or token

		//Given the code devoted to making the column order agnostic below we ought to search for these:

		//The first field is the fuel model number:
		std::getline(strStr, field, delimiter);
		int theModelNumber = stoi(field)

		//The second field is the fuel model code:
		std::getline(strStr, field, delimiter);

		//Add conditional for the selection method...

		if (theModelNumber == modelNumber)
		{
			found = true;
			break;
		}
	}

	//Extract fields from the matching row:
	if (found == true)
	{
		//This may ignore the fields pulled off above, which will cause the indexes to not match!!!!!
		std::vector<std::string> fields = SplitDelim(strStr, delimiter);
		
		this.number = theModelNumber;
		this.code = field;
		//All the remaining fields, other than the name, are numeric so they can be converted here?
		
		//The fuel categories have the same live dead status for all standard fuel models:
		//int default = {Dead, Dead, Dead, Live, Live}
		//this.liveDead.assign()
		std::vector<int> initial = {Dead, Dead, Dead, Live, Live};//Move type defs to header!!!!!
		this.liveDead = initial;
		
		//Resize: Make default?
		SAV_ij.resize(5);
		w_o_ij.resize(5);
		h_ij.resize(5);
		S_T_ij.resize(5);
		S_T = S_T_ij[0];
		S_e_ij.resize(5);
		rho_p_ij.resize(5);
		
		this.numClasses = 5;
		
		//Load the field values into the appropriate data members:
		//This is a bit of extra processing that allows us the not worry about the field order.
		for (j = 0; j < sizeof(feilds); j ++)
		{
			if (colNames[j].compare("Name"))
			{
				this.name = fields[j];
			}
			else if (colNames[j].compare("Type"))
			{
				if (fields[j].compare("Static"))
				{
					this.type = Static;
				}
				else if (fields[j].compare("Dynamic"))
				{
					this.type = Dynamic;
				}
				//else//Error handling:
				//{
				//	Error("Invalid value for fuel model type.")
				//}
			}
			else if (colNames[j].compare("SAV_11"))
			{
				this.SAV_ij[0] = stof(fields[j]);
			}
			else if (colNames[j].compare("SAV_12"))
			{
				this.SAV_ij[1] = stof(atof(fields[j]));
			}
			else if (colNames[j].compare("SAV_13"))
			{
				this.SAV_ij[2] = stof(fields[j]);
			}
			else if (colNames[j].compare("SAV_21"))
			{
				if fields[j].compare("NA")
				{
					this.SAV_ij[3] = 0.0;
				}
				else
				{
					this.SAV_ij[3] = stof(fields[j]);
				}
			}
			else if (colNames[j].compare("SAV_22"))
			{
				if fields[j].compare("NA")
				{
					this.SAV_ij[4] = 0.0;
				}
				else
				{
					this.SAV_ij[4] = stof(fields[j]);
				}
			}
			else if (colNames[j].compare("w_o_11"))
			{
				this.w_o_ij[0] = stof(fields[j]);
			}
			else if (colNames[j].compare("w_o_12"))
			{
				this.w_o_ij[1] = stof(fields[j]);
			}
			else if (colNames[j].compare("w_o_13"))
			{
				this.w_o_ij[2] = stof(fields[j]);
			}
			else if (colNames[j].compare("w_o_21"))
			{
				this.w_o_ij[3] = stof(fields[j]);
			}
			else if (colNames[j].compare("w_o_22"))
			{
				this.w_o_ij[4] = stof(fields[j]);
			}
			else if (colNames[j].compare("delta"))
			{
				this.delta = stof(fields[j]);
			}
			else if (colNames[j].compare("M_x"))
			{
				this.M_x_1 = stof(fields[j]);
			}
			else if (colNames[j].compare("h"))//This code assumes we expand!!!!!
			{
				std::fill(this.h_ij, h_ij.begin(), h_ij.end(), stof(fields[j]));
			}
			else if (colNames[j].compare("S_T"))
			{
				std::fill(this.S_T_ij, S_T_ij.begin(), S_T_ij.end(), stof(fields[j]));
			}
			else if (colNames[j].compare("S_e"))
			{
				std::fill(this.S_e_ij, S_e_ij.begin(), S_e_ij.end(), stof(fields[j]));
			}
			else if (colNames[j].compare("rho_p"))
			{
				std::fill(this.rho_p_ij, rho_p_ij.begin(), rho_p_ij.end(), stof(fields[j]));
			}
			else if (colNames[j].compare("CharacteristicSAV"))
			{
				this.cSAV = stof(fields[j]);
			}
			else if (colNames[j].compare("BulkDensity"))
			{
				this.bulkDensity = stof(fields[j]);
			}
			else if (colNames[j].compare("RelativePackingRatio"))
			{
				this.relativePackingRatio = stof(fields[j]);
			}
			else//Unrecognized field.  Report it?
			{
				//
			}
		}
	}
	else
	{
		//Not finding the fuel model is probably an error.  At the least we should warn that no
		//match was found.
	}

	//Further processing... ?
	
	//originalUnits
	
	//Already set in defaults:
	//this.cured = false;
	//units = ...

	fmCSV.close();
	//return fm;
}

//External functions:
//I think it makes the most sense to keep these functions outside the class.

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

/** Find a fuel model by number in the specified file and return it as a FuelModel object.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param originalUnits If true then the fuel model table file is in the original United States
 * customary units with loading in ton/acre..
 * @param expand If true expand properties provided as single values to the length of fuel classes. (Should always be true?)
 * 
 * Incomplete!!!!!
 */
FuelModel GetFuelModelFromCSV(int modelNumber, std::string fuelModelTableFile,//fuelModelPath = 
                              bool originalUnits, bool expand)//Neither yet used!
{
	FuelModel fm;//Return value.




	return fm;
}

/** Find a fuel model by alphanumeric code in the specified file and return it as a FuelModel object.
 *
 * @param modelCode The unique alphanumeric code of the fuel model requested.
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


/** Split a delimited string into a vector of substrings.
 *
 * This should probably be moved to a utility file somewhere.
 */
std::vector<std::string> SplitDelim(const std::string& str, char delimiter)
{
	std::vector<std::string> substrings;
	std::stringstream strStr(str);
	std::string substring;//Or token
	
	while(getline(strStr, substring, delimiter))
	{
		substrings.pushback(substring);
	}
	
	return substrings;
}
