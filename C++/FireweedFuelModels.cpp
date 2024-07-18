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

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "FireweedFuelModels.h"

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
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 * in the Rothermel & Albini spread model.
 */
FuelModel::FuelModel(const std::string& fuelModelTableFile, int modelNumber,
                     bool spreadModelUnits)
{
	Initialize();
	LoadFromCSV(fuelModelTableFile, modelNumber, "", spreadModelUnits);
}

/** File constructor: Initialize the fuel model specified by code from the specified file.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 * in the Rothermel & Albini spread model.
 */
FuelModel::FuelModel(const std::string& fuelModelTableFile, std::string modelCode,
                     bool spreadModelUnits)
{
	Initialize();
	LoadFromCSV(fuelModelTableFile, -1, modelCode, spreadModelUnits);
}

/** Print the fuel model data to an output stream.
 *
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
	//output << "Cured: " << cured << std::endl;
	if (cured)
	{
		output << "Cured: true" << std::endl;
	}
	else
	{
		output << "Cured: false" << std::endl;
	}

	output << "numClasses: " << numClasses << std::endl;

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
	//else error
	
	output << "M_x_Units: ";
	if (M_x_Units == Percent)
	{
		output << "Percent" << std::endl;
	}
	else if (M_x_Units == Fraction)
	{
		output << "Fraction" << std::endl;
	}
	//else error

	//I prefer lists separated by commas but that requires a more complex loop.
	output << "SAV_ij: ";
	//for (int i; i > numClasses; i++)
	for (double sav : SAV_ij)
	{
		output << sav << " ";
	}
	output << std::endl;

	output << "w_o_ij: ";
	PrintVector(output, w_o_ij);
	
	output << "Delta: " << delta << std::endl;

	output << "LiveDead: ";
	//PrintVector(output, liveDead);//Will cast to double?
	for (int ld : liveDead)
	{
		if (ld == Dead)
		{
			output << "Dead" << " ";
		}
		else//(ld == Live)
		{
			output << "Live" << " ";
		}
	}
	output << std::endl;

	output << "M_x / M_x_1: " << M_x_1 << std::endl;

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

	output << "cSAV: " << cSAV << std::endl;
	output << "bulkDensity: " << bulkDensity << std::endl;
	output << "relativePackingRatio: " << relativePackingRatio << std::endl;
	
	return output;
}

//Private functions:--------------------------------------------------------------------------------

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
}

/** Load a fuel model from the specified file.
 *
 * @param fuelModelTableFile The CSV file containing the table of fuel models.			Path?????
 * @param modelNumber The standard fuel model number of the fuel model requested.  -1 if not used.
 * @param modelCode The unique alphanumeric code of the fuel model requested.  Blank if not used.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 * in the Rothermel & Albini spread model.
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
void FuelModel::LoadFromCSV(const std::string& fuelModelTableFile,//fuelModelPath = 
                            int modelNumber, std::string modelCode,
                            bool spreadModelUnits)
{
	char delimiter = ',';
	std::string line;//To hold lines of the input file.
	bool found = false;	
	std::string theModelCode;
	int theModelNumber;

	//Open the file:
	std::ifstream fmCSV(fuelModelTableFile);
	//Error handling?????

	//Skip the first 3 lines, which are human readable column headers.
	for (int i = 3; i > 0; i--)
	{
		std::getline(fmCSV, line);
	}
	
	//Get the parsable header line:
	std::getline(fmCSV, line);
	//std::cout << "line: " << line << std::endl;//Temp debugging.
	
	//Extract the column names from the header:
	std::vector<std::string> colNames = SplitDelim(line, delimiter);
	
	//std::cout << "colNames" << colNames << std::endl;//Temp debugging.
	std::cout << "colNames: " << std::endl;//Temp debugging:
	for (std::string colName : colNames)
	{
		std::cout << colName << ", ";
	}
	std::cout << std::endl;//
	
	//Search rows until a match is found:
	while(std::getline(fmCSV, line))
	{
		std::stringstream lineStr(line);//This line as a stream.
		std::string field;//A field or token extracted from the line.

		//Given the code devoted to making the column order agnostic below we ought to search for these:

		//The first field is the fuel model number:
		std::getline(lineStr, field, delimiter);
		theModelNumber = stoi(field);
		std::cout << "field: " << field << std::endl;//Temp debugging.
		std::cout << "theModelNumber: " << theModelNumber << std::endl;//Temp debugging.

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
			//if (modelCode.compare(field))
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
		std::cout << "Match found!" << std::endl;//Temp debugging.
		
		//std::stringstream lineStr2(line);//The matching line as a stream.  We start over so all fields are included.
		//std::vector<std::string> fields = SplitDelim(lineStr2, delimiter);
		std::vector<std::string> fields = SplitDelim(line, delimiter);
		
		std::cout << "fields: " << std::endl;//Temp debugging:
		for (std::string theField : fields)
		{
			//std::cout << theField << ", ";
			std::cout << theField << ", " << std::endl;
		}
		//std::cout << std::endl;//

		this->number = theModelNumber;
		this->code = theModelCode;
		//The remaining fields, other than the name and type, are numeric so could be converted here?

		units = US;//Assumed to always be the case for now.  Should be determined or set.
		this->cured = false;

		//Load the field values into the appropriate data members:
		//This is a bit of extra processing that allows us the not worry about the field order.
		//for (int j = 0; j < sizeof(fields); j ++)
		for (int j = 0; j < fields.size(); j++)
		{
			std::cout << j << std::endl;//Temp debugging.
			
			if (!colNames[j].compare("Name"))
			{
				this->name = fields[j];
			}
			else if (!colNames[j].compare("Type"))
			{
				if (!fields[j].compare("Static"))
				{
					this->type = Static;
				}
				else if (!fields[j].compare("Dynamic"))
				{
					this->type = Dynamic;
				}
				//else//Error handling:
				//{
				//	Error("Invalid value for fuel model type.")
				//}
			}
			else if (!colNames[j].compare("SAV_11"))
			{
				this->SAV_ij[0] = stof(fields[j]);
			}
			else if (!colNames[j].compare("SAV_12"))
			{
				this->SAV_ij[1] = stof(fields[j]);
			}
			else if (!colNames[j].compare("SAV_13"))
			{
				this->SAV_ij[2] = stof(fields[j]);
			}
			else if (!colNames[j].compare("SAV_21"))
			{
				if (!fields[j].compare("NA"))
				{
					this->SAV_ij[3] = 0.0;
				}
				else
				{
					this->SAV_ij[3] = stof(fields[j]);
				}
			}
			else if (!colNames[j].compare("SAV_22"))
			{
				if (!fields[j].compare("NA"))
				{
					this->SAV_ij[4] = 0.0;
				}
				else
				{
					this->SAV_ij[4] = stof(fields[j]);
				}
			}
			else if (!colNames[j].compare("w_o_11"))
			{
				this->w_o_ij[0] = stof(fields[j]);
			}
			else if (!colNames[j].compare("w_o_12"))
			{
				this->w_o_ij[1] = stof(fields[j]);
			}
			else if (!colNames[j].compare("w_o_13"))
			{
				this->w_o_ij[2] = stof(fields[j]);
			}
			else if (!colNames[j].compare("w_o_21"))
			{
				this->w_o_ij[3] = stof(fields[j]);
			}
			else if (!colNames[j].compare("w_o_22"))
			{
				this->w_o_ij[4] = stof(fields[j]);
			}
			else if (!colNames[j].compare("delta"))
			{
				this->delta = stof(fields[j]);
			}
			else if (!colNames[j].compare("M_x"))
			{
				this->M_x_1 = stof(fields[j]);
			}
			else if (!colNames[j].compare("h"))
			{
				std::fill(h_ij.begin(), h_ij.end(), stof(fields[j]));
			}
			else if (!colNames[j].compare("S_T"))
			{
				std::fill(S_T_ij.begin(), S_T_ij.end(), stof(fields[j]));
			}
			else if (!colNames[j].compare("S_e"))
			{
				std::fill(S_e_ij.begin(), S_e_ij.end(), stof(fields[j]));
			}
			else if (!colNames[j].compare("rho_p"))
			{
				std::fill(rho_p_ij.begin(), rho_p_ij.end(), stof(fields[j]));
			}
			else if (!colNames[j].compare("CharacteristicSAV"))
			{
				this->cSAV = stof(fields[j]);
			}
			else if (!colNames[j].compare("BulkDensity"))
			{
				this->bulkDensity = stof(fields[j]);
			}
			else if (!colNames[j].compare("RelativePackingRatio"))
			{
				std::cout << "Match RelativePackingRatio" << std::endl;//Temp debugging.
				this->relativePackingRatio = stof(fields[j]);
				std::cout << "j" << j << std::endl;//Temp debugging.
				std::cout << "fields[j]" << fields[j] << std::endl;//Temp debugging.
				std::cout << "stof(fields[j])" << stof(fields[j]) << std::endl;//Temp debugging.
			}
			else//Unrecognized field.  Report it?
			{
				//
			}
		}

		//Typically the units of some parameters are different in the published tables than in the
		//model equations:
		//It would be an improvement to detect the units used the file.  That information is currently
		//in the header information but not in a form that would be ideal to parse.
		if (spreadModelUnits)
		{
			for (int i; i < numClasses; i++)
			{
				this->w_o_ij[i] = this->w_o_ij[i] * lbsPerTon / ft2PerAcre;//ton/acre to lb/ft^2
			}
			this->M_x_1 = this->M_x_1 / 100;//% to fraction

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
		//Not finding the fuel model is probably an error.  At the least we should warn that no
		//match was found.
	}

	fmCSV.close();
}

//External functions:-------------------------------------------------------------------------------

/** Find a fuel model by number in the specified file and return it as a FuelModel object.
 *
 * @param modelNumber The standard fuel model number of the fuel model requested.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 * in the Rothermel & Albini spread model.
 */
FuelModel GetFuelModelFromCSV(const std::string fuelModelTableFile, int modelNumber,
                              bool spreadModelUnits)
{
	FuelModel fm(fuelModelTableFile, modelNumber, spreadModelUnits);

	return fm;
}

/** Find a fuel model by alphanumeric code in the specified file and return it as a FuelModel object.
 *
 * @param modelCode The unique alphanumeric code of the fuel model requested.
 * @param fuelModelTableFile The CSV file containing the table of fuel models.
 * @param spreadModelUnits If true then convert units used in the file that differ from those used
 * in the Rothermel & Albini spread model.
 */
FuelModel GetFuelModelFromCSV(const std::string fuelModelTableFile, std::string modelCode, 
                              bool spreadModelUnits)
{
	FuelModel fm(fuelModelTableFile, modelCode, spreadModelUnits);

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
		substrings.push_back(substring);
	}
	
	return substrings;
}

/** Print a vector to an output stream on a single line with separators.
 *
 * Turn into a template?
 */
//std::ostream& PrintVector(std::ostream& output, std::vector <double> vec, std::string separator)
//std::ostream& PrintVector(std::ostream& output, const std::vector& <double> vec, const std::string separator)
//std::ostream& PrintVector(std::ostream& output, std::vector <double> vec, std::string separator)
//std::ostream& PrintVector(std::ostream& output, std::vector <double>& vec, std::string separator)
std::ostream& PrintVector(std::ostream& output, const std::vector <double>& vec, std::string separator)
{
	for (int i = 0; i < vec.size() - 1; i++)
	{
		output << vec[i] << separator;
	}
	output << vec[vec.size() - 1] << std::endl;

	return output;
}
