/***************************************************************************************************
FireweedFosbergFuelMoisture.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:----------------------------------------------------------------------------------------
  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
Fosberg dead fuel moisture model, specifically the Rothermel 1981 / NWCG variant.

References:-----------------------------------------------------------------------------------------

Fosberg, M.A. and Deeming, J.E.
Derivation of the 1-and 10-hour timelag fuel moisture calculations for fire-danger rating.
Research Note RM-RN-207. Rocky Mountain Forest and Range Experiment Station, Fort Collins, CO.
Forest Service, US Department of Agriculture, 1971.
  This is the first publication documenting the 'Fosberg' dead fuel moisture calculation.  It gives
calcuations for 1-hour and 10-hour fuels.  Subsequent publications present methods for 100-hour and
1000-hour fuels.  The Fosberg method was highly influential and is still used today, with some
modifications.

Rothermel, Richard C.
How to predict the spread and intensity of forest and range fires.
Gen. Tech. Report INT-143. Ogden, UT: U.S. Department of Agriculture, Forest Service, Intermountain
Forest and Range Experiment Station. 161 p., 1983.
  This presents a modified version of the tables presented in Fosberg and Deeming 1971 for use in
field estimation of 1-hour fuel moisture.  These tables remain in use today by the National
Wildfire Coordinating Group (NWCG) and are included in the NWCG Incident Response Pocket Guide
(IPTG, as of 2022) carried by wildland firefighters.

***************************************************************************************************/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include "FireweedFosbergFuelMoisture.h"
#include "FireweedMessaging.h"
#include "FireweedUtils.h"

//Code:---------------------------------------------------------------------------------------------

/** Convert a temperature in degrees Celsius to degrees Fahrenheit:
* 
* @param degreesC The temperature to convert (degrees Celsius).
* 
* @note This chould be moved either to a met utilities or to units.
 */
double CtoF(double degreesC)
{
  return degreesC * 9/5 + 32;
}

/** Return an estimate of the 1-hour dead fuel moisture based on the conditions passed.  The
 * modified Fosberg lookup table method presented in Rothermel 1981 and adopted by the NWCG is used.
 * 
 * @param tableA_Path - tableD_Path = Paths to tab delimited files holding the lookup table data.
 *   These tables are used in Rothermel 1981 and are also the NWCG website and the 2022 IPTG.
 *   These paths add to an already large number of parameters.  Their file location can not be safely
 *   assumed but there values could be consolidated in a list of placed in globals.
 * 
 * @param temp The air temperature at 4.5 feet / 1.37 meters above ground level (degrees F / C).
 * @param rh Relative humidity (percent).
 * @param monthOfYear The numeric month of year (1 - 12).
 * @param hourOfDay The numeric hour of day in 24 hour format (1 - 24).
 * @param slopePct Percent slope of the location (rise / run x 100%).
 * @param aspectCardinal The aspect of the slope as a cardinal direction (N, E, S, W).
 * @param shaded Does the location have 50% or greater canopy cover or is there full cloud cover.
 *               This could be converted to a percentage.
 * @param elevation A code indicating the slope position for the prediction relative to where the
 *                  weather conditions were taken.  Values are:
 *  B: The fire is 1000 - 2000 feet below the location where weather was recorded,
 *  L: The fire is within 1000 feet of elevation from the the location where weather was recorded,
 *  A: The fire is 1000 - 2000 feet above the location where weather was recorded,
 *  This is for field use only and may be omitted for other applications.
 * @param units: The units to use.  Only relevant to temp.
 *
 * @returns 1-hour fuel moisture (fraction: water weight/dry fuel weight).
 */
double FosbergNWCG_1HrFM(std::string tableA_Path, std::string tableB_Path, std::string tableC_Path,
                         std::string tableD_Path, double temp, double rh, int monthOfYear,
                         int hourOfDay, double slopePct, char aspectCardinal, bool shaded,
                         char elevation, UnitsType units)
{
	std::string correctionTablePath;
	double tempF, rfm, correction;

	//Convert the temperature if needed:
	if (units == US)
	{
		tempF = temp;
	}
	else//units == Metric
	{
		tempF = CtoF(temp);
	}

	rfm = FosbergNWCG_GetRFM(tableA_Path, tempF, rh);//Look up the reference fuel moisture.
	
	if (monthOfYear >= 5 && monthOfYear <= 7)//May - July
	{
		correctionTablePath = tableB_Path;
	}
	else if ((monthOfYear >= 2 && monthOfYear <= 4) || (monthOfYear >= 8 && monthOfYear <= 10))//Feb - April, Aug - Oct
	{
		correctionTablePath = tableC_Path;
	}
	else if (monthOfYear == 1 || monthOfYear == 11 || monthOfYear == 12)//Nov - Jan
	{
		correctionTablePath = tableD_Path;
	}
	else
	{
		Stop("Invalid month value. Please supply the integer month of year.");
	}
	//Or:
// 	switch (monthOfYear) {
// 		//May - July:
// 		case 5:
// 		case 6:
// 		case 7:
// 			correctionTablePath = tableB_Path;
// 			break;
// 		
// 		//Feb - April, Aug - Oct:
// 		case 2:
// 		case 3:
// 		case 4:
// 		case 8:
// 		case 9:
// 		case 10:
// 			correctionTablePath = tableC_Path;
// 			break;
// 		
// 		//Nov - Jan:
// 		case 1:
// 		case 11:
// 		case 12:
// 			correctionTablePath = tableD_Path;
// 			break;
// 		
// 		default:
// 			Stop("Invalid month value. Please supply the integer month of year.");
// 			break;
// 	}

	//Look up the correction factor for the conditions specified:
	correction = FosbergNWCG_GetCorrection(correctionTablePath, hourOfDay, slopePct, aspectCardinal,
	                                       shaded, elevation);

	return rfm + correction;//Combine and return.
}

/** Look up the Fosberg reference fuel moisture (NWCG variant) given the temperature and humidity:
 *
 * @param tableA_Path = A path to a tab delimited file holding lookup table A from Rothermel 1981 & the NWCG.
 *
 * @param temp The air temperature at 4.5 feet / 1.37 meters above ground level (degrees F / C).
 *             The results is undefined below 10 degrees Fahrenheit!!!!!
 * @param rh Relative humidity (percent).
 *
 * @returns Reference fuel moisture (fraction: water weight/dry fuel weight).
 */
double FosbergNWCG_GetRFM(std::string tableA_Path, double tempF, double rh)
{
	/*The size and shape of table A is known and fixed.  However, it could be better not to assume
	the shape.  The RH bins can be easily counted since they are on the first line.  The number
	of temperature bins on the other hand would require counting rows, which causes a bit more
	complication.*/
	const int numRHBins = 21;
	const int numTempBins = 6;
	
	//The temperatures and RHs are all integer but we could treat them as doubles?
	int rhRangeBottoms[numRHBins];
	int tempRangeBottoms[numTempBins];
	int tempIndex;//Used first as an counter and then as the row match index.
	int rhIndex;//The matching column index.

	std::string line;//To hold lines of the input file.

	int luMatrix[numTempBins][numRHBins];

	//Parameter checking:
	//I need to determine a valid temp range.  The tables have a bottom but not a top.
	if (!InRange(rh, 0, 100))
	{
		Stop("Relative humidity must be a valid percentage");
	}

	//Open the file:
	std::ifstream tableAFile(tableA_Path);
	//Error handling?????

	//Break into IDs and data:
	
	//The first line is the relative humidity row ID header with leading whitespace:
	std::getline(tableAFile, line);
	std::stringstream lineStr(line);//This line as a stream.

	//The first field will be empty whitespace and will be ignored, the rest are the RH bounds:
	for (int i = 0; i < numRHBins; i++)
	{
		lineStr >> rhRangeBottoms[i];
	}

	//Process he following lines:
	tempIndex = 0;
	while (std::getline(tableAFile, line))
	{
		lineStr.clear();//Needed?
		lineStr.str(line);

		//The first field is temperature row IDs:
		lineStr >> tempRangeBottoms[tempIndex];

		//The rest of the line is lookup table values:
		for (int j = 0; j < numRHBins; j++)
		{
			lineStr >> luMatrix[tempIndex][j];
		}

		tempIndex += 1;
	}
	tempIndex = -1;//Reset to use as the search index.

	//Search for the matching value:
	//numTempBins = tempRangeBottoms.size();
	for (int i = 0; i < numTempBins; i++)
	{
		if (i < (numTempBins - 1))
		{
			//The top bound for each bin is the bottom of the next:
			if (tempF >= tempRangeBottoms[i] && tempF < tempRangeBottoms[i + 1])
			{
				tempIndex = i;
				break;
			}
		}
		else
		{
			//The last bin only has a lower temperature bound:
			if (tempF >= tempRangeBottoms[i])
			{
				tempIndex = i;
			}
			else//This could happen if the temperature is below the lowest bin.  Add check!!!!
			{
				Stop("Match not found.");
			}
		}
	}

	//numRHBins = rhRangeBottoms.size();
	for (int j; j < j < numRHBins; j++)
	{
		if (j < numRHBins)
		{
			//The top bound for each bin is the bottom of the next:
			if (rh >= rhRangeBottoms[j] && rh < rhRangeBottoms[j + 1])
			{
				rhIndex = j;
				break;
			}
		}
		else
		{
			//The value must match the final bin but we check anyway:
			if (rh == rhRangeBottoms[j])
			{
				rhIndex = j;
			}
			//else//Should not possible with a valid values.
			//{
			//	Stop("Match not found.");
			//}
		}
	}

	return luMatrix[tempIndex][rhIndex];//Integer!!!!!
}

/** Look up the Fosberg fuel moisture correction (NWCG variant) for the conditions passed:
 * Note: This doesn't handle the nighttime values from table C or D.
 *
 * @param tableFilePath A path to a tab delimited file holding a correction lookup table (B - D)
 *                      from Rothermel 1981 & the NWCG.
 *
 * @param hourOfDay The numeric hour of day in 24 hour format (1 - 24).
 * @param slopePct Percent slope of the location (rise / run x 100%).
 * @param aspectCardinal The slope aspect as a cardinal direction (N, E, S, W).
 *   Note: We could add the ability to take the aspect as degrees.
 * @param shaded Does the location have 50% or greater canopy cover or is there full cloud cover.
 * @param elevation A code indicating the slope position for the prediction relative to where the
 *                  weather conditions were taken.  Values are:
 *  B: The fire is 1000 - 2000 feet below the location where weather was recorded,
 *  L: The fire is within 1000 feet of elevation from the the location where weather was recorded,
 *  A: The fire is 1000 - 2000 feet above the location where weather was recorded,
 *  This is for field use only and may be omitted for other applications.
 *
 * @returns Fuel moisture correction (fraction: water weight/dry fuel weight).
 */
double FosbergNWCG_GetCorrection(std::string tableFilePath, int hourOfDay, double slopePct,
                                 char aspectCardinal, bool shaded, char elevation)
{
	/*The size and shape of tables B-D are known, fixed, and the same for all currently.  However,
	it could be better not to assume the shape, especially since D could be modified with more
	columns for nighttime hours.*/
	const int numCols = 18;
	const int numRows = 12;

	//Column IDs:
	int hourStart[numCols];
	int hourEnd[numCols];
	char elevationCode[numCols];

	//Row IDs:
	std::string shadeID[numRows];
	char aspects[numRows];
	std::string slopeClass[numRows];

	//int theRow, theCol;//The matching row and column indexes.
	//The matching row and column indexes:
	int theRow = -1;
	int theCol = -1;

	std::string line;//To hold lines of the input file.

	int luMatrix[numRows][numCols];

	//Open the file:
	std::ifstream tableFile(tableFilePath);
	//Error handling?????

	//Break into IDs and data.  There are 3 column IDs and 3 row IDs:

	//Column IDs:
	//The column IDs give start and end times for each column in military time hours and minutes
	//without a colon.  Convert to hour of day:

	std::getline(tableFile, line);
	std::stringstream lineStr(line);//This line as a stream.

	//The first field will be empty whitespace and will be ignored, the rest are the hour bounds:
	for (int i = 0; i < numCols; i++)
	{
		lineStr >> hourStart[i];
		hourStart[i] = hourStart[i] / 100;
	}

	std::getline(tableFile, line);
	lineStr.clear();
	lineStr.str(line);
	for (int i = 0; i < numCols; i++)
	{
		lineStr >> hourEnd[i];
		hourEnd[i] = hourEnd[i] / 100;
	}

	//Parameter checking (after extracting hours):
	//if (!InRange(hourOfDay, min(hourStart), max(hourEnd)))
	if (!InRange(hourOfDay, *std::min_element(hourStart, hourStart + numCols),
	             *std::max_element(hourEnd, hourEnd + numCols)))
	{
		Stop("Invalid hour of day value or nighttime value passed.");
	}
	//The problem with percent slope is that is becomes huge as it approaches 90 degrees.  However,
	//realistically very high slopes should be very rare.  This check is minimal.
	if (!InRange(slopePct, 0, 1000))//0 - ~85 degrees.
	{
		Stop("Invalid percent slope");
	}
	if (!(aspectCardinal == 'N' || aspectCardinal == 'S' ||
	      aspectCardinal == 'E' || aspectCardinal == 'W'))
	{
		Stop("Invalid aspect.  Must be a cardinal direction.");
	}
	//Shaded should be logical but we could also accept a percentage.
	if (!(elevation == 'B' || elevation == 'L' || elevation == 'A'))
	{
		Stop("Invalid relative elevation code.");
	}

	for (int i = 0; i < numRows; i++)
	{
		lineStr >> elevationCode[i];
	}

	//Process he following lines:
	int counter = 0;
	while (std::getline(tableFile, line))
	{
		lineStr.clear();
		lineStr.str(line);

		//Row IDs:
		lineStr >> shadeID[counter] >> aspects[counter] >> slopeClass[counter];

		//The rest of the line is lookup table values:
		for (int j = 0; j < numRows; j++)
		{
			lineStr >> luMatrix[counter][j];
		}

		counter += 1;
	}

	//Find the matching row by ID:
	std::string theShadeID;
	if (shaded)
	{
		theShadeID = "Shaded";
	}
	else
	{
		theShadeID = "Unshaded";
	}
	
	std::string theSlopeID;
	if (slopePct <= 30)
	{
	  theSlopeID = "0-30%";
	}
	else
	{
	  theSlopeID = ">30%";
	}
	
	for (int j = 0; j < numRows; j++)
	{
		if (shadeID[j] == theShadeID)
		{
			if (slopeClass[j] == theSlopeID)
			{
				if (aspects[j] == aspectCardinal)
				{
					theRow = j;
					break;
				}
			}
		}
	}

	if (theRow = -1)
	{
		Stop("Failed to find a matching row.");
	}

	//Find the matching column:
	for (int j = 0; j < numRows; j += 3)//The hours ranges are repeated by threes (B, L, A).
	{
		if (hourOfDay >= hourStart[j] && hourOfDay < hourEnd[j])
		{
			for (int k = j; k < (j + 2); k++)
			{
				if (elevation == elevationCode[j])
				{
					theCol = k;
					break;
				}
			}
		}
	}

	if (theCol = -1)
	{
		Stop("Failed to find a matching column.");
	}

	return luMatrix[theRow][theCol];//Integer!!!!!
}

/* Make an estimate of the 100-hour dead fuel moisture:
 *
 * The estimation is based on the 1-hour Fosberg fuel moisture.  Rothermel 1981 suggests adding 1%.
 * NWCG suggests adding ~1-2%.  We use the middle of the latter recommendation.
 *
 * @param oneHrFM Fosberg NWCG prediction of 1-hr fuel moisture (fraction: water weight/dry fuel weight).
 *
 * @returns: 10-hour fuel moisture (fraction: water weight/dry fuel weight).
 */
double NWCG_10hrFM(double oneHrFM)
{
  return oneHrFM * 1.015;//Add 1.5%.
}

/* Make an estimate of the 100-hour dead fuel moisture:
 *
 * The estimation is based on the 1-hour Fosberg fuel moisture.  Rothermel 1981 suggests adding 2%.
 * WCG suggests adding ~2-4%.  We use the middle of the latter recommendation.
 *
 * @param oneHrFM Fosberg NWCG prediction of 1-hr fuel moisture (fraction: water weight/dry fuel weight).
 *
 * @eturns: 100-hour fuel moisture (fraction: water weight/dry fuel weight).
 */
double NWCG_100hrFM( double oneHrFM)
{
  return(oneHrFM * 1.03);//Add 3%.
}

//Unit Tests:---------------------------------------------------------------------------------------

/** Test FosbergNWCG_1HrFM() using the set of daytime scenarios outlined on Rothermel 1983 pg. 20:
 *
 * This is not an ideal unit test since it requires external file specificaitons.
 *
 * @param tableA_Path - tableD_Path = Paths to tab delimited files holding the lookup table data.
 *   These tables are used in Rothermel 1981 and are also the NWCG website and the 2022 IPTG.
 *
 * @returns True if the test passes, false otherwise.
 */
bool FosbergNWCG_1HrFM_UnitTest(std::string tableA_Path, std::string tableB_Path,
                                std::string tableC_Path, std::string tableD_Path)
{
	bool pass = true;
	//print("Daytime fuel moisture scenarios (Rothermel 1983, pg. 20):")//Verbose?????
	
	//Should be 6%:
	double valA = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
	                                80, 20, 8, 13, 35, 'N');
	if (valA != 6)
	{
		std::string msg = "a. Expected value = 6%, calculated: " + std::to_string(valA) + "%";
		Warning(msg);
		pass = false;
	}
	
	//Should be 4%:
	double valB = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
                                    80, 20, 8, 13, 35, 'E', false, 'B');
	if (valB != 4)
	{
		std::string msg = "b. Expected value = 4%, calculated: " + std::to_string(valB) + "%";
		Warning(msg);
		pass = false;
	}
	
	//Should be 7%:
	double valC = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
	                                80, 20, 8, 13, 35, 'N', true);
	if (valC != 7)
	{
		std::string msg = "c. Expected value = 7%, calculated: " + std::to_string(valC) + "%";
		Warning(msg);
		pass = false;
	}
	
	return pass;
}

