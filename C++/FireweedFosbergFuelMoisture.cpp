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

#include "FireweedFosbergFuelMoisture.h"

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
                         std::stringtableD_Path, double temp, double rh, int monthOfYear,
                         int hourOfDay, double slopePct, char aspectCardinal, bool shaded,
                         char elevation, UnitsType units)
{
	string correctionTablePath;
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
