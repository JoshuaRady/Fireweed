/***************************************************************************************************
FireweedUnits.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 19, 20

	This file is part of the Fireweed wildfire code library.
	This header file defines unit constants and conversion factors and declares unit functions used
in the library.

***************************************************************************************************/
#ifndef FIREWEEDUNITS_H
#define FIREWEEDUNITS_H

#include <cmath>//Include in the header for nonstandard M_PI.

/** UnitsType
 *
 * A type to specify the units class to be used or in use.  The options are 'US' for United States
 * customary units or 'Metric' for metric.
 */
enum UnitsType {US, Metric};

/** UnitCode
 *
 * Symbols for identifying specific units.
 *
 * Should there be undefined / no units / NA symbol?
 */
enum UnitCode {lbPer_ft2, tonPerAc, kgPer_m2, Percent, Fraction};

//--------------------------------------------------------------------------------------------------
/** Unit Conversion Factors
 *
 */

//Length: (exact per international yard and pound act)
const double cmPerIn = 2.54;
const double cmPerFt = 30.48;
const double mPerFt = 0.3048;
const double ftPerM = 3.28084;//1 / mPerFt
const double ftPerMi = 5280;//* for conversion of windspeed (U, MPH * ftPerMi / 60 = ft/min)

//SAV is in ft^2/ft^3 = 1/ft or cm^2/cm^3 = 1/cm
//Therefore units convert: ft^2/ft^3 * cmPerFt^2/cmPerFt^2 = 1/ft * 1/cmPerFt = 1/cm
//So: SAVft * 1/cmPerFt = SAVft / cmPerFt = SAVcm

//Area:
const double ft2PerAcre = 43560;

//Mass:
const double kgPerLb = 0.453592;
const int lbsPerTon = 2000;

//Density:
const double lbPerFtCuToKgPerMCu = 16.0185;//kgPerLb * (ftPerM)^3, 16.01846337396

//JPerBtu = 1055.06 or 1,054.35
/*The definition of a BTU can vary resulting in several different conversion factors.  Wilson 1980
seems to have used a value close to the themochemical value of 1.05435 J/BTU, based on his heat of
preignition conversion.  We will use that to be consistent with his converted constant values.
The IT value of 1.05506 would be a reasonable alternative.*/
const double kJPerBtu = 1.05435;

//tons/ac -> lb/ft^2: (See fuel loading note in FireweedFuelModels.cpp.)
const double tonsPerAcToLbPerSqFt = lbsPerTon / ft2PerAcre;

//pi:
#ifdef M_PI
	const double Pi = M_PI;
#else
	const double Pi = 3.14159265358979323846;
#endif 

//--------------------------------------------------------------------------------------------------
/** Units Functions
 *
 */
 
double CtoF(double degreesC);
double FtoC(double degreesF);

double SlopePctToSteepness(double slopePct);
double SlopeDegreesToSteepness(double slopeDegrees);
void CheckSlope(double slopeSteepness);

#endif //FIREWEEDUNITS_H
