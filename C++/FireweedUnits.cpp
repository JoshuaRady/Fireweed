/***************************************************************************************************
FireweedUnits.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 19, 20, 18

	This file is part of the Fireweed wildfire code library.
	This header file defines unit conversion functions used in the library.

***************************************************************************************************/

#include "FireweedUnits.h"
#include "FireweedMessaging.h"
#include "FireweedUtils.h"

/** Convert a temperature in degrees Celsius to degrees Fahrenheit.
 *
 * @param degreesC The temperature to convert (degrees Celsius).
 *
 * @returns	Temperature in degrees Fahrenheit.
 */
double CtoF(double degreesC)
{
  return (degreesC * 9/5) + 32;
}

/** Convert a temperature in degrees Fahrenheit to degrees Celsius.
 *
 * @param degreesF The temperature to convert (degrees Fahrenheit).
 *
 * @returns Temperature in degrees Celsius.
 */
double FtoC(double degreesF)
{
  return (degreesF - 32) * 5/9;
}

/** Convert from percent slope to slope steepness as used by the Rothermel Albini spread model.
 * 
 * @param slopePct The slope or grade etc. of the landscape in percent.
 * 
 * @returns Slope steepness (unitless fraction: vertical rise / horizontal distance), AKA tan ϕ.
 * Slope values are forced to postive values.
 */
double SlopePctToSteepness(double slopePct)
{
	//The percent slope is simply the rise / run x 100%:
	double slopeSteepness = std::abs(slopePct) / 100;
	
	CheckSlope(slopeSteepness);//Check that the value is reasonable.
	
	return slopeSteepness;
}

/** Convert from slope in degrees to slope steepness as used by the Rothermel Albini spread model.
 * 
 * @param slopeDegrees The slope or grade etc. of the landscape in degrees.
 * 
 * @returns Slope steepness (unitless fraction: vertical rise / horizontal distance), AKA tan ϕ.
 * Slope values are forced to postive values.
 */
double SlopeDegreesToSteepness(double slopeDegrees)
{
	//While any degree value may be mathematically valid values beyond 90 are not really used to
	//express slopes.  Large values may indicate the input value was miscomputed.  Negative slopes
	//could be valid too but are questionable.  Allow but warn:
	if (!InRange(slopeDegrees, 0, 90))
	{
		Warning("Slope is outside the expected range of 0 - 90 degrees.");
	}

	//tan() takes degrees:
	double slopeSteepness = std::tan(std::abs(slopeDegrees) * (Pi/180));

	CheckSlope(slopeSteepness);//Check that the value is reasonable.

	return slopeSteepness;
}

/** Check if a slope is reasonable:
 * 
 * The problem with fractional or percent slope is that is becomes huge as it approaches 90 degrees.
 * However, realistically very high slopes should be very rare on the landscape.
 * This check is minimal and might be improved.
 *
 * @param slopeSteepness Slope steepness (unitless fraction: vertical rise / horizontal distance).
 *   AKA tan ϕ.
 *
 * @returns Nothing.  Posts a warning.  Could return a status instead.  ValidSlope()?
 */
void CheckSlope(double slopeSteepness)
{
	if (!InRange(slopeSteepness, 0, 11.43005))//0 - ~85 degrees.
	{
		Warning("Questionable slope value.");
	}
}
