/***************************************************************************************************
FireweedUnits.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 9/16/2024
Reference: Proj. 11 Exp. 19

	This file is part of the Fireweed wildfire code library.
	This header file defines unit conversion functions used in the library.

***************************************************************************************************/

#include "FireweedUnits.h"

/** Convert a temperature in degrees Celsius to degrees Fahrenheit.
 *
 * @param degreesC The temperature to convert (degrees Celsius).
 *
 * @returns	Temperature in degrees Fahrenheit.
 */
double CtoF(double degreesC)
{
  return degreesC * 9/5 + 32;
}
