/***************************************************************************************************
FireweedUtils.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 8/12/2024
Reference: Proj. 11 Exp. 19

	This file is part of the Fireweed wildfire code library.
	This file defines a set of (general) utilities.

***************************************************************************************************/

#include <cmath>//fabs()

#include "FireweedUtils.h"

//See the header for SameLengths().

//This utility checks that a value falls in a valid range.
bool InRange(double value, double low, double high)
{
	return (value >= low && value <= high);
}

bool InRange(std::vector<double> value, double low, double high)
{
	for (int i = 0; i < value.size(); i++)
	{
		if (value[i] < low || value[i] > high)
		{
			return false;
		}
	}

	return true;
}

/** Check if a value is from 0 to 1, a valid proportion.
 *
 * @param value A number to check.
 *
 * @returns Is the value from 0 to 1.
 */
bool ValidProportion(double value)
{
	return InRange(value, 0, 1);
}

/** Check if a value is from 0 to 1, a valid proportion.
 *
 * @param value A vector of numbers to check.
 *
 * @returns Are all the values from 0 to 1.
 */
bool ValidProportion(std::vector<double> value)
{
	return InRange(value, 0, 1);
}

/** Compare two floating point values for near/effective equality:
 *
 * It is easy for math operations to result in small floating point errors that may make two
 * variables you expect to be identical to be slightly different.  Use this fucntion in place of ==
 * in such cases.
 *
 * @param val1 The first value to compare.
 * @param val2 The second value to compare.
 * @param precision How close must the values be to be "equal".
 *
 * @returns True if the two values are nearly equal.
 *
 * @note I'm not sure what the default should be for the precision of this comparison (see header for default value).
 * @note This function is currently C++ only.
 */
bool FloatCompare(double val1, double val2, double precision)//Or epsilon?
{
	if (std::fabs(val1 - val2) < precision)
	{
		return true;
	}
	else
	{
		return false;
	}
}
