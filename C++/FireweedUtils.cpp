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

/** This utility checks that the parameters (vectors) passed have the same length.  2-4 and 7
 * arguments are accepted.
 *
 * Overloading has been used to reproduce the behavior of the R function, although that function can
 * handle arguments of different types in any order.  We only handle a subset of possible type
 * combinations.  This function is currently used primarily for variable vectors of type double and
 * liveDead, which is currently an integer but could be a boolean.  To reduce the number of
 * possibilities we require liveDead to be last.  It might be possible to do this more compactly
 * with template functions or ariadic arguments.
 *
 * @param arg1, arg2, ... A series of vector argument to compare.
 *
 * @returns Do the vectors have the same number of elements.
 *
 * Could change arguments to const &?
 */
// bool SameLengths(std::vector<double> arg1, std::vector<double> arg2)
// {
// 	return (arg1.size() == arg2.size());
// }

// template <typename T1, typename T2>
// bool SameLengthsT(const std::vector<T1>& arg1, const std::vector<T2>& arg2)
// {
// 	return (arg1.size() == arg2.size());
// }

// bool SameLengths(std::vector<double> arg1, std::vector<int> arg2)
// {
//  	return (arg1.size() == arg2.size());
// }

// bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3)
// {
// 	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3));
// }

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3));
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3) && SameLengths(arg1, arg4));
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3) && SameLengths(arg1, arg4));
}

//7 doubles:
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4, std::vector<double> arg5, std::vector<double> arg6,
                 std::vector<double> arg7)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3) && SameLengths(arg1, arg4) &&
			SameLengths(arg1, arg5) && SameLengths(arg1, arg6) && SameLengths(arg1, arg7));
}

//7 doubles and an int vector:
bool SameLengths(std::vector<double> dbl1, std::vector<double> dbl2, std::vector<double> dbl3,
                 std::vector<double> dbl4, std::vector<double> dbl5, std::vector<double> dbl6,
                 std::vector<double> dbl7, std::vector<int> int8)
{
	return (SameLengths(dbl1, dbl2) && SameLengths(dbl1, dbl3) && SameLengths(dbl1, dbl4) &&
			SameLengths(dbl1, dbl5) && SameLengths(dbl1, dbl6) && SameLengths(dbl1, dbl7) &&
			SameLengths(dbl1, int8));
}

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
