/***************************************************************************************************
FireweedUtils.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 8/12/2024
Reference: Proj. 11 Exp. 19

	This file is part of the Fireweed wildfire code library.
	This header file declares a set of (general) utilities.

***************************************************************************************************/
#ifndef FIREWEEDUTILS_H
#define FIREWEEDUTILS_H

#include <vector>

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<int> arg2);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4);
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4, std::vector<double> arg5, std::vector<double> arg6,
                 std::vector<double> arg7);
bool SameLengths(std::vector<double> dbl1, std::vector<double> dbl2, std::vector<double> dbl3,
                 std::vector<double> dbl4, std::vector<double> dbl5, std::vector<double> dbl6,
                 std::vector<double> dbl7, std::vector<int> int8);

bool InRange(double value, double low, double high);
bool InRange(std::vector<double> value, double low, double high);

bool ValidProportion(double value);
bool ValidProportion(std::vector<double> value);

bool FloatCompare(double val1, double val2, double precision = 0.0001);

#endif
