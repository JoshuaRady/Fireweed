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

//bool SameLengths(std::vector<double> arg1, std::vector<double> arg2);

template <typename T1, typename T2>
bool SameLengths(const std::vector<T1>& arg1, const std::vector<T2>& arg2)
{
	return (arg1.size() == arg2.size());
}

//bool SameLengths(std::vector<double> arg1, std::vector<int> arg2);
//bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3);
template <typename T1, typename T2, typename T3>
bool SameLengths(const std::vector<T1>& arg1, const std::vector<T2>& arg2, const std::vector<T3>& arg3)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3));
}

//bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3);

// bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
//                  std::vector<double> arg4);
// bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
//                  std::vector<int> arg4);
template <typename T1, typename T2, typename T3, typename T4>
bool SameLengths(const std::vector<T1>& arg1, const std::vector<T2>& arg2,
                 const std::vector<T3>& arg3, const std::vector<T4>& arg4)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3) && SameLengths(arg1, arg4));
}

// bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
//                  std::vector<double> arg4, std::vector<double> arg5, std::vector<double> arg6,
//                  std::vector<double> arg7);
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
bool SameLengths(const std::vector<T1>& arg1, const std::vector<T2>& arg2,
                 const std::vector<T3>& arg3, const std::vector<T4>& arg4,
                 const std::vector<T5>& arg5, const std::vector<T6>& arg6,
                 const std::vector<T7>& arg7)
{
	return (SameLengths(arg1, arg2) && SameLengths(arg1, arg3) && SameLengths(arg1, arg4) &&
			SameLengths(arg1, arg5) && SameLengths(arg1, arg6) && SameLengths(arg1, arg7));
}

// bool SameLengths(std::vector<double> dbl1, std::vector<double> dbl2, std::vector<double> dbl3,
//                  std::vector<double> dbl4, std::vector<double> dbl5, std::vector<double> dbl6,
//                  std::vector<double> dbl7, std::vector<int> int8);
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7, const std::vector<T8>& vec8)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7) &&
			SameLengths(vec1, vec8));
}

bool InRange(double value, double low, double high);
bool InRange(std::vector<double> value, double low, double high);

bool ValidProportion(double value);
bool ValidProportion(std::vector<double> value);

bool FloatCompare(double val1, double val2, double precision = 0.0001);

#endif
