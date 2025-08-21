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

/** The SameLengths() utility checks that the vectors passed have the same length. 2-12 vectors
 * are accepted.
 *
 * Overloading and templates have been used to handle vector vecuments of different types in any
 * order reproducing the behavior of the R function.
 *
 * @param vec1, vec2, ... A series of vector vecument to compare.
 *
 * @returns Do the vectors have the same number of elements?
 */
template <typename T1, typename T2>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2)
{
	return (vec1.size() == vec2.size());
}

template <typename T1, typename T2, typename T3>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2, const std::vector<T3>& vec3)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3));
}

template <typename T1, typename T2, typename T3, typename T4>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7));
}

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

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7, const std::vector<T8>& vec8,
                 const std::vector<T9>& vec9)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7) &&
			SameLengths(vec1, vec8) && SameLengths(vec1, vec9));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9, typename T10>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7, const std::vector<T8>& vec8,
                 const std::vector<T9>& vec9, const std::vector<T10>& vec10)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7) &&
			SameLengths(vec1, vec8) && SameLengths(vec1, vec9) && SameLengths(vec1, vec10));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9, typename T10, typename T11>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7, const std::vector<T8>& vec8,
                 const std::vector<T9>& vec9, const std::vector<T10>& vec10,
                 const std::vector<T11>& vec11)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7) &&
			SameLengths(vec1, vec8) && SameLengths(vec1, vec9) && SameLengths(vec1, vec10) &&
			SameLengths(vec1, vec11));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7,
          typename T8, typename T9, typename T10, typename T11, typename T12>
bool SameLengths(const std::vector<T1>& vec1, const std::vector<T2>& vec2,
                 const std::vector<T3>& vec3, const std::vector<T4>& vec4,
                 const std::vector<T5>& vec5, const std::vector<T6>& vec6,
                 const std::vector<T7>& vec7, const std::vector<T8>& vec8,
                 const std::vector<T9>& vec9, const std::vector<T10>& vec10,
                 const std::vector<T11>& vec11, const std::vector<T12>& vec12)
{
	return (SameLengths(vec1, vec2) && SameLengths(vec1, vec3) && SameLengths(vec1, vec4) &&
			SameLengths(vec1, vec5) && SameLengths(vec1, vec6) && SameLengths(vec1, vec7) &&
			SameLengths(vec1, vec8) && SameLengths(vec1, vec9) && SameLengths(vec1, vec10) &&
			SameLengths(vec1, vec11) && SameLengths(vec1, vec12));
}

bool InRange(double value, double low, double high);
bool InRange(std::vector<double> value, double low, double high);

bool ValidProportion(double value);
bool ValidProportion(std::vector<double> value);

bool FloatCompare(double val1, double val2, double precision = 0.0001);

#endif
