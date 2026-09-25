/***************************************************************************************************
FireweedStringUtils.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 8/29/2024
Reference: Proj. 11 Exp. 20

	This file is part of the Fireweed wildfire code library.  This header file declares a set of
string utilities.

***************************************************************************************************/
#ifndef FIREWEEDSTRINGUTILS_H
#define FIREWEEDSTRINGUTILS_H

#include <string>
#include <type_traits>
#include <vector>

std::vector<std::string> SplitDelim(const std::string& str, char delimiter);
std::vector<std::string> SplitDelim(const std::string& str, char delimiter, bool allowQuotes);
std::ostream& PrintVector(std::ostream& output, const std::vector <double>& vec, std::string separator = ", ");

/** Convert a numeric or boolean vector to a string with separators between elements and return it.
 *
 * @param vec The string vector to print.
 * @param separator The string to separate vector elements.  Defaults to a comma and space.
 */
template <typename T>
std::string VectorToStr(const std::vector<T>& vec, std::string separator = ", ")
{
	static_assert(std::is_arithmetic_v<T>, "VectorToStr() expects numeric or boolean vectors.");

	std::string str;

	for (int i = 0; i < vec.size() - 1; i++)
	{
		str += std::to_string(vec[i]) + separator;
	}
	str += std::to_string(vec[vec.size() - 1]);
	
	return str;
}

#endif //FIREWEEDSTRINGUTILS_H
