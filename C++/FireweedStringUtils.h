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
#include <vector>

std::vector<std::string> SplitDelim(const std::string& str, char delimiter);
std::vector<std::string> SplitDelim(const std::string& str, char delimiter, bool allowQuotes);
std::ostream& PrintVector(std::ostream& output, const std::vector <double>& vec, std::string separator = ", ");
std::string VectorToStr(const std::vector <double> vec, std::string separator = ", ");

#endif //FIREWEEDSTRINGUTILS_H
