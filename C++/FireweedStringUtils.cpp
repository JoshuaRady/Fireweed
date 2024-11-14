/***************************************************************************************************
FireweedStringUtils.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 8/29/2024
Reference: Proj. 11 Exp. 20

	This file is part of the Fireweed wildfire code library.  It defines a set of string utilities.

***************************************************************************************************/

#include <algorithm>//For std::all_of().
#include <iostream>
#include <sstream>

#include "FireweedStringUtils.h"
#include "FireweedMessaging.h"

/** Split a delimited string into a vector of substrings.
 *
 * @param str A string containing delimited text fields to split.
 * @param delimiter The delimiter character.
 *
 * This should probably be moved to a utility file somewhere.
 */
std::vector<std::string> SplitDelim(const std::string& str, char delimiter)
{
	std::stringstream strStrm(str);//The string as a stream.
	std::string substring;//To hold extracted substrings (or tokens).
	std::vector<std::string> substrings;//Return value: A vector to hold the split string.
	
	while(getline(strStrm, substring, delimiter))
	{
		substrings.push_back(substring);
	}
	
	return substrings;
}

/** Split a delimited string into a vector of substrings.
 *
 * @param str A string containing delimited text fields to split.
 * @param delimiter The delimiter character.
 * @param stripQuotes If true allow double quoted fields.  Delimiter characters inside quoted fields
 *                    will ignored.  Quotes will be stripped from the output.
 *
 * This should probably be moved to a utility file.
 */
std::vector<std::string> SplitDelim(const std::string& str, char delimiter, bool allowQuotes)
{
	if (!allowQuotes)
	{
		 return SplitDelim(str, delimiter);
	}
	else
	{
		char quote = '\"';
		std::stringstream strStrm(str);//The string as a stream.
		bool content;//Was content extracted from the stream.
		std::string substring;//To hold extracted substrings (or tokens).
		std::vector<std::string> substrings;//Return value: A vector to hold the split string.

		do//Split the string while there is content to parse:
		{
			if (strStrm.peek() == quote)
			{
				//Discard the quote:
				char discard = strStrm.get();
				
				//Look for the next quote:
				if (getline(strStrm, substring, quote))
				{
					content = true;
					substrings.push_back(substring);
				}
				else
				{
					//If the closing quote is missing we can't determine where the field ends.
					//We could ignore the opening quote and extract or discard to the next delimiter.
					//These are all risky.  It is better to throw an error. !!!!!
					Stop("Closing quote not found.");
				}

				//Discard the delimiter that follows:
				std::string trash;
				if (getline(strStrm, trash, delimiter))
				{
					//Check trash to make sure it is empty:
					//The next character should be a delimiter, whitespace followed by a delimiter
					//(tolerable but bad form), or nothing if we are at the end of the string.
					//Anything else means that we could be discarding data.
					if (!std::all_of(trash.begin(), trash.end(), isspace))
					{
						//Warn or throw error: !!!!!
						Warning("Unexpected content when splitting.");
					}
				}
			}
			else//Look for the next delimiter:
			{
				if (getline(strStrm, substring, delimiter))
				{
					content = true;
					substrings.push_back(substring);
				}
				else
				{
					content = false;
				}
			}
		} while (content);

		return substrings;
	}
}

/** Print a numeric vector to an output stream on a single line with separators.
 *
 * @param output The output stream to print to.
 * @param vec The string vector to print.
 * @param separator The string to separate vector elements.  Default to a comma and space.
 *
 * Turn into a template?
 */
std::ostream& PrintVector(std::ostream& output, const std::vector <double>& vec, std::string separator)
{
	for (int i = 0; i < vec.size() - 1; i++)
	{
		output << vec[i] << separator;
	}
	output << vec[vec.size() - 1] << std::endl;

	return output;
}

/** Convert a numeric vector to a string with separators between elements and return.
 *
 * @param vec The string vector to print.
 * @param separator The string to separate vector elements.  Default to a comma and space.
 *
 * Turn into a template?
 */
std::string VectorToStr(const std::vector <double> vec, std::string separator)
{
	std::string str;

	for (int i = 0; i < vec.size() - 1; i++)
	{
		str = str + std::to_string(vec[i]) + separator;
	}
	str = str + std::to_string(vec[vec.size() - 1]);
	
	return str;
}
