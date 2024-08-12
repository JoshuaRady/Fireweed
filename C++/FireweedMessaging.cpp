/***************************************************************************************************
FireweedMessaging.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  It included a set of functions for
simple error messaging and logging.

***************************************************************************************************/

#include <iostream>
#include "FireweedMessaging.h"

/** Log a neutral message.
 *
 * @param message A message to log.
 */
void LogMsg(const char* message)
{
	std::cout << message << "\n";
}

/** Log a neutral message with a numeric value.
 *
 * @param message A message to log.
 * @param value A numeric value to appended after the message.  A space is added between them.
 */
void LogMsg(const char* message, double value)
{
	std::cout << message << " " << value << "\n";
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void LogMsg(const char* message, std::vector<double> value)
{
	std::cout << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		std::cout << value[i] << ", ";
	}

	std::cout << value[value.size() - 1] << "\n";
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void LogMsg(const char* message, std::vector<int> value)
{
	std::cout << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		std::cout << value[i] << ", ";
	}

	std::cout << value[value.size() - 1] << "\n";
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void Warning(const char* message)
{
	std::cout << message << "\n";
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void Warning(const std::string& message)
{
	std::cout << message << "\n";
}

/** Post the passed message and shutdown (not yet implemented!!!!!).
 *
 * @param message An error message.
 */
void Stop(const char* message)
{
	std::cout << message << "\n";
	//Add error throwing.
}

/** Post the passed message and shutdown (not yet implemented!!!!!).
 *
 * @param message An error message.
 */
void Stop(const std::string& message)
{
	std::cout << message << "\n";
	//Add error throwing.
}