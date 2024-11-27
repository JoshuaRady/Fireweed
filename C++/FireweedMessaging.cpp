/***************************************************************************************************
FireweedMessaging.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  It includes a set of functions for
simple error messaging and logging.

***************************************************************************************************/

#include <exception>
//#include <iostream>
#include "FireweedMessaging.h"

FWMessenger Msg;

/** Constructor
 *
 * Default all streams to cout.  These can be customized later.
 */
FWMessenger::FWMessenger()
{
	logStream = &std::cout;
	warnStream = &std::cout;
	errorStream = &std::cout;
}

/** Set the stream used to post log messages.
 *
 */
void FWMessenger::SetLogStream(std::ostream* streamPtr)
{
	logStream = streamPtr;
}

/** Set the stream used to post warning messages.
 *
 */
void FWMessenger::SetWarnStream(std::ostream* streamPtr)
{
	warnStream = streamPtr;
}

/** Set the stream used to post fatal error messages.
 *
 */
void FWMessenger::SetErrorStream(std::ostream* streamPtr)
{
	errorStream = streamPtr;
}

/** Log a neutral message.
 *
 * @param message A message to log.
 */
void FWMessenger::LogMsg(const char* message)
{
	*logStream << message << "\n";
}

/** Log a neutral message with a numeric value.
 *
 * @param message A message to log.
 * @param value A numeric value to appended after the message.  A space is added between them.
 */
void FWMessenger::LogMsg(const char* message, double value)
{
	*logStream << message << " " << value << "\n";
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void FWMessenger::LogMsg(const char* message, std::vector<double> value)
{
	*logStream << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		*logStream << value[i] << ", ";
	}

	*logStream << value[value.size() - 1] << "\n";
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void FWMessenger::LogMsg(const char* message, std::vector<int> value)
{
	*logStream << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		*logStream << value[i] << ", ";
	}

	*logStream << value[value.size() - 1] << "\n";
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void FWMessenger::Warning(const char* message)
{
	*warnStream << message << "\n";
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void FWMessenger::Warning(const std::string& message)
{
	*warnStream << message << "\n";
}

/** Post the passed message and shutdown (not yet implemented!!!!!).
 *
 * @param message An error message.
 */
void FWMessenger::Stop(const char* message)
{
	*errorStream << message << "\n";
	//Without a termination handler this will just result in abort() being called, which is probably
	//acceptable since we only expect Stop() to be called when we can;t recover.
	//Perhaps there is a more robust error throwing approach.
	std::terminate();
}

/** Post the passed message and shutdown (not yet implemented!!!!!).
 *
 * @param message An error message.
 */
void FWMessenger::Stop(const std::string& message)
{
	*errorStream << message << "\n";
	//Without a termination handler this will just result in abort() being called, which is probably
	//acceptable since we only expect Stop() to be called when we can;t recover.
	//Perhaps there is a more robust error throwing approach.
	std::terminate();
}
