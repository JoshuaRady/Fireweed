/***************************************************************************************************
FireweedMessaging.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  It defined a class for error
messaging and logging.

	The Fireweed library is modular and may be deployed as part of standalone applications or may be
integrated with existing programs that have there own messaging routines.  This class provides
functions for sending messages to the user or a log.  The default behavior is to send all messages
to standard output.  This behavior can be altered by setting the output streams for each message
type.  In this way messages can be routed into a host program's messaging API.

	Access to the functions is through a global instantiation of the class called Msg or through the
provided convenience functions.

	The R version of the library uses base R messaging functions.  We parallel the function naming
used in R.

***************************************************************************************************/

#include <exception>
#include "FireweedMessaging.h"

//Globals:------------------------------------------------------------------------------------------

/** Use this global instantiation of the FWMessenger class for all Fireweed messaging.
 *
 */
FWMessenger Msg;

//Public Functions:---------------------------------------------------------------------------------

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
void FWMessenger::Log(const char* message)
{
	*logStream << message << std::endl;
}

/** Log a neutral message with a numeric value.
 *
 * @param message A message to log.
 * @param value A numeric value to appended after the message.  A space is added between them.
 */
void FWMessenger::Log(const char* message, double value)
{
	*logStream << message << " " << value << std::endl;
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void FWMessenger::Log(const char* message, std::vector<double> value)
{
	*logStream << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		*logStream << value[i] << ", ";
	}

	*logStream << value[value.size() - 1] << std::endl;
}

/** Log a neutral message with a numeric vector.
 *
 * @param message A message to log.
 * @param value A numeric vector to appended after the message, separated by commas.
 */
void FWMessenger::Log(const char* message, std::vector<int> value)
{
	*logStream << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		*logStream << value[i] << ", ";
	}

	*logStream << value[value.size() - 1] << std::endl;
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void FWMessenger::Warning(const char* message)
{
	*warnStream << message << std::endl;
}

/** Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void FWMessenger::Warning(const std::string& message)
{
	*warnStream << message << std::endl;
}

/** Post the passed message and shutdown.
 *
 * @param message An error message.
 */
void FWMessenger::Stop(const char* message)
{
	*errorStream << message << std::endl;
	//Without a termination handler this will just result in abort() being called, which is probably
	//acceptable since we only expect Stop() to be called when we can;t recover.
	//Perhaps there is a more robust error throwing approach.
	std::terminate();
}

/** Post the passed message and shutdown.
 *
 * @param message An error message.
 */
void FWMessenger::Stop(const std::string& message)
{
	*errorStream << message << std::endl;
	//Without a termination handler this will just result in abort() being called, which is probably
	//acceptable since we only expect Stop() to be called when we can;t recover.
	//Perhaps there is a more robust error throwing approach.
	std::terminate();
}

//External Convenience Functions:-------------------------------------------------------------------

/** Convenience wrapper for FWMessenger::Warning(). Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void Warning(const char* message)
{
	Msg.Warning(message);
}

/** Convenience wrapper for FWMessenger::Warning(). Post a non-fatal warning.
 *
 * @param message A warning message.
 */
void Warning(const std::string& message)
{
	Msg.Warning(message);
}

/** Convenience wrapper for FWMessenger::Stop(). Post the passed message and shutdown.
 *
 * @param message An error message.
 */
void Stop(const char* message)
{
	Msg.Stop(message);
}

/** Convenience wrapper for FWMessenger::Stop(). Post the passed message and shutdown.
 *
 * @param message An error message.
 */
void Stop(const std::string& message)
{
	Msg.Stop(message);
}
