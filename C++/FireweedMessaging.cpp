/***************************************************************************************************
FireweedMessaging.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  It included a set of functions for
simple error messaging and logging.

***************************************************************************************************/

#include "FireweedMessaging.h"

/*Logging:------------------------------------------------------------------------------------------
This code may be deployed in multiple ways so the available infrastructure for logging and error
messaging may vary.  These functions provided an interface for basic log messages.  This is an
initial simple implementation that will likely revised or replaced soon.  Right now the plan is to
send messages to standard out by default with means to alter that behavior to be added later.
--------------------------------------------------------------------------------------------------*/

//Log a neutral message:
void LogMsg(const char* message)
{
	std::cout << message << "\n";
}

//Log a neutral message with a numeric value:
void LogMsg(const char* message, double value)
{
	std::cout << message << " " << value << "\n";
}

//Log a neutral message with a numeric vector:
void LogMsg(const char* message, std::vector<double> value)
{
	std::cout << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		std::cout << value[i] << ", ";
	}

	std::cout << value[value.size() - 1] << "\n";
}

//Log a neutral message with a numeric vector:
void LogMsg(const char* message, std::vector<int> value)
{
	std::cout << message << " ";

	for (int i = 0; i < (value.size() - 1); i++)
	{
		std::cout << value[i] << ", ";
	}

	std::cout << value[value.size() - 1] << "\n";
}

//Post a non-fatal warning:
void Warning(const char* message)
{
	std::cout << message << "\n";
}

//Log the passed message and shutdown (not yet implemented):
void Stop(const char* message)
{
	std::cout << message << "\n";
	//Add error throwing.
}
