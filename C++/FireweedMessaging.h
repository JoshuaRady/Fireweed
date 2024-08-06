/***************************************************************************************************
FireweedMessaging.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  This header file declares a set of
functions for error messaging and logging.

	The Fireweed code may be deployed in multiple ways so the available infrastructure for logging
and error messaging may vary.  These functions provided an interface for basic log messages.  This
is an initial simple implementation that will likely revised or replaced soon.  Right now the plan
is to send messages to standard out by default with means to alter that behavior to be added later.

***************************************************************************************************/
#ifndef FIREWEEDMESSAGING_H
#define FIREWEEDMESSAGING_H

#include <vector>

void LogMsg(const char* message);
void LogMsg(const char* message, double value);
void LogMsg(const char* message, std::vector<double> value);
void LogMsg(const char* message, std::vector<int> value);
void Warning(const char* message);
void Stop(const char* message);

#endif FIREWEEDMESSAGING_H
