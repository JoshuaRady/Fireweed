/***************************************************************************************************
FireweedMessaging.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/9/2024
Reference: Proj. 11 Exp. 14 and 11 Exp. 19

	This file is part of the Fireweed wildfire code library.  This header file declares a class for
error messaging and logging.

***************************************************************************************************/
#ifndef FIREWEEDMESSAGING_H
#define FIREWEEDMESSAGING_H

#include <iostream>//Or just <ostream>?
#include <string>
#include <vector>

/** @class FWMessenger
 * @brief A class for messaging (log, warning, error) in Fireweed.
 *
 */
class FWMessenger {
	public:
		FWMessenger();

		void SetLogStream(std::ostream* streamPtr);
		void SetWarnStream(std::ostream* streamPtr);
		void SetErrorStream(std::ostream* streamPtr);

		void Log(const char* message) const;
		virtual void Log(const std::string& message) const;
		void Log(const char* message, double value) const;
		void Log(const char* message, std::vector<double> value) const;
		void Log(const char* message, std::vector<int> value) const;
		void Warning(const char* message) const;
		virtual void Warning(const std::string& message) const;
		void Stop(const char* message) const;
		virtual void Stop(const std::string& message) const;

	protected:
		std::ostream* logStream;
		std::ostream* warnStream;
		std::ostream* errorStream;
};

extern FWMessenger* Msg;//Global interface.

//External functions:

void Warning(const char* message);
void Warning(const std::string& message);
void Stop(const char* message);
void Stop(const std::string& message);

#endif //FIREWEEDMESSAGING_H
