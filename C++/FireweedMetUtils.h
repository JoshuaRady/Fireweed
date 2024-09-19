/***************************************************************************************************
FireweedMetUtils.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:---------------------------------------------------------------------------------------
	This is part of the Fireweed wildfire code library.
	This file declares utility functions for calculatina and converting meteorological properties.

***************************************************************************************************/
#ifndef FIREWEEDMETUTILS_H
#define FIREWEEDMETUTILS_H

double SaturationVaporPressureTetens(double tempC);
double SaturationVaporPressureBuck(double tempC, double p_hPa = 1013);
double RHfromVP(double P, double P_s);
double RHfromDewPointBuck(double tempC, double T_d, double p_hPa = 1013);
double VPDfromRH(double tempC, double rhPct, double P_s);
double VPDfromRHBuck(double tempC, double rhPct, double p_hPa = 1013);

#endif //FIREWEEDMETUTILS_H
