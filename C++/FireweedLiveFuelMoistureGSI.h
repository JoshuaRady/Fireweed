/***************************************************************************************************
FireweedLiveFuelMoistureGSI.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:----------------------------------------------------------------------------------------
  This file is part of the Fireweed wildfire code library.  It declares functions for calculating
GSI live fuel moisture as used in the National Danger Rating System (NFDRS) 2016.

***************************************************************************************************/
#ifndef FIREWEEDLIVEFUELMOISTUREGSI_H
#define FIREWEEDLIVEFUELMOISTUREGSI_H

double GrowingSeasonIndex(double tempCMin, double vpdPa, double dayLength);
double GSI_LiveFuelMoisture(double gsi, double lfmMin, double lfmMax, double gu = 0.5);
double HerbaceousLiveFuelMoisture(double gsi);
double WoodyLiveFuelMoisture(double gsi);

#endif //FIREWEEDLIVEFUELMOISTUREGSI_H
