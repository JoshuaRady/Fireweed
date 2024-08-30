/***************************************************************************************************
FireweedFosbergFuelMoisture.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:----------------------------------------------------------------------------------------
  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
Fosberg dead fuel moisture model, specifically the Rothermel 1981 / NWCG variant.

***************************************************************************************************/
#ifndef FIREWEEDFOSBERGFUELMOISTURE_H
#define FIREWEEDFOSBERGFUELMOISTURE_H

#include <string>
#include <vector>

#include <FireweedUnits.h>

double CtoF(double degreesC);

double FosbergNWCG_1HrFM(std::string tableA_Path, std::string tableB_Path, std::string tableC_Path,
                         std::string tableD_Path, double temp, double rh, int monthOfYear,
                         int hourOfDay, double slopePct, char aspectCardinal, bool shaded = false,
                         char elevation = 'L', UnitsType units = US);
double FosbergNWCG_GetRFM(std::string tableA_Path, double tempF, double rh);
double FosbergNWCG_GetCorrection(std::string tableFilePath, int hourOfDay, double slopePct,
                                 char aspectCardinal, bool shaded, char elevation = 'L');
double NWCG_10hrFM(double oneHrFM);
double NWCG_100hrFM( double oneHrFM);

bool FosbergNWCG_1HrFM_UnitTest(std::string tableA_Path, std::string tableB_Path,
                                std::string tableC_Path, std::string tableD_Path);

#endif //FIREWEEDFOSBERGFUELMOISTURE_H
