/***************************************************************************************************
FireweedCrownFireScottReinhardt.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2026
Reference: Proj. 11 Exp. 26

	This file is part of the Fireweed wildfire code library.  It contains an implementation of the
crown fire equations of Scott & Reinhardt 2001.

***************************************************************************************************/

#ifndef FIREWEEDCROWNFIRESCOTTREINHARDT_H
#define FIREWEEDCROWNFIRESCOTTREINHARDT_H

#include "FireweedFuelModels.h"
#include "FireweedRAFireSpread.h"

FuelModel ConvertToFuelModel10(const FuelModel& fuelModel, FuelModel fuelModel10);
void CheckFuelModel(FuelModel& fuelModel);
void CheckFuelModel(FuelModel& fuelModel, FuelModel fuelModel10);
std::pair<double, double> CheckConvertWindSpeed(const double windSpeed, const double WRF,
                                                const char windType = 'O');

double SpreadRateCrownRothermel(FuelModel fuelModel, const double O, const double slopeSteepness,
                                FuelModel fuelModel10);
double SpreadRateCrownSR(FuelModel fuelModel, const double windSpeed, const double WRF,
                         const double slopeSteepness, const double CBD, const double CBH,
                         const double FMC, FuelModel fuelModel10, const char windType = 'O');

double CriticalCrowningIntensityVanWagner(const double CBH, const double FMC);
double CriticalCrowningROSVanWagner(const double IPrime_initiation, const double HPA);
double CriticalActiveROSVanWagner(const double CBD);

double TorchingIndex(const SpreadCalcs spreadCalcs, const double WRF, const double CBH,
                     const double FMC);
double CrowningIndex(const SpreadCalcs spreadCalcs, const double CBD);
double CrownFractionBurned(FuelModel fuelModel, const double windSpeed, const double WRF, 
                           const double slopeSteepness, const double CBD, const double CBH,
                           const double FMC, FuelModel fuelModel10, const char windType = 'O');

double CrownFireIntensity(FuelModel fuelModel, const double windSpeed, const double WRF,
                          const double slopeSteepness, const double CBD, const double CBH,
                          const double FMC, const double W_canopy, FuelModel fuelModel10,
                          const char windType = 'O', const double H_canopy = 18000);
double CrownHeatPerArea(FuelModel fuelModel, const double windSpeed, const double WRF,
                        const double slopeSteepness, const double CBD, const double CBH,
                        const double FMC, const double W_canopy, FuelModel fuelModel10,
                        const char windType = 'O', const double H_canopy = 18000);

#endif //FIREWEEDCROWNFIRESCOTTREINHARDT_H
