/***************************************************************************************************
FireweedFuelTools.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 15

	This file is part of the Fireweed wildfire code library.
	This header file declares functions for manipulating the fuel amounts and converting them to
fuel model loadings.

***************************************************************************************************/
#ifndef FIREWEEDFUELTOOLS_H
#define FIREWEEDFUELTOOLS_H

#include <vector>

//Constants and Flags:------------------------------------------------------------------------------
enum DistribMethod {Nearest, Proportional};//Is this overkill?

//Functions:----------------------------------------------------------------------------------------

std::vector <double> RedistributeFuelNearest(std::vector <double> inputSizes,
                                             std::vector <double> loadings,
                                             std::vector <double> outputSizes);
std::vector <double> RedistributeFuelProportional(std::vector <double> inputSizes,
                                                  std::vector <double> loadings,
                                                  std::vector <double> outputSizes);
                                                  
std::vector <double> DistributeFuel(std::vector <double> distribSizes,
                                    std::vector <double> distribWts,
                                    double totalLoading,
                                    std::vector <double> outputSizes,
                                    DistribMethod method = Proportional)

#endif //FIREWEEDFUELTOOLS_H
