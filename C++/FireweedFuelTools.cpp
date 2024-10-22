/***************************************************************************************************
FireweedFuelTools.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 15

	This file is part of the Fireweed wildfire code library.
	This header file defines functions for manipulating the fuel amounts and converting them to
fuel model loadings.

***************************************************************************************************/

#include <algorithm>//For min_element(), is_sorted(), reverse(), sort().
//#include <array>
#include <iterator>//std::distance().
#include <cmath>//std::abs().
#include <numeric>//accumulate()

#include "FireweedFuelTools.h"
#include "FireweedMessaging.h"
#include "FireweedUtils.h"
//#include "FireweedMessaging.h"

/* Take a set of fuel loadings and redistribute the fuel to a second set of size classes.  All the
 * fuel for each input class is placed into the nearest output size class.
 * 
 * Parameters:
 * @param inputSizes A vector of fuel sizes (SAVs, diameters, etc.) for the input fuel loadings.
 * @param loadings A vector of fuel loading masses for each size class in inputSizes.
 * @param outputSizes The size classes to redistribute to.  Units much match inputSizes.
 * 
 * @returns A vector of fuel loadings for each size class specified with outputSizes.  Units will
 * match loadings.
 */
std::vector <double> RedistributeFuelNearest(std::vector <double> inputSizes,
                                             std::vector <double> loadings,
                                             std::vector <double> outputSizes)
{
	std::vector <double> w_o_final(outputSizes.size(), 0);//Return value.

	//Parameter checking:
	if (!SameLengths(inputSizes, loadings))
	{
		Stop("distribSAVs and distribWts must have equal lengths.");
	}

	//If the size classes are the same we are done:
	//if (length(inputSizes) == length(outputSizes) && all(inputSizes == outputSizes))//!!!!!
	if ((inputSizes.size() == outputSizes.size()) && (inputSizes == outputSizes))
	{
		w_o_final = loadings;//w_o_initial
	}
	else//Otherwise redistribute the mass to the requested size classes:
	{
		//Find the closest output size class for each input class:
		for (int i = 0; i < loadings.size(); i++)
		{
			//double diffs[outputSizes.size()];
			//std::array<double, outputSizes.size()> diffs;//!!!!!
			std::vector<double> diffs(outputSizes.size(), 0);
			
			for (int j = 0; j < diffs.size(); j++)
			{
				diffs[j] = std::abs(inputSizes[i] - outputSizes[j]);
			}

			int closest = std::distance(std::begin(diffs),
			                            std::min_element(std::begin(diffs), std::end(diffs)));
			w_o_final[closest] += loadings[i];
		}
	}

	return w_o_final;
}

/** Take a set of fuel loadings and redistribute the fuel to a second set of size classes.  When a
 * fuel size class from the original falls between two of the output set it's fuel is split
 * proportionally based on proximity to its neighbors.
 * 
 * This make sense when there are about one or fewer input classes fall between sets of output
 * classes.  If multiple classes fall between two output classes some fuel from the smaller input
 * class will be contributed to an output size class larger than it adjacent input class.  This isn't
 * logical.  In that case the nearest bin approach will be better.
 * 
 * Parameters:
 * @param inputSizes A vector of fuel sizes (SAVs, diameters, etc.) for the input fuel loadings.
 * @param loadings A vector of fuel loading masses for each size class in inputSizes.
 * @param outputSizes The size classes to redistribute to.  Units much match inputSizes.
 * 
 * @returns A vector of fuel loadings for each size class specified with outputSizes.  Units will
 * match inputSizes.
 */
std::vector <double> RedistributeFuelProportional(std::vector <double> inputSizes,
                                                  std::vector <double> loadings,
                                                  std::vector <double> outputSizes)
{
	std::vector <double> w_o_output(outputSizes.size(), 0);//Return value.

	//Parameter checking:
	if (!SameLengths(inputSizes, loadings))
	{
		Stop("distribSAVs and distribWts must have equal lengths.");
	}

	//The output size classes need to be in order for the algorithm below to work properly.  As a
	//matter of convention we expect SAVs to be in descending order and diameters to be in ascending
	//order. The following code was written primarily for SAV in descending order.  In the case of
	//ascending sizes we reverse the size bins and then reverse the output before returning it below.
	bool sizesAscending = false;

	//if (std::is_sorted(outputSizes))//Sizes are in ascending order.
	if (std::is_sorted(outputSizes.begin(), outputSizes.end()))//Sizes are in ascending order.
	{
		//outputSizes = std::reverse(outputSizes.begin(), outputSizes.end());//In place.
		std::reverse(outputSizes.begin(), outputSizes.end());
		sizesAscending = true;
	}
	//else if (!std::is_sorted(std::reverse(outputSizes.begin(), outputSizes.end()))))//!!!!!
// 	else if (!std::is_sorted(std::reverse(outputSizes.begin(), outputSizes.end())))
// 	{
// 		Stop("outputSizes must be in ascending or descending order.");
// 	}
	else
	{
		//std::vector <double> outputSizesRev = std::reverse(outputSizes.begin(), outputSizes.end());
		std::vector <double> outputSizesRev(outputSizes);
		std::reverse(outputSizesRev.begin(), outputSizesRev.end());
		if (!std::is_sorted(outputSizesRev.begin(), outputSizesRev.end()))
		{
			Stop("outputSizes must be in ascending or descending order.");
		}
	}
	
	//It is not strictly necessary for the input sizes to be ordered to be sorted.  They do need to be
	//in the same order as the output for comparison purposes though:
	//std::vector <double> inputSizesSorted(inputSizes.size(), 0);//!!!!!
	std::vector <double> inputSizesSorted(inputSizes);
	if (sizesAscending)
	{
		//inputSizesSorted = std::sort(inputSizes.begin(), inputSizes.end());//!!!!!
		std::sort(inputSizesSorted.begin(), inputSizesSorted.end());
	}
	else
	{
		//inputSizesSorted =  std::reverse(std::sort(inputSizes.begin(), inputSizes.end()));//!!!!!
		//std::vector <double> sizesSortedAscending = std::sort(inputSizes.begin(), inputSizes.end());
		//inputSizesSorted = std::reverse(sizesSortedAscending.begin(), sizesSortedAscending.end());
		//inputSizesSorted = std::sort(inputSizes.begin(), inputSizes.end());
		//inputSizesSorted = std::reverse(inputSizesSorted.begin(), inputSizesSorted.end());
		std::sort(inputSizesSorted.begin(), inputSizesSorted.end());
		std::reverse(inputSizesSorted.begin(), inputSizesSorted.end());
	}
	
	//If the size classes are the same we are done:
	if ((inputSizes.size() == outputSizes.size()) && (inputSizesSorted == outputSizes))
	{
		w_o_output = loadings;
	}
	else//Otherwise redistribute the mass to the requested size classes:
	{
		//Redistribute the mass for each initial size class to the proper recipient classes://!!!!!
		for (int i = 0; i < loadings.size(); i++)
		{
			for (int j = 0; j < outputSizes.size(); j++)
			{
				if (inputSizes[i] > outputSizes[j])//Greater = to the left of.
				{
					if (j == 1)//The SAV is greater than any in the output.  Put it into the the largest class:
					{
						w_o_output[j] += loadings[i];
						break;
					}
					else
					{
						//Calculate the distance to the adjacent sizes and use to split between them:
						int range = outputSizes[j - 1] - outputSizes[j];
						int diffLeft = outputSizes[j - 1] - inputSizes[i];
						int diffRight = inputSizes[i] - outputSizes[j];
						
						w_o_output[j - 1] += loadings[i] * (diffRight / range);
						w_o_output[j] += loadings[i] * (diffLeft / range);
						break;
					}
				}
				else if (inputSizes[i] == outputSizes[j])//If it is a perfect match put it all in this class.
				{
					w_o_output[j] += loadings[i];
					break;
				}
				else//if (inputSizes[i] < outputSizes[j])
				{
					//If the fuel is smaller than any in the output put it into the the smallest size class:
					if (j == outputSizes.size())
					{
						w_o_output[j] += loadings[i];
						break;
					}
					//Otherwise continue to the next size bin.
				}
			}
		}
	}
	
	if (sizesAscending)
	{
		//w_o_output = std::reverse(w_o_output);
		//w_o_output = std::reverse(w_o_output.begin(), w_o_output.end());
		std::reverse(w_o_output.begin(), w_o_output.end());
	}
	
	return w_o_output;//We could name the values with sizes.
}

/** Distribute a pool of fuel across a set of size classes.
 * 
 *   This function takes an estimated distribution of fuel sizes for a fuel bed.  The distribution is
 * provided as a discreet list of fuel sizes and weighing factors.  The fuel mass input is
 * distributed according to this distribution and is then redistributed to the set of requested
 * output size classes.
 * 
 *   The purpose of this function is to make it easy to convert a total dead fuel estimate to a fuel
 * model loading using any fuel distribution estimate from the literature.  Fuel bed descriptions can
 * take different forms and may not match the size bins for the desired fuel model.  The function
 * handles this conversion in the process.  The function expects weighting factors (fractions) but
 * for convenience fuel loadings may be input and they will be converted to weights.
 * 
 * Parameters:
 * @param distribSizes A vector of fuel sizes (SAVs, diameters, etc.) describing the fuel distribution.
 * @param distribWts A vector of weightings (fractions summing to 1) for each size class in distribSizes.
 *   Alternatively typical fuel loadings (i.e lb/ft^2 | kg/m^2) for each size class in distribSizes.
 * @param totalLoading The total fuel mass to distribute. (lb/ft^2 | kg/m^2)
 * @param outputSizes The SAV size classes to redistribute to.
 * @param method The method for sorting: "Nearest" or "Proportional".
 * 
 * @returns A vector of fuel loadings for each size class in outputSizes.  Units will match
 * totalLoading.
 */
std::vector <double> DistributeFuel(std::vector <double> distribSizes,
                                    std::vector <double> distribWts,
                                    double totalLoading,//std::vector <double> totalLoading,
                                    std::vector <double> outputSizes,
                                    DistribMethod method = Proportional)
{
	std::vector <double> w_o_final(outputSizes.size(), 0);//Return value.

	//Parameter checking:
	if (!SameLengths(distribSizes, distribWts))
	{
		Stop("distribSizes and distribWts must have equal lengths.");
	}
	//Further checking will occur in the called fuctions.
	
	//Make sure the weights are valid and add up to 1 (allowing for floating point error):
	//if (!isTRUE(all.equal(sum(distribWts), 1)))//!!!!!
	//if (!FloatCompare(sum(distribWts), 1)))
	double totalWt = std::accumulate(distribWts.begin(), distribWts.end(), 0);
	if (totalWt != 1)//Temporary until FloatCompare() is moved!!!!!
	{
	//If they aren't we convert them:
		Warning("Adjusting weights to sum to zero.");

		//double totalWt = std::accumulate(distribWts.begin(), distribWts.end(), 0);//!!!!!
		for (int i = 0; i < distribWts.size(); i++)
		{
			distribWts[i] = distribWts[i] / totalWt;
		}
	}
	
	//Distribute the litter mass to the provided distribution. Multiply the mass by the weights:
	std::vector <double> w_o_initial(distribSizes.size(), 0);
	for (int j = 0; j < w_o_initial.size(); j++)
	{
		w_o_initial[j] = totalLoading * distribWts[j];
	}
	
	//Redistribute the mass to the requested size classes:
	if (method == Nearest)
	{
		w_o_final = RedistributeFuelNearest(distribSizes, w_o_initial, outputSizes);
	}
	else if (method == Proportional)
	{
		w_o_final = RedistributeFuelProportional(distribSizes, w_o_initial, outputSizes);
	}
	else//Unnecessary?
	{
		Stop("Unknown method.");
	}
	
	return w_o_final;
}
