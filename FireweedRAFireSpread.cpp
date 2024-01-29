//FireweedRAFireSpread.cpp
//Programmed by: Joshua M. Rady
//Woodwell Climate Research Center
//Started: 1/12/2024
//Reference: Proj. 11 Exp. 14
//
//Description:--------------------------------------------------------------------------------------
//	This file is part of the Fireweed fire code library.  It contains an R implementation of the
//Rothermel fire spread model (Rothermel 1972) with the modifications of Albini (Albini 1976).
//

//----- This is an ongoing port of the R code.  It is a work in progress.  Please stand by. ------//

//[More header here...]


#include "FireweedRAFireSpread.h"

//Globals:------------------------------------------------------------------------------------------
//enum UnitsType {US, Metric};//Move to header?????

//Specify the units to use.  The default is United States customary units.
//This should not be set directly.  Use SetModelUnits().
//private 
//UnitsType ModelUnits = US;

//Code:---------------------------------------------------------------------------------------------
//SECTION TO BE PORTED!!!!!


//Bulk Density:--------------------------------------------------------------------------------------
//  The bulk density is the mass/wt. of oven dry surface fuel per volume of fuel bed (fuel mass per
//area divided by the fuel bed depth).
//
//Rothermel 1972 equation 40:
//ρb = w_o/δ
//
//Input variables / parameters:
//w_o = Oven dry fuel load (lb/ft^2 | kg/m^2).  This includes combustible and mineral fractions.
//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
//
//Output units: lb/ft^3 | kg/m^3
//The inputs carry the units.  No metric conversions are needed.
double BulkDensity(double w_o, double fuelBedDepth)
{
  double rho_b;//Bulk density.
  
  rho_b = w_o / fuelBedDepth;
  return rho_b;
}

//Mean Bulk Density:
//  The heterogeneous fuel version of the spread equation requires a mean bulk density for the 
//fuels.
//
//Rothermel 1972 equation 74:
//ρb = 1/δ Σi Σj (wo)ij
//For i = 1 to m fuel categories (live vs. dead) and j = 1 to n fuel size classes.
//The original notation includes from and to sum subscripts and bars over rho and w.
//
//Input variables / parameters:
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
//
//Output units: lb/ft^3 | kg/m^3
//The inputs carry the units.  No metric conversions are needed.
double MeanBulkDensity(std::vector<double> w_o_ij, double fuelBedDepth)
{
  
  double totalLoading;//Sum of w_o_ij.
  double rho_b_bar;//Return value.
  
  //Sum the individual fuel loadings across elements:
  //No weights are needed since  w_o is expressed as mass per area.
  //The fuel loading could be a 2D array / matrix but the positions have no significance in the
  //calculation.  Only a sum of all elements needs to be computed.
  //If we know the size of dimensions we could apply error checking.  The 53 standard fire behavior
  //fuel models have 3 dead and 2 live classes.  That is not a fixed requirement of the Rothermel
  //model in theory, but is in practice [I think].
  
  if (w_o_ij.size() < 2)
  {
    Warning("More than one fuel class expected.");
  }
  
  totalLoading = std::accumulate(w_o_ij.begin(), w_o_ij.end(), 0.0);
  rho_b_bar = totalLoading / fuelBedDepth;
  
  return rho_b_bar;
}

//Packing Ratio:-------------------------------------------------------------------------------------
// The packing ratio (beta) is the fraction of the (surface) fuel bed volume occupied by fuel, aka
//compactness.
//
//Rothermel 1972 equation 31:
//β = ρb/ρp
//
//Input variables / parameters:
//rho_b = Fuel array bulk density (lb/ft^3 | kg/m^3).
//rho_p = Fuel particle density (lb/ft^3 | kg/m^3).
//  For the 53 standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
//
//Output units: Dimensionless ratio.
//Input units cancel out.  No metric conversion needed.
double PackingRatio(double rho_b, double rho_p)
{
  return (rho_b / rho_p);
}

//Mean Packing Ratio:
//For heterogeneous fuelbeds the mean packing ratio (beta_bar) must be calculated.
//
//Rothermel 1972 equation 74:
//β = 1/δ Σi Σj (wo)ij/(ρp)ij
//For i = 1 to m fuel categories (live vs. dead) and j = 1 to n fuel size classes.
//The original notation includes from and to sum subscripts and bars over beta, rho, and w.
//
//Input variables / parameters:
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//rho_p_ij = Fuel particle density for each fuel type (lb/ft^3 | kg/m^3).
//  For the 53 standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
//
//Output units: Dimensionless ratio (scalar)
//Input units cancel out.
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth)
{
  int numLoadings, numDensities;
  double meanPackingRatio;//Return value.
  std::vector<double> x;//Intermediate calculation.
  
  //Parameter checking:
  numLoadings = w_o_ij.size();
  numDensities = rho_p_ij.size();
  
  //If only one particle density is provided assume that is it the same for all fuel classes:
  if (numDensities == 1)
  {
    //rho_p_ij.resize(numLoadings, val = rho_p_ij[0])
    rho_p_ij.resize(numLoadings, rho_p_ij[0]);
  }
  else//Otherwise one should be provided for each fuel class.
  {
    if (numDensities != numLoadings)
    {
      Stop("The number of fuel loadings and particle densities do not match.");
    }
  }
  
  //Confirm delta is a scalar: Not needed in C++!
//   if (length(fuelBedDepth) != 1)
//   {
//     //This can be caused if the arguments are out of order.
//     stop("A single fuelbed depth must be provided.")
//     //Add value checking? > 0, < ?
//   }
  
  //Calculate w_o / rho_p for each fuel element:
  //x = w_o_ij / rho_p_ij
  for (int i = 0; i < x.size(); i++)
  {
  	x[i] = w_o_ij[i] / rho_p_ij[i];
  }
  
  meanPackingRatio = std::accumulate(x.begin(), x.end(), 0.0) / fuelBedDepth;//AKA beta_bar
  return meanPackingRatio;
}

//Optimum Packing Ratio:
// This is the packing ratio (compactness) at which the maximum reaction intensity will occur.
//
//Rothermel 1972 equations 37,69:
//βop = 3.348(σ)^-0.8189
//
//Input variables / parameters:
//SAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//  For heterogeneous fuels the SAV of the fuel bed / complex is used.
//
//Output units: Dimensionless ratio
double OptimumPackingRatiofunction(double SAV, UnitsType units)// = ModelUnits
{
  double optPackingRatio;//Return value.
  
  if (units == US)
  {
    optPackingRatio = 3.348 * pow(SAV, -0.8189);
  }
  else// if (units == "Metric")
  {
    //Wilson 1980 gives:
    //optPackingRatio = 0.219685 * SAV^-0.8189
    //Tests show that this is pretty good but I was able to calculate a better conversion as follows:
    
    // SAV is in 1/ft so:
    // SAVcm = SAVft * 1/cmPerFt
    // and
    // SAVft = SAVcm * cmPerFt
    
    //Solve:
    // x = 3.348 * SAV^-0.8189
    // x / 3.348 = SAV^-0.8189
    // (x / 3.348)^(1/-0.8189) = SAV
    // (x / 3.348)^(1/-0.8189) = SAVcm * cmPerFt
    // x / 3.348 = (SAVcm * cmPerFt)^-0.8189
    // x / 3.348 = SAVcm^-0.8189 * cmPerFt^-0.8189
    // x = SAVcm^-0.8189 * cmPerFt^-0.8189 * 3.348
    // x = cmPerFt^-0.8189 * 3.348 * SAVcm^-0.8189
    // x = 0.2039509 * SAVcm^-0.8189
    
    //Same thing:
    // x = 3.348 * SAV^-0.8189
    // x / 3.348 = SAV^-0.8189
    // (x / 3.348)^(1/-0.8189) = SAV
    // (x / 3.348)^(1/-0.8189) = SAVcm * cmPerFt
    // (x / 3.348)^(1/-0.8189) / cmPerFt = SAVcm
    // (x^(1/-0.8189) / 3.348^(1/-0.8189)) / cmPerFt = SAVcm
    // x^(1/-0.8189) / (3.348^(1/-0.8189) * cmPerFt) = SAVcm
    // x^(1/-0.8189) = (3.348^(1/-0.8189) * cmPerFt) * SAVcm
    // x = ((3.348^(1/-0.8189) * cmPerFt) * SAVcm)^-0.8189
    // x = (3.348^(1/-0.8189) * cmPerFt)^-0.8189 * SAVcm^-0.8189
    // x = 0.2039509 * SAVcm^-0.8189
    
    optPackingRatio = 0.2039509 * pow(SAV, -0.8189);
  }
  
  return optPackingRatio;
}

//Heat Source Components:---------------------------------------------------------------------------
//MORE CODE HERE!!!!!



//Propagating Flux Ratio:
//The propagating flux ratio, represented as lower case xi, is the proportion of the reaction
//intensity that heats fuels adjacent to the fire front.
//
//The equation is the same for Rothermel 1972 (eq. 42/76) and Albini 1976 (pg. 4):
//ξ = (192 + 0.2595σ)^-1 exp[(0.792 + 0.681σ^0.5)(β + 0.1)]
//
//Input variables / parameters:
//packingRatio = (Mean) Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//  (dimensionless). For heterogeneous fuels the mean packing ratio is passed in.
//SAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//  For heterogeneous fuels the fuel bed level SAV is used.
//
//Output units: Dimensionless proportion
//R: PropagatingFluxRatio <- function(packingRatio, SAV, units = ModelUnits)
double PropagatingFluxRatio(double packingRatio, double SAV, UnitsType units = ModelUnits)
{
  double xi = 0;//Output
  
  if (units == US)
  {
    xi = pow((192 + 0.2595 * SAV), -1) * exp((0.792 + 0.681 * pow(SAV, 0.5)) * (packingRatio + 0.1));
  }
  else
  {
    xi = pow((192 + 7.90956 * SAV), -1) * exp((0.792 + 3.759712 * pow(SAV, 0.5)) * (packingRatio + 0.1));
    //Wilson 1980 uses:
    //xi = (192 + 7.9095 * SAV)^-1 * exp((0.792 + 3.7597 * SAV^0.5) * (packingRatio + 0.1))
  }
  
  return xi;
}

//This is a wrapper for PropagatingFluxRatio() that allows it to be called from R:
//The primary reason for providing this wrapper is for verification testing.
extern "C" void PropagatingFluxRatioR(const double* packingRatio, const double* SAV, double* xi)
{
	*xi = PropagatingFluxRatio(*packingRatio, *SAV);
}

//Heat Sink Components:-----------------------------------------------------------------------------
//SECTION TO BE PORTED!!!!!


/*Logging:------------------------------------------------------------------------------------------
This code may be deployed in multiple ways so the available infrastructure for logging and error
messaging may vary.  These functions provided an interface for basic log messages.  This is an
initial simple implementation that will likely revised or replaced soon.  Right now the plan is to
send messages to standard out by default with means to alter that behavior to be added later.
--------------------------------------------------------------------------------------------------*/

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
