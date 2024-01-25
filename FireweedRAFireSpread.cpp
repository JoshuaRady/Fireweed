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


#include <math.h>
#include <vedtor.h>


//Globals:------------------------------------------------------------------------------------------
enum UnitsType {US, Metric};//Move to header?????

//Specify the units to use.  The default is United States customary units.
//This should not be set directly.  Use SetModelUnits().
//private 
UnitsType ModelUnits = US;

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
  return(rho_b);
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
double MeanBulkDensity(vector<double> w_o_ij, double fuelBedDepth)
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
    WARNING << "More than one fuel class expected.";
  }
  
  totalLoading = std::accumulate(w_o_ij.begin(), w_o_ij.end(), 0.0);
  rho_b_bar = totalLoading / fuelBedDepth;
  
  return rho_b_bar;
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
  
  return(xi);
}

//This is a wrapper for PropagatingFluxRatio() that allows it to be called from R:
//The primary reason for providing this wrapper is for verification testing.
extern "C" void PropagatingFluxRatioR(const double* packingRatio, const double* SAV, double* xi)
{
	*xi = PropagatingFluxRatio(*packingRatio, *SAV);
}

//Heat Sink Components:-----------------------------------------------------------------------------
//SECTION TO BE PORTED!!!!!




