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


//Globals:------------------------------------------------------------------------------------------
enum UnitsType {US, Metric};//Move to header?????

//Specify the units to use.  The default is United States customary units.
//This should not be set directly.  Use SetModelUnits().
//private 
UnitsType ModelUnits = US;

//Code:---------------------------------------------------------------------------------------------
//SECTION TO BE PORTED!!!!!



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
extern "C"
{
	void PropagatingFluxRatioR(double* packingRatio, double* SAV, double* xi)
	{
		*xi = PropagatingFluxRatio(*packingRatio, *SAV);
	}
}

//Heat Sink Components:-----------------------------------------------------------------------------
//SECTION TO BE PORTED!!!!!




