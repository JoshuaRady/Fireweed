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

//Unit Conversion Factors:---------------------------------------------------------------------------
//Conversion factors marked with an asterisk are not used in this file.  They are provided for use by
//calling code.
//It might be better to have a interface of some sort to request conversion factors from.
//Move to header?????

//Length: (exact per international yard and pound act)
const double cmPerIn = 2.54;
const double cmPerFt = 30.48;
const double mPerFt = 0.3048;
const double ftPerM = 3.28084;//1 / mPerFt
const int ftPerMi = 5280;//* for conversion of windspeed (U, MPH * ftPerMi / 60 = ft/min)

//SAV is in ft^2/ft^3 = 1/ft or cm^2/cm^3 = 1/cm
//Therefore units convert: ft^2/ft^3 * cmPerFt^2/cmPerFt^2 = 1/ft * 1/cmPerFt = 1/cm
//So: SAVft * 1/cmPerFt = SAVft / cmPerFt = SAVcm

//Area:
const double ft2PerAcre = 43560;//*

//Mass:
const double kgPerLb = 0.453592;
const int lbsPerTon = 2000;//*

//Density:
const double lbPerFtCuToKgPerMCu = 16.0185;//kgPerLb * (ftPerM)^3, 16.01846337396

//JPerBtu = 1055.06 or 1,054.35
//The definition of a BTU can vary resulting in several different conversion factors.  Wilson 1980
//seems to have used a value close to the themochemical value of 1.05435 J/BTU, based on his heat of
//preignition conversion.  We will use that to be consistent with his converted constant values.
//The IT value of 1.05506 would be a reasonable alternative.
const double kJPerBtu = 1.05435;

//tons/ac -> lb/ft^2: (See fuel loading note above.)
const double tonsPerAcToLbPerSqFt = lbsPerTon / ft2PerAcre;//*

//Code:---------------------------------------------------------------------------------------------

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
double OptimumPackingRatio(double SAV, UnitsType units)// = ModelUnits
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

//Weighting Factors:---------------------------------------------------------------------------------
//The spread model uses weights in many of the calculations for heterogeneous fuels.  These need to
//be calculated from the surface areas of the different fuel classes.  This function includes the
//changes to Rothermel 1972 by Albini 1976.
//
//Input variables / parameters:
//SAV_ij = Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//rho_p_ij = Fuel particle density for each fuel type (lb/ft^3 | kg/m^3).
//liveDead = An array indicating if each index in each of the other input variables represents a
//  dead (1) or live (2) fuel category.
//
//Output units: unitless weighting factors
//Input units cancel out for calculations.  Metric conversion only needed for SAV sorting.
//
//Note: It makes sense to calculate these together and they only need to be calculated once for a
//give spread rate scenario.
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           std::vector<double> rho_p_ij, std::vector<int> liveDead,
                           UnitsType units)
{
	int numFuelTypes;
	std::vector<double> A_ij(SAV_ij.size(), 0.0);
	double A_i[2] = {0, 0};
	double A_T;//Sum of A_i.
	FuelWeights wts;//Return value.
	double unitFactor;
	//Mapping of fuels to size subclasses:
	//Subclasses are identified as indexes, with -1 indicating an unmapped value, which occurs when
	//there is an empty / undefined fuel class, indicated by an invalid SAV value.
	std::vector<int> subclass_ij(SAV_ij.size(), -1);
	
  //Validity checking:
  //Are arguments the same length?
  if (!SameLengths(SAV_ij, w_o_ij, rho_p_ij, liveDead))
  {
    Stop("CalcWeightings() expects arguments of the same length.");
  }
  
  numFuelTypes = SAV_ij.size();//Types = sum of size classes in both categories.
  
  //Set the size of weight vectors:
  wts.f_ij.resize(numFuelTypes, 0);
  wts.f_i.resize(2, 0);//This should always be length two.  Keeping as a vector for now.
  wts.g_ij.resize(numFuelTypes, 0);
  
  //Calculate the (mean) total surface area for each fuel component:
  //Rothermel equation 53:
  //Aij = (σ)ij (wo)ij ⁄(ρp)ij
  for (int i = 0; i < numFuelTypes; i++)
  {
  	A_ij[i] = SAV_ij[i] * w_o_ij[i] / rho_p_ij[i];
  }
  
  //Mean total surface area by live / dead fuel categories:
  //Rothermel equation 54:
  //Ai = ΣjAij
  for (int k = 0; k < numFuelTypes; k++)
  {
    A_i[liveDead[k]] += A_ij[k];
  }
  
  //Mean total surface area of the fuel:
  //Rothermel equation 55:
  //AT = ΣiAi
  A_T = A_i[0] + A_i[1];//Single scalar value.
  
  //f_ij fuel class weighting factor:
  //Rothermel equation 56:
  //fij = Aij/Ai
  for (int l = 0; l < numFuelTypes; l++)
  {
    if (A_i[liveDead[l]] != 0)
    {
      wts.f_ij[l] = A_ij[l] / A_i[liveDead[l]];
    }
    //A_i can be 0 if there in no fuel in either the live or dead fuel category.  If A_i = 0 then
    //A_ij for this fuel should also be 0.  In this case we avoid a divide by 0 and give an
    //appropriate weight of 0.  Given the initialization we don't have to do this explicitly.
  }
  
  //fi fuel category (live/dead) weighting factor:
  //Rothermel equation 57:
  //fi = Ai/AT
  wts.f_i[0] = A_i[0] / A_T;
  wts.f_i[1] = A_i[1] / A_T;
  
  //g_ij weighting factor:
  //The final set of weights was added in Albini 1976 to get around a logical problem of using f_ij
  //for fuel loading.
  //Albini 1972 pg. 15:
  //g_ij = Σ_(subclass to which j belongs) f_ij
  //This notation is a bit dense (and potentially confusing).  We accomplish this in three steps:
  //1. Determine the size subclass for each fuel type.
  //2. Compute the total weight (from f_ij) in each subclass by live/dead category.
  //3. Set g_ij equal the total weight for the corresponding size subclass.
  
  //What size subclass is each fuel type in?
  //Note: This maps fuel types to size subclasses even when there is no fuel present (loading = 0).
  //Also missing classes, i.e. classes where no SAV is provided, will not be mapped to a subclass.
  //Both these conditions have to be handled below.
  
  //The size subclass ranges are defined by SAV so the units do matter here:
  if (units == Metric)
  {
    unitFactor = 1 / cmPerFt;
  }
  else
  {
    unitFactor = 1;
  }
  //Alternatively we could use an array of range edges with the appropriate units.
  
  for (int n = 0; n < numFuelTypes; n++)
  {
    if (SAV_ij[n] >= 1200 * unitFactor)
    {
      subclass_ij[n] = 0;
    }
    else if (SAV_ij[n] >= 192 * unitFactor)
    {
      subclass_ij[n] = 1;
    }
    else if (SAV_ij[n] >= 96 * unitFactor)
    {
      subclass_ij[n] = 2;
    }
    else if (SAV_ij[n] >= 48 * unitFactor)
    {
      subclass_ij[n] = 3;
    }
    else if (SAV_ij[n] >= 16 * unitFactor)
    {
      subclass_ij[n] = 4;
    }
    else if (SAV_ij[n] > 0)//SAV_ij[n] < 16
    {
      subclass_ij[n] = 5;
    }
    //A value of 0 indicates an empty / undefined SAV value.  Note that this undefined value is
    //specific to this implementation.  How this is indicated in publications varied.  In the
    //original publication of "the 40" 9999 is used.
  }
  
  for (int i = 0; i < 2; i++)//i reused.
  {
    //Calculate the total weight for each size subclass (bin them) for this live/dead category:
    double subclassTotal[6] = {0};//Explicitly initialized to reset for each iteration of the loop.
    
    //The weight of the sixth and largest subclass (index 5) is always 0 so we skip it.
    for (int o = 0; o < 5; o++)
    {
      //Combine the weights of all classes in this size subclass:
      for (int k = 0; k < wts.f_ij.size(); k++)
      {
      	if (liveDead[k] == i && subclass_ij[k] == o)
      	{
      		subclassTotal[o] += wts.f_ij[k];
      	}
      }
    }

    //Assign the subclass weights to each size class.  Some may share the same weight:
    for (int k = 0; k < wts.f_ij.size(); k++)
    {
    	//If a fuel class is not fully specified, i.e. has an invalid SAV of 0, it will not be mapped
    	//to a size subclass.  In that case leave g_ij[k] = 0.  Also don't assign weights to classes
    	//that have no fuel loading.
    	if (liveDead[k] == i && subclass_ij[k] != -1 && w_o_ij[k] != 0)
    	{
    		wts.g_ij[k] = subclassTotal[subclass_ij[k]];
    	}
    	//A value of NA might be more logical but a 0 weight makes the math simpler.
    }
  }
  
  //Return value error checking:
  //Note: if (sum(X) != 1) these comparisons can fail due to small floating point differences
  //when we reassemble the weights.  FloatCompare() handle this problem.
  
  //The dead fuel components of f_ij should always sum to 1:
  if (!FloatCompare(SumByClass(wts.f_ij, liveDead, Dead), 1))
  {
    Stop("f_ij dead fuels do not sum to 1.");
  }
  
  //The live fuel components of f_ij will sum to 1 if present or 0 if not present:
  if (!(FloatCompare(SumByClass(wts.f_ij, liveDead, Live), 0) ||
        FloatCompare(SumByClass(wts.f_ij, liveDead, Live), 1)))
  {
    Stop("Invalid f_ij weights for live fuels.");
  }
  
  //f_i should always sum to 1:
  if (!FloatCompare((wts.f_i[0] + wts.f_i[1]), 1))
  {
    Stop("f_i does not sum to 1.");
  }
  
  //The dead fuel components of g_ij should always sum to 1:
  if (!FloatCompare(SumByClass(wts.g_ij, liveDead, Dead), 1))
  {
    Stop("g_ij dead fuels do not sum to 1.");
  }
  
  //For static models the live fuel components of f_ij will sum to 1 if present or 0 if not present.
  //However, for dynamic fuel models both live classes may be have values of 0 or 1, so sums of 0, 1,
  //and 2 are possible:
  if (!(FloatCompare(SumByClass(wts.g_ij, liveDead, Live), 0) ||
        FloatCompare(SumByClass(wts.g_ij, liveDead, Live), 1) ||
        FloatCompare(SumByClass(wts.g_ij, liveDead, Live), 2)))
  {
    Stop("Invalid g_ij weights for live fuels.");
  }
  
  return wts;
}

//The fuel particle density is often the same across fuel classes.  Allow a single value to be used.
FuelWeights CalcWeightings(std::vector<double> SAV_ij, std::vector<double> w_o_ij,
                           double rho_p, std::vector<int> liveDead, UnitsType units)
{
	std::vector<double> rho_p_ij(SAV_ij.size(), rho_p);//Expand.
	
	return CalcWeightings(SAV_ij, w_o_ij, rho_p_ij, liveDead, units);
}

//This is a wrapper for CalcWeightings() that allows it to be called from R:
//The primary reason for providing this wrapper is for verification testing.
extern "C" void CalcWeightingsR(const double* SAV_ij, const double* w_o_ij, const double* rho_p_ij,
                                const int* liveDead, const int* numFuelTypes, const int* units,
                                double* f_ij, double* f_i, double* g_ij)
{
	FuelWeights wts;
	UnitsType cUnits;
	
	//Convert input arrays to vectors:
	std::vector<double> SAV_ijVec(SAV_ij, SAV_ij + *numFuelTypes);
	std::vector<double> w_o_ijVec(w_o_ij, w_o_ij + *numFuelTypes);
	std::vector<double> rho_p_ijVec(rho_p_ij, rho_p_ij + *numFuelTypes);
	std::vector<int> liveDeadVec(liveDead, liveDead + *numFuelTypes);
	
	//The values used for live and dead differ from R to C++ (due to indexes).  Convert them:
	for (int k = 0; k < *numFuelTypes; k++)
	{
		liveDeadVec[k] -= 1;
	}
	
	//The R code uses a string for units, which can't be passed via .C().  We have to use something
	//as an intermediate translation.  Using the numerical order of options seems as good as any.
	if (*units == 1)
	{
		cUnits = US;
	}
	else if (*units == 2)
	{
		cUnits = Metric;
	}
	else
	{
		Stop("Invalid value passed for units.");//This may not be a R-safe way to abort.  Return an error?
	}
	
	wts = CalcWeightings(SAV_ijVec, w_o_ijVec, rho_p_ijVec, liveDeadVec, cUnits);
	
	//Copy output weights into the return arguments:
	std::copy(wts.f_ij.begin(), wts.f_ij.end(), f_ij);
	std::copy(wts.f_i.begin(), wts.f_i.end(), f_i);
	std::copy(wts.g_ij.begin(), wts.g_ij.end(), g_ij);
}

//Fuel Bed Surface-area-to-volume Ratio:-------------------------------------------------------------
//  For heterogeneous fuels a SAV for the entire fuel bed must be calculated.  This is frequently
//referred to as the characteristic SAV.  It is a weighted average of the fuel component SAVs.
//
//Input variables / parameters:
//SAV_ij =	Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
//f_ij = Weighting factors for each fuel type (dimensionless).
//f_i = Weighting factors for each fuel live/dead category (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//  dead (1) or live (2) fuel category.
//
//Output units: ft^2/ft^3 | cm^2/cm^3
//The inputs carry the units.  No metric conversions are needed.
//
//Note: We could pass in a FuelWeights object rather than f_ij and f_i but this keeps the interface
//the same as the R version.
double FuelBedSAV(std::vector<double> SAV_ij, std::vector<double> f_ij, std::vector<double> f_i,
                  std::vector<int> liveDead)
{
  int numFuelTypes;
  double SAV_i[2];//SAV by live / dead category.
  double fuelBedSAV;//Return value.
  
  //Argument checking:
  if (!SameLengths(SAV_ij, f_ij, liveDead))
  {
    Stop("FuelBedSAV() expects arguments of the same length.");
  }
  
  numFuelTypes = SAV_ij.size();//Types = sum of size classes in both categories.
  
  //Characteristic live and dead SAVs:
  //Rothermel 1972 equation 72:
  //σi = Σj fijσij (~ over sigma sub i and bar over sigma sub ij in original)
  
  for (int k = 0; k < numFuelTypes; k++)
  {
    //SAV_i[liveDead[k]] = SAV_i[liveDead[k]] + (f_ij[k] * SAV_ij[k]);
    SAV_i[liveDead[k]] += f_ij[k] * SAV_ij[k];
  }
  
  //Sum the live and dead components to get the final value:
  //Rothermel 1972 equation 71:
  //σ = Σi fiσi (~ over sigma and sigma sub i in original)
  fuelBedSAV = (f_i[0] * SAV_i[0]) + (f_i[1] * SAV_i[1]);//Or:
  //fuelBedSAV = (f_i[Dead] * SAV_i[Dead]) + (f_i[Live] * SAV_i[Live]);
  
  return fuelBedSAV;
}

//Net Fuel Load:-------------------------------------------------------------------------------------
// The net fuel load (w_n) is mass per ground area of dry combustible fuel removing the
//noncombustible mineral mass.  This is the revised calculation developed by Albini.
//
//Albini 1976 pg. 14:
//wn = wo(1 - ST)
//
//Input variables / parameters:
//w_o = Oven dry fuel load (lb/ft^2 | kg/m^2).  This includes combustible and mineral fractions.
//S_T = Total mineral content (unitless fraction: mineral mass / total dry mass).
//  For all standard fuel models this is 5.55% (0.0555).
//
//Output units: lb/ft^2 | kg/m^3
//The inputs carry the units.  No metric conversions are needed.
double NetFuelLoad_Homo(double w_o, double S_T)
{
	double w_n;//Return value.
	
	w_n = w_o * (1 - S_T);
	
	return w_n;
}

//For heterogeneous fuels the net fuel load for each fuel category (live/dead) is calculated using
//weights.
//
//Input variables / parameters:
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//S_T_ij = An array of total mineral content for each fuel type (unitless fraction).
//g_ij = Net fuel load weights (Albini 1976).
//
//Returns: w_n_i = Net fuel load for live/dead fuel categories.
//Output units: lb/ft^2 | kg/m^3
//The inputs carry the units.  No metric conversions are needed.
std::vector <double> NetFuelLoad_Het(std::vector <double> w_o_ij, std::vector <double> S_T_ij,
                                     std::vector <double> g_ij, std::vector <int> liveDead)
{
	int numFuelTypes;
	std::vector <double> w_n_ij(w_o_ij.size(), 0);//Intermediate
	std::vector <double> w_n_i(2, 0);//Return value.
	
	//Argument checking:
	if (!SameLengths(w_o_ij, S_T_ij, g_ij, liveDead))
	{
		Stop("NetFuelLoad_Het() expects arguments of the same length.");
	}
	
	numFuelTypes = w_o_ij.size();
	
	//w_n_ij.resize(numFuelTypes, 0);
	
	for (int k = 0; k < numFuelTypes; k++)
	{
		//Calculate the net fuel load for each fuel class:
		//Rothermel equation 60 modified by Albini 1976 pg. 14:
		//(wn)ij = (wo)ij (1 – (ST)ij)
		w_n_ij[k] = w_o_ij[k] * (1 - S_T_ij[k]);
		
		//Sum the net fuel load for each fuel category:
		//Albini 1976 pg. 15:
		//(wn)i = Σjgij(wn)ij
		w_n_i[liveDead[k]] += g_ij[k] * w_n_ij[k];
	}
	
	return w_n_i;
}

//Damping Coefficients:-----------------------------------------------------------------------------

//Moisture Damping Coefficient (homogeneous fuels):
// This returns the extent to which fuel moisture reduces combustion for one fuel component.
//
//Input variables / parameters:
//M_f = Fuel moisture content (fraction: water weight/dry fuel weight).
//M_x = Moisture of extinction (fraction: water weight/dry fuel weight).
//
//Output units: Dimensionless coefficient
//Input units cancel out.  No metric conversion needed.
double MoistureDampingCoefficient_Homo(double M_f, double M_x)
{
	double r_M;//Ratio of fuel moisture content to moisture of extinction.
	double eta_M;//Return value.
	
	//Moisture content can well exceed 1 for live fuels (at least to 300%):
	if (!InRange(M_f, 0, 3.5))
	{
		Stop("Suspect moisture content (M_f).");
	}
	//See MoistureDampingCoefficient_Het() for notes on valid moistures of extinction:
	if (!InRange(M_x, 0, 8))
	{
		Stop("Suspect moisture of extinction.");
	}
	
	//Calculate the ratio of fuel moisture content to moisture of extinction:
	//Rothermel 1972 equation 29,65 with maximum added in Albini 1976:
	//Note: This uses the notation of Andrews 2018, the original is a bit different.
	//rM = Mf/Mx (max = 1.0)
	r_M = M_f / M_x;
	
	if (r_M > 1.0)
	{
		r_M = 1.0;
	}
	
	//Use the ratio to calculate the damping coefficient:
	//Rothermel 1972 equations 29,64:
	//ηM = 1 - 2.59rM + 5.11(rM)^2 - 3.52(rM)^3
	eta_M = 1 - 2.59 * r_M + 5.11 * pow(r_M, 2) - 3.52 * pow(r_M, 3);
	
	return eta_M;
}

//Moisture Damping Coefficient (heterogeneous fuels):
//  For heterogeneous fuel beds the moisture damping coefficient is calculated for each fuel category
//(live/dead).
//
//Input variables / parameters:
//M_f_ij = Fuel moisture content for for each fuel type (fraction: water weight/dry fuel weight).
//M_x_i = Moisture of extinction each fuel category (fraction: water weight/dry fuel weight).
//f_ij = Weighting factors for each fuel type (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//  dead (1) or live (2) fuel category.
//
//Output units: Dimensionless coefficient (array length 2)
//Input units cancel out.  No metric conversion needed.
std::vector <double> MoistureDampingCoefficient_Het(std::vector <double> M_f_ij,
                                                    std::vector <double> M_x_i,
                                                    std::vector <double> f_ij,
                                                    std::vector <double> liveDead)
{
	int numFuelTypes;
	std::vector <double> M_f_i(2, 0);//Weighted moisture content.
	std::vector <double> eta_m_i(2, 0);//Return value.
	
	if (!SameLengths(M_f_ij, f_ij, liveDead))
	{
		Stop("MoistureDampingCoefficient_Het() expects arguments of the same length.");
	}
	if (!InRange(M_f_ij, 0, 3.5))
	{
		Stop("Suspect moisture content (M_f_ij).");
	}
	//Dead fuels have moisture of extinction values with a range of 12-40% for the standard models,
	//though higher might be possible so we add a little wiggle room:
	if (!InRange(M_x_i[Dead], 0, 0.5))
	{
		Stop("Invalid dead fuel moisture of extinction.");
	}
	//Calculated live fuel moisture of extinction can reach over 700%, even though that moisture level
	//is not physiologic:
	if (!InRange(M_x_i[Live], 0, 8))
	{
		Stop("Suspect live fuel moisture of extinction.");
	}
	
	numFuelTypes = M_f_ij.size();
	
	//Calculate the weighted moisture content for each fuel category:
	//Rothermel 1972 equations 66:
	//(Mf)i = Σj fij (Mf)ij
	for (int k = 0; k < numFuelTypes; k++)
	{
		M_f_i[liveDead[k]] += f_ij[k] * M_f_ij[k];
	}
	
	//Calculate the moisture damping coefficient for each fuel category:
	//Rothermel 1972 equations 64,65:
	//(ηM)i = 1 – 2.59(rM)i + 5.11(rM)i2 – 3.52(rM)i3 (max = 1)
	eta_m_i[Dead] = MoistureDampingCoefficient_Homo(M_f_i[Dead], M_x_i[Dead]);
	eta_m_i[Live] = MoistureDampingCoefficient_Homo(M_f_i[Live], M_x_i[Live]);
	
	return eta_m_i;
}

//Live Fuel Moisture of Extinction:
//  The live fuel moisture of extinction determines if live fuels will burn and contribute to the
//heat source term in the calculations.  The live fuel moisture of extinction is calculated from the
//dead fuel moisture of extinction and in relation to the ratio of live to dead fuels. The value can
//be quite variable from near zero to over 700 percent (see Andrews 2018 section 5.3.2.2).
//
//Input variables / parameters:
//M_f_ij = Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
//M_x_1 = Dead fuel moisture of extinction (fraction: water weight/dry fuel weight).
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//SAV_ij =	Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
//liveDead = An array indicating if each index in each of the other input variables represents a
//  dead (1) or live (2) fuel category.
//
//Output units: fraction, water weight/dry fuel weight
double LiveFuelMoistureOfExtinction(std::vector <double> M_f_ij, double M_x_1,
                                    std::vector <double> w_o_ij, std::vector <double> SAV_ij,
                                    std::vector <int> liveDead, UnitsType units)
{
	int numFuelTypes;
	double liveSum = 0;
	double deadSum = 0;
	double W;//Live / dead ratio.
	double common;
	double top = 0;//The numerator sum.
	double bottom = 0;//The denominator sum.
	double M_f_dead;
	double M_x_2;//Return value.
	
	if (!SameLengths(M_f_ij, w_o_ij, SAV_ij, liveDead))
	{
		Stop("LiveFuelMoistureOfExtinction() expects arguments of the same length.");
	}
	if (!InRange(M_f_ij, 0, 3.5))
	{
		Stop("Suspect moisture content (M_f_ij).");
	}
	if (!InRange(M_x_1, 0, 0.5))
	{
		Stop("Invalid dead fuel moisture of extinction.");
	}
	
	numFuelTypes = M_f_ij.size();
	
	//Changing the equations is more complicated than changing the inputs.
	if (units == Metric)
	{
		for (int k = 0; k < numFuelTypes; k++)
		{
			w_o_ij[k] = w_o_ij[k] / lbPerFtCuToKgPerMCu;
			SAV_ij[k] = SAV_ij[k] * cmPerFt;//1/cm to 1/ft
		}
	}
	
	//Calculate dead:live loading ratio, notated W:
	//Albini 1976 pg. 16:
	//W = Σj(wo)1jexp(-138/σ1j) / Σj(wo)2jexp(-500/σ2j)
	for (int k = 0; k < numFuelTypes; k++)
	{
		if (liveDead[k] == Live)
		{
			liveSum += w_o_ij[k] * exp(-138 / SAV_ij[k]);
		}
		else//(liveDead[k] == Dead)
		{
			deadSum += w_o_ij[k] * exp(-500 / SAV_ij[k]);
		}
	}
	
	//If the loading for the live fuel categories are all 0 or live categories are missing liveSum will
	//be zero.  The ratio W will be therefore also be zero.  This in turn will result in M_x_2 = M_x_1.
	//While not conceptual meaningful this has no has no mathematical consequence downstream.  Forcing
	//the value to NA or 0 would cause mathematical problems downstream..
	
	W = liveSum / deadSum;//Unitless ratio.
	
	//Calculate fine dead fuel moisture as:
	//Albini 1976 pg. 16:
	//Mf,dead = Σj(Mf)1j(wo)1jexp(–138/σ1j) / Σj(wo)1jexp(–138/σ1j)
	for (int k = 0; k < numFuelTypes; k++)
	{
		if (liveDead[k] == Dead)
		{
			common = w_o_ij[k] * exp(-138 / SAV_ij[k]);
			top += M_f_ij[k] * common;
			bottom += common;
		}
	}
	
	M_f_dead = top / bottom;//Moisture fraction / unitless.
	
	//Calculate the live fuel moisture of extinction ((Mx)2):
	//Rothermel 1972 equation 88 with Albini 1976 pg. 16 modifications:
	//(Mx)2 = 2.9W[1 – Mf,dead⁄(Mx)1] - 0.226, (min = (Mx)1)
	M_x_2 = 2.9 * W * (1 - M_f_dead / M_x_1) - 0.226;
	
	if (M_x_2 < M_x_1)
	{
		M_x_2 = M_x_1;
	}
	
	return M_x_2;
}

//Mineral Damping Coefficient (homogeneous fuels):
//
//Albini 1976 pg. 14 adds a maximum to Rothermel 1972 equation 30/62:
//ηs = 0.174Se^-0.19 (max = 1.0)
//
//Input variables / parameters:
//S_e = Effective mineral content (unitless fraction: (mineral mass – mass silica) / total dry mass).
//  For all standard fuel models this is 1% (0.01).
//
//Output units: Dimensionless coefficient
//Unitless inputs and outputs.  No metric conversion needed.
double MineralDampingCoefficient_Homo(double S_e)
{
	double eta_s;//Return value.
	
	//Parameter checking:
	if (!ValidProportion(S_e))
	{
		Stop("Effective mineral content must be from 0-1.");
	}
	
	eta_s = 0.174 * pow(S_e, -0.19);
	
	if (eta_s < 1.0)
	{
		return eta_s;
	}
	else
	{
		return 1.0;
	}
}

//Mineral Damping Coefficient (heterogeneous fuels):
//  For heterogeneous fuels the mineral damping coefficient is calculated for each fuel category
//(live/dead).
//
//Input variables / parameters:
//S_e_ij = Effective mineral content for each fuel type (unitless fraction:
//  (mineral mass – mass silica) / total dry mass).
//f_ij = Weighting factors for each fuel type (dimensionless).
//
//Output units: Dimensionless coefficient (array of 2)
//Unitless inputs and outputs.  No metric conversion needed.
std::vector <double> MineralDampingCoefficient_Het(std::vector <double> S_e_ij,
                                                   std::vector <double> f_ij,
                                                   std::vector <int> liveDead)
{
	double numFuelTypes;
	double S_e_i[2] = {0, 0};//Effective mineral content by live / dead category.
	std::vector <double> eta_s_i(2, 0);//Return value.
	
	//Parameter checking:
	if (!SameLengths(S_e_ij, f_ij, liveDead))
	{
		Stop("MineralDampingCoefficient_Het() expects arguments of the same length.");
	}
	if (!ValidProportion(S_e_ij))
	{
		Stop("Effective mineral content must be from 0-1.");
	}
	
	numFuelTypes = S_e_ij.size();//Types = sum of size classes in both categories.
	
	//Calculate the weighted effective mineral content for each fuel category:
	//(Se)i = Σj fij (Se)ij
	for (int k; k < numFuelTypes; k++)
	{
		S_e_i[liveDead[k]] += f_ij[k] * S_e_ij[k];
	}
	
	//Caculate the mineral damping coefficient for each fuel category:
	//(ηs)i = 0.174(Se)i^–0.19 (max = 1)
	eta_s_i[Dead] = MineralDampingCoefficient_Homo(S_e_i[Dead]);
	eta_s_i[Live] = MineralDampingCoefficient_Homo(S_e_i[Live]);
	
	return eta_s_i;
}

//Slope and Wind Factors:---------------------------------------------------------------------------

//Slope factor:
//Dimensionless multiplier that accounts for the effect of slope on spread behavior.  Same for
//homogeneous and heterogeneous fuels.
//
//Rothermel 1972 equation 51,80:
//ϕs = 5.275β^-0.3(tan ϕ)^2
//
//Input variables / parameters:
//packingRatio = Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//  (dimensionless).
//slopeSteepness = Slope steepness maximum (unitless fraction: vertical rise / horizontal distance).
//  AKA tan ϕ.
//
//Output units: Dimensionless adjustment factor
//Inputs are fractions which do not change with units.  Mo metric conversion required.
double SlopeFactor(double packingRatio, double slopeSteepness)
{
	double phi_s;

	phi_s = 5.275 * pow(packingRatio, -0.3) * pow(slopeSteepness, 2);
	return phi_s;
}

//Wind factor:
//Dimensionless multiplier that accounts for the effect of wind speed on spread behavior
//(propagating flux ratio specifically).  Same for homegenous and heterogeneous fuels.
//
//Input variables / parameters:
//SAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//packingRatio = Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//  (dimensionless).
//optPackingRatio = Optimum packing ratio (dimensionless).
//   Note: optimum packing ratio is a function of SAV.
//U = Wind speed at midflame height (ft/min | m/min).
//
//Output units: Dimensionless
//
//Note: It is not possible to calculate if a wind limit is indicated internal to this function
//and not all authors agree that a wind limit should be used.  U should be capped, if deemed
//appropriate prior to passing it in to this function.
double WindFactor(double SAV, double packingRatio, double optPackingRatio, double U,
                  UnitsType units)
{
	double C, B, E;//Partial terms.
	double phi_wl;//Return value.

	C = WindFactorC(SAV, units);
	B = WindFactorB(SAV, units);
	E = WindFactorE(SAV, units);

	if (units == US)
	{
		//Rothermel 1972 equation 47,79:
		//ϕw = CU^B(β/βop)^-E
		phi_w = C * pow(U, B) * pow((packingRatio / optPackingRatio), -E);
	}
	else
	{
		phi_w = C * pow((ftPerM * U), B) * pow((packingRatio / optPackingRatio), -E);
		//Wilson 1980 uses:
		//phi_w = C * (0.3048 * U)^B * (packingRatio / optPackingRatio)^-E
	}

	return phi_w;
}

//The following three functions represent the equations internal to the calculation of the wind
//factor.  They have been broken out because they are also used in EffectiveWindSpeed().
//The metric conversions agree with Andrews 2018 except for the number of digits.
//Should significant digits be observed for the conversions here?
//Having a default for units argument is probably unnecessary since these will probably never be
//called directly.

double WindFactorC(double SAV, UnitsType units)
{
	double C;//Return value.

	if (units == US)
	{
		//C = unnamed term
		//Rothermel 1972 equation 48,82:
		//C = 7.47exp(-0.133σ^0.55)
		C = 7.47 * exp(-0.133 * pow(SAV, 0.55));
	}
	else
	{
		C = 7.47 * exp(-0.8710837 * pow(SAV, 0.55));//-0.133 * cmPerFt^0.55 = -0.8710837
		//Wilson 1980 uses:
		//C = 7.47 * exp(-0.8711 * SAV^0.55)
	}
	return C;
}

double WindFactorB(SAV, UnitsType units)
{
	double B;//Return value.

	if (units == US)
	{
		//B = unnamed term
		//Rothermel 1972 equation 49,83:
		//B = 0.02526σ^0.54
		B = 0.02526 * pow(SAV, 0.54);
	}
	else
	{
		B = 0.1598827 * pow(SAV, 0.54);//0.02526 * cmPerFt^0.54 = 0.1598827
		//Wilson 1980 uses:
		//B = 0.15988 * SAV^0.54
	}
	return B;
}

double WindFactorE(SAV, UnitsType units)
{
	double E;//Return value.

	if (units == US)
	{
		//E = unnamed term
		//Rothermel 1972 equation 50,84:
		//E = 0.715exp(-3.59×10^-4 σ)
		E = 0.715 * exp(-0.000359 * SAV)
	}
	else
	{
		E = 0.715 * exp(-0.01094232 * SAV)//-0.000359 * cmPerFt = -0.01094232
		//Wilson 1980 uses:
		//E = 0.715 * exp(-0.01094 * SAV)
	}
	return E;
}

//In the discussion of the wind factor in Andrews 2018 the equation is simplified to AU^B, combining
//the multiplicative terms into a single factor A.  This function is provided to return this value
//for verification purposes and is not currently needed to compute model outputs.
//Note: This currently does not include the unit conversion of U when U is metric in A.  Since A is
//used for diagnostic purposes that seems the correct approach.
double WindFactorA(double SAV, double packingRatio, double optPackingRatio, UnitsType units)
{
	double C, E, A;//Partial terms.

	C = WindFactorC(SAV, units);
	E = WindFactorE(SAV, units);
	A = C * pow((packingRatio / optPackingRatio), -E);
	return A;
}

//Wind limit:
//  The "wind limit" or "maximum reliable wind" is used to limit the effect of wind on the fire
//spread rate as wind speed gets high.  It caps the wind speed at a value that is a function of the
//reaction intensity.
//   There is not agreement on whether the wind limit should be used.  Albini chose to not use it,
//but his code reports if the limit was reached (Albini 1976, pg 26).  More recent work finds the
//original calculation to be flawed and presents an alternate formulation from (Andrews et. al 2013).
//However, They conclude that in general neither should be used.  They state a better alternative is
//to cap the spread rate at the “effective wind speed”.
//  We implement the original formulation as an option to be able to reproduce results that do use
//the wind limit.
//
//Input variables / parameters:
//U = Wind speed at midflame height (ft/min | m/min).
//I_R = Reaction intensity (Btu/ft^2/min | kJ/m^2/min).
//
//Output units: adjusted wind speed (U) at midflame height (ft/min | m/min)
double WindLimit(double U, double I_R, FormUnits units = ModelUnits)
{
	double threshold;//Windlimit threshold.

	if (units == US)
	{
		threshold = 0.9;
	}
	else
	{
		threshold = 0.02417144;
	}
	
	//Rothermel 1972 Equation 87:
	if (U/I_R > threshold)
	{
		U = threshold * I_R;
		
		//Or Andrews et al. 2013 equation 21:
		//U = 96.8 * I_R^(1/3)
	}
	//Otherwise return U unchanged.
	
	return U;
}

//Heat Source Components:---------------------------------------------------------------------------

//Optimum (Potential) Reaction Velocity:
//  This is a measure of the optimum (potential) fuel consumption rate (fire efficiency / reaction
//time).  The 'optimum' rate is the ideal rate that would occur for alpha cellulose in the absence
//of minerals and moisture.
//Notation: Γ' (Gamma prime)
//
//  This version includes the Albini 1976 modification.
//
//Input variables / parameters:
//β = (Mean) Packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
//packingRatio = (Mean) Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//  (dimensionless). For heterogeneous fuels the mean packing ratio is passed in.
//SAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//  For heterogeneous fuels the SAV of the fuel bed / complex is used.
//
//Output units: min^-1
double OptimumReactionVelocity(double packingRatio, double SAV, UnitsType units)
{
	double optPackingRatio;
	double A;//Intermediate.
	double GammaPrime;//Return value.

	//Calculate the maximum reaction velocity (min^-1):
	//This is the rate for moisture free fuel with mineral composition of alpha cellulose.
	//Rothermel 1972 equations 36,68:
	//Γ'max = σ^1.5/(495 + 0.0594σ^1.5)
	if (units == US)
	{
		GammaPrimeMax = pow(SAV, 1.5) / (495 + 0.0594 * pow(SAV, 1.5));
		//Or equivalently:
		//GammaPrimeMax = 1 / (0.0594 * 495 / SAV^1.5)
	}
	else
	{
		GammaPrimeMax = 1 / (0.0594 + 2.941594 / pow(SAV, 1.5));
		//Wilson 1980 uses:
		//GammaPrimeMax = (0.0591 + 2.926 * SAVcm^-1.5)^-1 = 1 / (0.0591 + 2.926 / SAVcm^1.5)
	}
	
	optPackingRatio = OptimumPackingRatio(SAV, units);
	
	//"Arbitrary" variable:
	//Albini 1976 pg. 15:
	//A = 133σ^-0.7913
	if (units == "US")
	{
		A = 133 * pow(SAV, -0.7913);
	}
	else
	{
		A = 8.903291 * pow(SAV, -0.7913);
		//Wilson 1980 uses:
		//A = 8.9033 * SAV^-0.7913
	}
	
	//These are combined to produce the optimal reaction velocity (min^-1):
	//Rothermel 1972 equation 38:
	//Γ' = Γ'max(β/βop)^A exp[A(1 - β/βop)]
	GammaPrime = GammaPrimeMax * pow((packingRatio/optPackingRatio), A)  *
		exp(A * (1 - packingRatio/optPackingRatio));
	
	return GammaPrime;
}

//Live / Dead Heat Content:
//  Calculate the weighted (low?) heat content of the live / dead fuel categories.
//Only used by reaction intensity calculation.
//
//Input variables / parameters:
//h_ij = Heat content of the fuel types (Btu/lb | kJ/kg).
//f_ij = Weighting factors for each fuel type (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//  dead (1) or live (2) fuel category.
//
//Output units: btu/lb | kJ/kg
//Whatever units are input, the same will come out.  No unit conversions needed.
double LiveDeadHeatContent(std::vector <double> h_ij, std::vector <double> f_ij,
                           std::vector <double> liveDead)
{
	std::vector <double> h_i(2, 0);//Could this cause issues if only live or dead fuel is present?
	int numFuelTypes;

	numFuelTypes = f_ij.size();

	if (!SameLengths(h_ij, f_ij, liveDead))
	{
		Stop("LiveDeadHeatContent() expects arguments of the same length.");
	}

	//Rothermel 1972 equation 61:
	//hi = Σj fij hij
	for (int k = 0; k < numFuelTypes; k++)
	{
		h_i[liveDead[k]] = h_i[liveDead[k]] + f_ij[k] * h_ij[k];
	}

	return h_i;
}

//Reaction Intensity:
//  The reaction intensity (I_R) is the total energy released by the fire front in Btu/ft^2/min in
//all forms (radiation, conduction, and convection).
//  This it not the same as fireline intensity!

//Reaction Intensity, Rothermel version for homogeneous fuels:
//
//Rothermel equation 27:
//IR = Γ'wnhηMηs
//
//Input variables / parameters:
//GammaPrime = Optimum reaction velocity (min^-1).
//w_n = Net fuel load for a single fuel component (lb/ft^2 | kg/m^2).
//h = Heat content of the fuel type (Btu/lb | kJ/kg).
//eta_M = Moisture damping coefficient (unitless).
//eta_s = Mineral damping coefficient (unitless).
//
//Output units: Btu/ft^2/min | kJ/m^2/min
//Inputs carry units.  No unit conversions are needed.
//
//Note: The alternate parameters (GammaPrime, w_n, h, M_f, M_x, Se) could be used.
double ReactionIntensity_Homo (double GammaPrime, double w_n, double h, double eta_M, double eta_s)
{
	double I_R;//Return value.

	I_R = GammaPrime * w_n * h * eta_M * eta_s;

	return I_R;
}

//Reaction Intensity for Heterogeneous Fuels:
//
//Rothermel equation 58 modified by Albini 1976 pg. 17:
//IR = Γ' Σi (wn)ihi(ηM)i(ηs)i
//
//Input variables / parameters:
//GammaPrime = Optimum reaction velocity (min^-1).
//w_n_i = Net fuel load for live/dead fuel categories (lb/ft^2 | kg/m^2).
//h_i = Heat content for live/dead fuel (Btu/lb | kJ/kg).
//eta_M_i = Moisture damping coefficient for live/dead fuel categories (unitless).
//eta_s_i = Mineral damping coefficient for live/dead fuel categories (unitless).
//
//Output units: Btu/ft^2/min | kJ/m^2/min
//Inputs carry units.  No unit conversions are needed.
double ReactionIntensity_Het(double GammaPrime, std::vector <double> w_n_i,
                             std::vector <double> h_i, std::vector <double> eta_M_i,
                             std::vector <double> eta_s_i)
{
	double I_R;//Return value.

	if (!SameLengths(w_n_i, h_i, eta_M_i, eta_s_i))
	{
		Stop("ReactionIntensity_Het() expects arguments of the same length.");
	}

	I_R = GammaPrime * ((w_n_i[0] * h_i[0] * eta_M_i[0] * eta_s_i[0]) +
	                    (w_n_i[1] * h_i[1] * eta_M_i[1] * eta_s_i[1]));

	return I_R;
}

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
double PropagatingFluxRatio(double packingRatio, double SAV, UnitsType units = ModelUnits)//Default to header!!!!!
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


//Utilities:----------------------------------------------------------------------------------------

//Return the heat content (h) used in the 53 standard fuel models in the appropriate units:
double StdHeatContent(UnitsType units)
{
  double h;//Return value.
  
  if (units == US)
  {
    h = 8000;//Btu/lb
  }
  else//if (units == Metric)
  {
    h = 18595.57;//kJ/kg, 8000 * kJPerBtu / kgPerLb
  }
  return h;
}

//Return the fuel particle density (rho_p) used in the 53 standard fuel models in the appropriate
//units:
double StdRho_p(UnitsType units)
{
  double rho_p;//Return value.
  
  if (units == US)
  {
    rho_p = 32;//lb/ft^3
  }
  else//if (units == Metric)
  {
    rho_p = 512.592;//kg/m^3, (32 * lbPerFtCuToKgPerMCu)
  }
  return rho_p;
}

//[InitSpreadParam()...]

//SameLengths():
//This utility checks that the parameters (vectors) passed have the same length.  Between 2 and 4
//arguments are accepted.
/*
Overloading has been used to reproduce the behavior of the R function, although that function can
handle arguments of different types in any order.  We only handle a subset of possible type
combinations.  This function is currently used primarily variable vectors of type double and
liveDead, which is currently an integer but could be a boolean.  To reduce the number of
possibilities we require liveDead to be last.  It might be possible to do this more compactly with
template functions or ariadic arguments.

Could change arguments to const &?
 */
//SameLengths <- function(arg1, arg2, arg3 = NULL, arg4 = NULL)
bool SameLengths(std::vector<double> arg1, std::vector<double> arg2)
{
  //Put the argments in a list removing NULL elements:
  // argList = list(arg1, arg2, arg3, arg4)
//   argList = argList[!sapply(argList, is.null)]
//   //This would work too for omitted arguments but would ignore any zero length vectors passed in:
//   //argList = argList[length(argList) != 0]
//   
//   //Are arguments the same length?
//   if (all(sapply(argList, length) == length(arg1)))
//   {
//     //return(length(arg1))
//     return(TRUE)
//   }
//   else
//   {
//     //return(-1)
//     return(FALSE)
//   }
	
	return (arg1.size() == arg2.size());
}

bool SameLengths(std::vector<double> arg1, std::vector<int> arg2)
{
	return (arg1.size() == arg2.size());
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3)
{
	return (SameLengths(arg1, arg2) || SameLengths(arg1, arg3));
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<int> arg3)
{
	return (SameLengths(arg1, arg2) || SameLengths(arg1, arg3));
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<double> arg4)
{
	return (SameLengths(arg1, arg2) || SameLengths(arg1, arg3) || SameLengths(arg1, arg4));
}

bool SameLengths(std::vector<double> arg1, std::vector<double> arg2, std::vector<double> arg3,
                 std::vector<int> arg4)
{
	return (SameLengths(arg1, arg2) || SameLengths(arg1, arg3) || SameLengths(arg1, arg4));
}

//This utility checks that a value falls in a valid range.
bool InRange(double value, double low, double high)
{
	return (value >= low && value <= high);
}

bool InRange(std::vector<double> value, double low, double high)
{
	for (int i = 0; i < value.size(); i++)
	{
		if (value[i] < low || value[i] > high)
		{
			return false;
		}
	}

	return true;
}

//Check if a value is from 0 to 1, a valid proportion.
bool ValidProportion(double value)
{
  return InRange(value, 0, 1);
}

bool ValidProportion(std::vector<double> value)
{
  return InRange(value, 0, 1);
}

//Calculate the sum of a variable array of the form X_ij by the specified live/dead class:
//This is a draft.  It could return the sum (an array) for each class rather than specifying one.
//C++ only.
//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat)
//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, FuelCategory liveDeadCat)
double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat)
{
	double sum = 0;//Return value.
	
	for (int k = 0; k < x_ij.size(); k++)
	{
		if (liveDead[k] == liveDeadCat)
		{
			sum += x_ij[k];
		}
	}
	
	return sum;
}

//Compare two floating point values for near/effective equality:
//Note: I'm not sure what the default should be for the precision of this comparison (see header).
//C++ only.
bool FloatCompare(double val1, double val2, double precision)//Or epsilon?
{
	if (std::fabs(val1 - val2) < precision)
	{
		return true;
	}
	else
	{
		return false;
	}
}

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
