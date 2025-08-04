/***************************************************************************************************
FireweedRAFireSpread.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/12/2024
Reference: Proj. 11 Exp. 14

Description:--------------------------------------------------------------------------------------
	This file is part of the Fireweed fire code library.  It contains a C++ implementation of the
Rothermel fire spread model (Rothermel 1972) with the modifications of Albini (Albini 1976).

	The Rothermel/Albini equations have been implemented in a form as close to the originals as was
possible.  The modifications of Albini 1972 are considered official parts of the Rothermel fire
spread model and all equations include these modifications where they apply.  The Rothermel model
consists of a nested set of equations. The library presents functions for each equation with useful
output.  Some functions are very simple and could be inlined but a focus on clarity and modularity
have been prioritized over compactness.

References:-----------------------------------------------------------------------------------------
Richard C. Rothermel.
A mathematical model for predicting fire spread in wildland fuels.
Res. Pap. INT-115. Ogden, UT: U.S. Department of Agriculture, Intermountain Forest and Range
Experiment Station. 40 p. 1972.
	This is the original article describing the Rothermel fire spread model.

Frank A. Albini.
Computer-based models of wildland fire behavior: a user's manual.
Intermountain Forest and Range Experiment Station, Forest Service, U.S. Department of Agriculture.
1976.
	This paper documents the changes that Albini made to the equations of Rothermel 1972 as well as
a set of functions implementing the equations in FORTRAN (sometimes referred to collectively as
FIREMOD, the name of the spread rate routine).  The code here was developed directly from the
equations in the paper.  I have not been able to find the FORTRAN code itself.

Patricia L. Andrews, Miguel G. Cruz, and Richard C. Rothermel.
Examination of the wind speed limit function in the Rothermel surface fire spread model.
International Journal of Wildland Fire 22(7): 959-69, 2013. http://dx.doi.org/10.1071/WF12122
	This review summarizes and contextualizes the equations of Rothermel and Albini as well as
related work and was an important reference for preparing this code.

Ralph Wilson.
Reformulation of forest fire spread equations in SI units.
Research Note INT-292.  U.S. Department of Agriculture, Forest Service, Intermountain Range and
Forest Experiment Station. Ogden, UT. 5 pages, 1980. https://doi.org/10.2737/INT-RN-292
	This report provides SI conversions of some of the spread and related equations.	These were used
to check the conversions performed here.

Anderson, Hal E.
Heat transfer and fire spread.
Res. Pap. INT-69. Ogden, UT: U.S. Department of Agriculture, Forest Service, Intermountain Forest
and Range Experiment Station. 20 p. 1969.
	This report is the source of the residency time calculation used here.

George M. Byram.
Combustion of forest fuels.
In Forest fire: control and use. Davis, K. P. editor. Pages 61-89. New York, NY: McGraw-Hill. 1959.
	This chapter defines Byram's fireline intensity and Byram's equaiton for flame length.

Notation:-------------------------------------------------------------------------------------------
	The equations from these papers contain mathematical notation and variables with characters that
cannot be directly represented in R/C++.  To represent the variables the following translations
were used.

Greek letters in variable names:
	Many model variables use Greek characters.  In most cases these are translated in the code using
their English phonetic names with the case indicating the case of the character, e.g. β -> Beta and
σ -> sigma.  In a few cases Greek variable names have been changed to abbreviations or descriptive
names.  Greek is used where the original equations are shown in the comments.

Subscripts:
	Variables with subscripts are represented with underscores, e.g. Ab (A sub b) -> A_b.  A number
of variables have two levels of subscript, the second representing fuel type indexes (i and j, see
fuel classes below).  These are represented with underscores as well, e.g (Ab)ij ((A sub b) sub ij)
-> A_b_ij.
Note: This notation is used throughout the code but is not fully consistent in comments yet.

Diacritical marks:
	In Rothermel 1972 some variables for the heterogeneous fuels equations are marked with either
bars to indicate a mean (across all fuel classes) or tildes for characteristic values of a fuel
category (live/dead).  Most reprints ignore these.  We have left them out in most cases although
a couple variables of the form x_bar are used.

	For further information consult the accompanying variable notation cross-reference document.
This contains all the input parameters and variables, many of which are common to multiple
functions.

Variable notes:
 - The surface-area-to-volume ratio for fuels is notated as σ (sigma) and abbreviated as sav or
SAV in different places in the papers.  We use SAV in the code.
 - It is unclear if fuel loading is w0 or wo.  In Rothermel 1972 it is not clear and in reprints it
varies.  We use w sub o (w_o).
 - Fuel load is specified as lb/ft^2 in the equations but most table of values use tons/acre.
Conversion needs to be done prior to passing in the data.  Use tonsPerAcToLbPerSqFt.
 - Total mineral content is occasionally notated S sub t rather than S sub T (e.g. Rothermel
1972, pg. 36, Table 1).  We use S sub T (S_T).
 - Slope steepness is notated as tan ϕ in the original model.  We use slopeSteepness instead.
Slope steepness is a fraction, not degrees.

Fuel classes:---------------------------------------------------------------------------------------
	Heterogeneous fuels fuel types are distinguished with subscripts i = 1 to m categories
(live vs. dead) and j = 1 to n fuel size classes.  j = 1 for dead fuels and j = 2 for live fuels.
All the standard fuels have three dead fuel size classes and two live classes (herbaceous vs.
woody).  Therefore storing fuel class properties in a matrix / 2D array would result in an
empty/undefined position in the live row.  We avoid this by representing fuel properties in the
code as vectors.  Dead and live size classes are stored contiguously and we maintain a liveDead
vector that holds the live/dead status for each vector index.  Where the original paper operations
operate over indexes i and j we use index k, where k = 1 to (n + m) consistently to iterate over
all positions.

Units:----------------------------------------------------------------------------------------------
	The original Rothermel/Albini model equations used United States customary units.  The units for
inputs and outputs are given in the comments for each function.  Metric/SI conversions have been
added for all functions that need them.  The units to use can be specified explicitly for some
functions, including the main spread rate functions.  If a set of intermediate calculations are to
be performed SetModelUnits() can be called at the start of the session.  The valid unit specifiers
are "US" (USCU is too hard to remember) and "Metric".  See below for more on constants and
functions used to manage units in the code.
[Note: Not entirely accurate.  See below!!!!!]
***************************************************************************************************/

#include "FireweedRAFireSpread.h"
#include "FireweedMessaging.h"
#include "FireweedStringUtils.h"
#include "FireweedUtils.h"

//Globals:------------------------------------------------------------------------------------------

/** Specify the default units to use.  The default is United States customary units.
 * This should not be set directly.  [Use SetModelUnits().]
 * 
 * Note*: This global is not integrated in the same way it is in the R code.  It establishes default
 * units for functions but the code does not curently change.  The R code provides SetModelUnits()
 * to set the units for a session.  It is unclear how useful this would be for the C++ use cases.
 */
UnitsType ModelUnits = US;

//Code:---------------------------------------------------------------------------------------------

//Bulk Density:--------------------------------------------------------------------------------------
//	The bulk density is the mass/wt. of oven dry surface fuel per volume of fuel bed (fuel mass per
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
//	The heterogeneous fuel version of the spread equation requires a mean bulk density for the 
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
//	The packing ratio (beta) is the fraction of the (surface) fuel bed volume occupied by fuel, aka
//compactness.
//
//Rothermel 1972 equation 31:
//β = ρb/ρp
//
//Input variables / parameters:
//rho_b = Fuel array bulk density (lb/ft^3 | kg/m^3).
//rho_p = Fuel particle density (lb/ft^3 | kg/m^3).
//	For the 53 standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
//
//Output units: Dimensionless ratio.
//Input units cancel out.  No metric conversion needed.
double PackingRatio(double rho_b, double rho_p)
{
	return (rho_b / rho_p);
}

//Mean Packing Ratio:
//For heterogeneous fuel beds the mean packing ratio (beta_bar) must be calculated.
//
//Rothermel 1972 equation 74:
//β = 1/δ Σi Σj (wo)ij/(ρp)ij
//For i = 1 to m fuel categories (live vs. dead) and j = 1 to n fuel size classes.
//The original notation includes from and to sum subscripts and bars over beta, rho, and w.
//
//Input variables / parameters:
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//rho_p_ij = Fuel particle density for each fuel type (lb/ft^3 | kg/m^3).
//	For the 53 standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
//fuelBedDepth = Fuel bed depth, AKA delta (ft | m).
//
//Output units: Dimensionless ratio (scalar)
//Input units cancel out.
double MeanPackingRatio(std::vector<double> w_o_ij, std::vector<double> rho_p_ij, double fuelBedDepth)
{
	int numLoadings, numDensities;
	double meanPackingRatio;//Return value.
	double x;//Intermediate accumulator.

	//Parameter checking:
	numLoadings = w_o_ij.size();
	numDensities = rho_p_ij.size();

	//If only one particle density is provided assume that it is the same for all fuel classes:
	if (numDensities == 1)
	{
		rho_p_ij.resize(numLoadings, rho_p_ij[0]);
	}
	else//Otherwise one should be provided for each fuel class.
	{
		if (numDensities != numLoadings)
		{
			Stop("The number of fuel loadings and particle densities do not match.");
		}
	}

	//Calculate w_o / rho_p for each fuel element and accumulate:
	//x = w_o_ij / rho_p_ij
	for (int i = 0; i < numLoadings; i++)
	{
		x += w_o_ij[i] / rho_p_ij[i];
	}

	meanPackingRatio = x / fuelBedDepth;//AKA beta_bar
	return meanPackingRatio;
}

//Optimum Packing Ratio:
// This is the packing ratio (compactness) at which the maximum reaction intensity will occur.
//
//Rothermel 1972 equations 37,69:
//βop = 3.348(σ)^-0.8189
//
//Input variables / parameters:
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//	For heterogeneous fuels the SAV of the fuel bed / complex is used.
//
//Output units: Dimensionless ratio
double OptimumPackingRatio(double cSAV, UnitsType units)
{
	double optPackingRatio;//Return value.

	if (units == US)
	{
		optPackingRatio = 3.348 * pow(cSAV, -0.8189);
	}
	else// if (units == "Metric")
	{
		//Wilson 1980 gives:
		//optPackingRatio = 0.219685 * SAV^-0.8189
		//Tests show that this is pretty good but I was able to calculate a better conversion as follows:
		optPackingRatio = 0.2039509 * pow(cSAV, -0.8189);
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
//	dead or live fuel category.
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
	A_T = A_i[Dead] + A_i[Live];//Single scalar value.

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
	wts.f_i[Dead] = A_i[Dead] / A_T;
	wts.f_i[Live] = A_i[Live] / A_T;

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
		double subclassTotal[6] = {};//Explicitly initialize to zero for each iteration of the loop.
		
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
	//when we reassemble the weights.  FloatCompare() handles this problem.

	//The dead fuel components of f_ij should always sum to 1:
	if (!FloatCompare(SumByFuelCat(wts.f_ij, liveDead, Dead), 1))
	{
		Stop("f_ij dead fuels do not sum to 1.");
	}

	//The live fuel components of f_ij will sum to 1 if present or 0 if not present:
	if (!(FloatCompare(SumByFuelCat(wts.f_ij, liveDead, Live), 0) ||
			  FloatCompare(SumByFuelCat(wts.f_ij, liveDead, Live), 1)))
	{
		Stop("Invalid f_ij weights for live fuels.");
	}

	/*The dead fuel components of g_ij should sum to 1, except when curing has been applied.  When
	curing results in the addition of a dead herbaceous class the weights will still sum to 1 if the
	dead SAVs all fall in different bins.  However, frequently the dead herbaceous SAV is close to
	the 1hr class SAV and they both fall in the same bin.  In such cases they get the same weight.
	If they are the only fuels present they will both have a weight of 1 for a total of 2.  If there
	is another class present the total will be greater than 1 but less than 2.
	
	If we had access the fuel model object here we could check if curing had been applied.  Since we
	don't we have to be clever.  The standard fuel models all have five fuel classes.  Having six
	classes is pretty good evidence of curing, but we need to be tolerant of custom fuel models.
	We currently always put the dead herbaceous class at the second position (1,2).  If we compare
	the SAVs where we expect the herbaceous classes to be that is a pretty good check.
	It might be safer to pass in the curing status.*/
	//int oneHrPos = 0;//Might be good to abstract these assumed indexes!!!!!
	//int herbPos = 0;
	if (SAV_ij[FuelClassIndex(liveDead, Dead, 1)] == SAV_ij[FuelClassIndex(liveDead, Live, 0)])//Implied cured:
	{
		if (!FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Dead), 1.0))
		{
			//The 1hr fuel (0,0) and dead herbaceous (0,1) classes should have the same weights,
			//unless the curing is 0.  In that case w_o_01 is 0 and the total will be 1:
			if (wts.g_ij[FuelClassIndex(liveDead, Dead, 0)] !=
			    wts.g_ij[FuelClassIndex(liveDead, Dead, 1)])
			{
				//Stop("g_ij 1hr and dead herbaceous fuels do not have the same weights.");
				Stop("g_ij 1hr and dead herbaceous fuels do not have the same weights: " +
				     VectorToStr(wts.g_ij));
			}

			//Excluding the dead herbaceous the dead fuels should sum to 1:
			if (!FloatCompare((SumByFuelCat(wts.g_ij, liveDead, Dead) -
			                   wts.g_ij[FuelClassIndex(liveDead, Dead, 1)]), 1.0))
			{
				//Stop("g_ij dead fuels do not have values expected for implied cured status.");
				Stop("g_ij dead fuels do not have values expected for implied cured status: " +
				     VectorToStr(wts.g_ij));
			}
		}
	}
	else//Uncured:
	{
		if (!FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Dead), 1.0))
		{
			Stop("g_ij dead fuels do not sum to 1.");
		}
	}

	//The dead fuel components of g_ij should always sum to 1:
	//Actually see above?????
	//if (!FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Dead), 1))
	//{
	//	Stop("g_ij dead fuels do not sum to 1.");
	//}

	//For static models the live fuel components of f_ij will sum to 1 if present or 0 if not present.
	//However, for dynamic fuel models both live classes may be have values of 0 or 1, so sums of 0, 1,
	//and 2 are possible:
	if (!(FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Live), 0.0) ||
		  FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Live), 1.0) ||
		  FloatCompare(SumByFuelCat(wts.g_ij, liveDead, Live), 2.0)))
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

	return (CalcWeightings(SAV_ij, w_o_ij, rho_p_ij, liveDead, units));
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
//	For heterogeneous fuels a SAV for the entire fuel bed must be calculated.  This is frequently
//referred to as the characteristic SAV (cSAV).  It is a weighted average of the fuel component SAVs.
//
//Input variables / parameters:
//SAV_ij = Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
//f_ij = Weighting factors for each fuel type (dimensionless).
//f_i = Weighting factors for each fuel live/dead category (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//	dead or live fuel category.
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
	double SAV_i[2] = {0, 0};;//SAV by live / dead category.
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
		SAV_i[liveDead[k]] += f_ij[k] * SAV_ij[k];
	}

	//Sum the live and dead components to get the final value:
	//Rothermel 1972 equation 71:
	//σ = Σi fiσi (~ over sigma and sigma sub i in original)
	fuelBedSAV = (f_i[Dead] * SAV_i[Dead]) + (f_i[Live] * SAV_i[Live]);

	//Value checking:
	if (fuelBedSAV <= 0.0)
	{
		Warning("FuelBedSAV() output is invalid.");
		Msg.Log("fuelBedSAV =", fuelBedSAV);
		Msg.Log("FuelBedSAV() inputs:");
		Msg.Log("SAV_ij =", SAV_ij);
		Msg.Log("f_ij =", f_ij);
		Msg.Log("f_i =", f_i);
		Msg.Log("liveDead =", liveDead);
	}
	//Checking positve values wouold require knowing the units and a reasonable value ranges.

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
//	For all standard fuel models this is 5.55% (0.0555).
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
//	For heterogeneous fuel beds the moisture damping coefficient is calculated for each fuel category
//(live/dead).
//
//Input variables / parameters:
//M_f_ij = Fuel moisture content for for each fuel type (fraction: water weight/dry fuel weight).
//M_x_i = Moisture of extinction each fuel category (fraction: water weight/dry fuel weight).
//f_ij = Weighting factors for each fuel type (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//	dead or live fuel category.
//
//Output units: Dimensionless coefficient (array length 2)
//Input units cancel out.  No metric conversion needed.
std::vector <double> MoistureDampingCoefficient_Het(std::vector <double> M_f_ij,
                                                    std::vector <double> M_x_i,
                                                    std::vector <double> f_ij,
                                                    std::vector <int> liveDead)
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
//	The live fuel moisture of extinction determines if live fuels will burn and contribute to the
//heat source term in the calculations.  The live fuel moisture of extinction is calculated from the
//dead fuel moisture of extinction and in relation to the ratio of live to dead fuels. The value can
//be quite variable from near zero to over 700 percent (see Andrews 2018 section 5.3.2.2).
//
//Input variables / parameters:
//M_f_ij = Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
//M_x_1 = Dead fuel moisture of extinction (fraction: water weight/dry fuel weight).
//w_o_ij = An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
//SAV_ij = Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
//liveDead = An array indicating if each index in each of the other input variables represents a
//	dead or live fuel category.
//
//Output units: fraction, water weight/dry fuel weight
double LiveFuelMoistureOfExtinction(std::vector <double> M_f_ij, double M_x_1,
                                    std::vector <double> w_o_ij, std::vector <double> SAV_ij,
                                    std::vector <int> liveDead, UnitsType units)
{
	int numFuelTypes;
	double deadSum = 0;
	double liveSum = 0;
	double W;//Dead / live ratio.
	double common;
	double top = 0;//The numerator sum.
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
		if (liveDead[k] == Dead)
		{
			deadSum += w_o_ij[k] * exp(-138 / SAV_ij[k]);
		}
		else//(liveDead[k] == Live)
		{
			liveSum += w_o_ij[k] * exp(-500 / SAV_ij[k]);
		}
	}

	//If the loading for the live fuel categories are all 0 or live categories are missing liveSum will
	//be zero.  The ratio W will be therefore also be zero.  This in turn will result in M_x_2 = M_x_1.
	//While not conceptual meaningful this has no has no mathematical consequence downstream.  Forcing
	//the value to NA or 0 would cause mathematical problems downstream.

	W = deadSum / liveSum;//Unitless ratio.

	//Calculate fine dead fuel moisture as:
	//Albini 1976 pg. 16:
	//Mf,dead = Σj(Mf)1j(wo)1jexp(–138/σ1j) / Σj(wo)1jexp(–138/σ1j)
	//The denominator sum is the same as deadSum above.
	for (int k = 0; k < numFuelTypes; k++)
	{
		if (liveDead[k] == Dead)
		{
			top += M_f_ij[k] * w_o_ij[k] * exp(-138 / SAV_ij[k]);
		}
	}

	M_f_dead = top / deadSum;//Moisture fraction / unitless.

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
//	For all standard fuel models this is 1% (0.01).
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
//	For heterogeneous fuels the mineral damping coefficient is calculated for each fuel category
//(live/dead).
//
//Input variables / parameters:
//S_e_ij = Effective mineral content for each fuel type (unitless fraction:
//	(mineral mass – mass silica) / total dry mass).
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
	for (int k = 0; k < numFuelTypes; k++)
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
//	(dimensionless).
//slopeSteepness = Slope steepness maximum (unitless fraction: vertical rise / horizontal distance).
//	AKA tan ϕ.
//
//Output units: Dimensionless adjustment factor
//Inputs are fractions which do not change with units.  Mo metric conversion required.
double SlopeFactor(double packingRatio, double slopeSteepness)
{
	double phi_s;

	//Parameter checking:
	CheckSlope(slopeSteepness);

	phi_s = 5.275 * pow(packingRatio, -0.3) * pow(slopeSteepness, 2);
	return phi_s;
}

//Wind factor:
//Dimensionless multiplier that accounts for the effect of wind speed on spread behavior
//(propagating flux ratio specifically).  Same for homegenous and heterogeneous fuels.
//
//Input variables / parameters:
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//packingRatio = Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//	(dimensionless).
//optPackingRatio = Optimum packing ratio (dimensionless).
//	 Note: optimum packing ratio is a function of SAV.
//U = Wind speed at midflame height (ft/min | m/min).
//
//Output units: Dimensionless
//
//Note: It is not possible to calculate if a wind limit is indicated internal to this function
//and not all authors agree that a wind limit should be used.  U should be capped, if deemed
//appropriate prior to passing it in to this function.
double WindFactor(double cSAV, double packingRatio, double optPackingRatio, double U,
                  UnitsType units)
{
	double C, B, E;//Partial terms.
	double phi_w;//Return value.

	C = WindFactorC(cSAV, units);
	B = WindFactorB(cSAV, units);
	E = WindFactorE(cSAV, units);

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

double WindFactorC(double cSAV, UnitsType units)
{
	double C;//Return value.

	if (units == US)
	{
		//C = unnamed term
		//Rothermel 1972 equation 48,82:
		//C = 7.47exp(-0.133σ^0.55)
		C = 7.47 * exp(-0.133 * pow(cSAV, 0.55));
	}
	else
	{
		C = 7.47 * exp(-0.8710837 * pow(cSAV, 0.55));//-0.133 * cmPerFt^0.55 = -0.8710837
		//Wilson 1980 uses:
		//C = 7.47 * exp(-0.8711 * SAV^0.55)
	}
	return C;
}

double WindFactorB(double cSAV, UnitsType units)
{
	double B;//Return value.

	if (units == US)
	{
		//B = unnamed term
		//Rothermel 1972 equation 49,83:
		//B = 0.02526σ^0.54
		B = 0.02526 * pow(cSAV, 0.54);
	}
	else
	{
		B = 0.1598827 * pow(cSAV, 0.54);//0.02526 * cmPerFt^0.54 = 0.1598827
		//Wilson 1980 uses:
		//B = 0.15988 * SAV^0.54
	}
	return B;
}

double WindFactorE(double cSAV, UnitsType units)
{
	double E;//Return value.

	if (units == US)
	{
		//E = unnamed term
		//Rothermel 1972 equation 50,84:
		//E = 0.715exp(-3.59×10^-4 σ)
		E = 0.715 * exp(-0.000359 * cSAV);
	}
	else
	{
		E = 0.715 * exp(-0.01094232 * cSAV);//-0.000359 * cmPerFt = -0.01094232
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
double WindFactorA(double cSAV, double packingRatio, double optPackingRatio, UnitsType units)
{
	double C, E, A;//Partial terms.

	C = WindFactorC(cSAV, units);
	E = WindFactorE(cSAV, units);
	A = C * pow((packingRatio / optPackingRatio), -E);
	return A;
}

//Wind limit:
//	The "wind limit" or "maximum reliable wind" is used to limit the effect of wind on the fire
//spread rate as wind speed gets high.  It caps the wind speed at a value that is a function of the
//reaction intensity.
//	 There is not agreement on whether the wind limit should be used.  Albini chose to not use it,
//but his code reports if the limit was reached (Albini 1976, pg 26).  More recent work finds the
//original calculation to be flawed and presents an alternate formulation from (Andrews et. al 2013).
//However, They conclude that in general neither should be used.  They state a better alternative is
//to cap the spread rate at the “effective wind speed”.
//	We implement the original formulation as an option to be able to reproduce results that do use
//the wind limit.
//
//Input variables / parameters:
//U = Wind speed at midflame height (ft/min | m/min).
//I_R = Reaction intensity (Btu/ft^2/min | kJ/m^2/min).
//
//Output units: adjusted wind speed (U) at midflame height (ft/min | m/min)
double WindLimit(double U, double I_R, UnitsType units)
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
//	This is a measure of the optimum (potential) fuel consumption rate (fire efficiency / reaction
//time).  The 'optimum' rate is the ideal rate that would occur for alpha cellulose in the absence
//of minerals and moisture.
//Notation: Γ' (Gamma prime)
//
//	This version includes the Albini 1976 modification.
//
//Input variables / parameters:
//β = (Mean) Packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
//packingRatio = (Mean) Packing ratio (β), the fraction of the fuel bed volume occupied by fuel
//	(dimensionless). For heterogeneous fuels the mean packing ratio is passed in.
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//	For heterogeneous fuels the SAV of the fuel bed / complex is used.
//
//Output units: min^-1
double OptimumReactionVelocity(double packingRatio, double cSAV, UnitsType units)
{
	double optPackingRatio;
	double A;//Intermediate.
	double GammaPrimeMax;//Intermediate.
	double GammaPrime;//Return value.

	//Calculate the maximum reaction velocity (min^-1):
	//This is the rate for moisture free fuel with mineral composition of alpha cellulose.
	//Rothermel 1972 equations 36,68:
	//Γ'max = σ^1.5/(495 + 0.0594σ^1.5)
	if (units == US)
	{
		GammaPrimeMax = pow(cSAV, 1.5) / (495 + 0.0594 * pow(cSAV, 1.5));
		//Or equivalently:
		//GammaPrimeMax = 1 / (0.0594 * 495 / cSAV^1.5)
	}
	else
	{
		GammaPrimeMax = 1 / (0.0594 + 2.941594 / pow(cSAV, 1.5));
		//Wilson 1980 uses:
		//GammaPrimeMax = (0.0591 + 2.926 * SAVcm^-1.5)^-1 = 1 / (0.0591 + 2.926 / SAVcm^1.5)
	}

	optPackingRatio = OptimumPackingRatio(cSAV, units);

	//"Arbitrary" variable:
	//Albini 1976 pg. 15:
	//A = 133σ^-0.7913
	if (units == US)
	{
		A = 133 * pow(cSAV, -0.7913);
	}
	else
	{
		A = 8.903291 * pow(cSAV, -0.7913);
		//Wilson 1980 uses:
		//A = 8.9033 * SAV^-0.7913
	}

	//These are combined to produce the optimal reaction velocity (min^-1):
	//Rothermel 1972 equation 38:
	//Γ' = Γ'max(β/βop)^A exp[A(1 - β/βop)]
	GammaPrime = GammaPrimeMax * pow((packingRatio/optPackingRatio), A) *
		exp(A * (1 - packingRatio/optPackingRatio));

	return GammaPrime;
}

//Live / Dead Heat Content:
//	Calculate the weighted (low?) heat content of the live / dead fuel categories.
//Only used by reaction intensity calculation.
//
//Input variables / parameters:
//h_ij = Heat content of the fuel types (Btu/lb | kJ/kg).
//f_ij = Weighting factors for each fuel type (dimensionless).
//liveDead = An array indicating if each index in each of the other input variables represents a
//	dead or live fuel category.
//
//Output units: btu/lb | kJ/kg
//Whatever units are input, the same will come out.  No unit conversions needed.
std::vector <double> LiveDeadHeatContent(std::vector <double> h_ij, std::vector <double> f_ij,
                                         std::vector <int> liveDead)
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
//	The reaction intensity (I_R) is the total energy released by the fire front in Btu/ft^2/min in
//all forms (radiation, conduction, and convection).
//	This it not the same as fireline intensity!

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
double ReactionIntensity_Homo(double GammaPrime, double w_n, double h, double eta_M, double eta_s)
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

	I_R = GammaPrime * ((w_n_i[Dead] * h_i[Dead] * eta_M_i[Dead] * eta_s_i[Dead]) +
	                    (w_n_i[Live] * h_i[Live] * eta_M_i[Live] * eta_s_i[Live]));

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
//	(dimensionless). For heterogeneous fuels the mean packing ratio is passed in.
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//	For heterogeneous fuels the fuel bed level SAV is used.
//
//Output units: Dimensionless proportion
double PropagatingFluxRatio(double packingRatio, double cSAV, UnitsType units)
{
	double xi = 0;//Output

	if (units == US)
	{
	  xi = pow((192 + 0.2595 * cSAV), -1) * exp((0.792 + 0.681 * pow(cSAV, 0.5)) * (packingRatio + 0.1));
	}
	else
	{
		xi = pow((192 + 7.90956 * cSAV), -1) * exp((0.792 + 3.759712 * pow(cSAV, 0.5)) * (packingRatio + 0.1));
		//Wilson 1980 uses:
		//xi = (192 + 7.9095 * SAV)^-1 * exp((0.792 + 3.7597 * SAV^0.5) * (packingRatio + 0.1))
	}

	return xi;
}

//This is a wrapper for PropagatingFluxRatio() that allows it to be called from R:
//The primary reason for providing this wrapper is for verification testing.
extern "C" void PropagatingFluxRatioR(const double* packingRatio, const double* cSAV, double* xi)
{
	*xi = PropagatingFluxRatio(*packingRatio, *cSAV);
}

//Heat Sink Components:-----------------------------------------------------------------------------

//Effective Heating Number:
//	This represents the proportion of a fuel type that is heated to ignition temperature in advance
//of the fire front.  It is a function of SAV.  For example, with the same heating fine fuels might
//be brought fully to combustion while for a large stick only surface would be dried and heated to
//burn.
// This function is only used in the homogeneous fuel calculations.
//
//Rothermel 1972 equation 14
//ε = exp(-138/σ)
//
//Input variables / parameters:
//SAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//
//Units: Dimensionless.
double EffectiveHeatingNumber(double SAV, UnitsType units)
{
	double epsilon;//Return value.

	if (units == US)
	{
		epsilon = exp(-138/SAV);
	}
	else// if (units == Metric)
	{
		epsilon = exp(-4.527559/SAV);//-138 * cmPerFt = -4.527559.  Wilson 1980 uses −4.528.
	}

	return epsilon;
}

//Heat Of Preignition:
//
//Rothermel 1972 equations 12,78:
//Qig = 250 + 1,116Mf
//
//Input variables / parameters:
//M_f = Fuel moisture content (fraction: water weight/dry fuel weight).
//
//Output units: btu/lb or kJ/kg
//
//For heterogeneous fuels this is calculated for each fuel type.
double HeatOfPreignition(double M_f, UnitsType units)
{
	double Q_ig;//Return value.

	//Validity checking:
	if (!InRange(M_f, 0, 3.5))
	{
		Stop("Suspect moisture content (M_f).");
	}

	if (units == US)
	{
		Q_ig = 250 + 1116 * M_f;
	}
	else// if (units == "Metric")
	{
		Q_ig = 581.1114 + 2594.081 * M_f;//Constants * kJPerBtu / kgPerLb
		//Wilson 1980 uses:
		//Qig = 581 + 2594 * M_f
		//This implies Wilson was using the thermochemical BTU conversion which is ~1,054.35 J/BTU.
	}

	return Q_ig;
}

//A vector aware version of HeatOfPreignition().
std::vector <double> HeatOfPreignition(std::vector <double> M_f_ij, UnitsType units)
{
	int numFuelTypes;
	std::vector <double> Q_ig_ij(M_f_ij.size(), 0);//Return value.

	numFuelTypes = M_f_ij.size();

	for (int i = 0; i < numFuelTypes; i++)
	{
		Q_ig_ij[i] = HeatOfPreignition(M_f_ij[i]);
	}

	return Q_ig_ij;
}

//Spread Rate Calculations:-------------------------------------------------------------------------

/** Albini 1976 modified Rothermel spread model for homogeneous fuels:
 *
 * 	Calculate the steady state spread rate for surface fuels and environmental conditions passed in.
 *
 * @par Input variables / parameters:
 * There are 11 input variables in total (see Andrews 2018 table 11), 4 fuel array characteristics,
 * 4 fuel particle characteristics, and three environmental.
 *
 * Fuel array:
 * @param SAV Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
 * @param w_o Oven dry fuel load (lb/ft^2 | kg/m^2).  This includes combustible and mineral fractions.
 * @param fuelBedDepth Fuel bed depth, AKA delta (ft | m).
 * @param M_x Moisture of extinction (fraction: water weight/dry fuel weight).
 *
 * Environmental:
 * @param M_f Fuel moisture content (fraction: water weight/dry fuel weight).
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 *
 * Fuel particle properties: 
 * @param h     Heat content of the fuel type (Btu/lb | kJ/kg).
 *              All the 53 standard fuel models use 8,000 Btu/lb.
 * @param S_T   Total mineral content (unitless fraction: mineral mass / total dry mass).
 *              For all standard fuel models this is 5.55% (0.0555).
 * @param S_e   Effective mineral content (unitless fraction: (mineral mass – mass silica) / total dry mass).
 *              For all standard fuel models this is 1% (0.01).
 * @param rho_p Fuel particle density (lb/ft^3 | kg/m^3).
 *              All the 53 standard fuel models use 32 lb/ft^3.
 *
 * Settings:
 * @param useWindLimit Use the wind limit calculation or not.
 * @param units Specify the class of units for the inputs.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 */
double SpreadRateRothermelAlbini_Homo(double SAV, double w_o, double fuelBedDepth, double M_x,
                                      double M_f, double U, double slopeSteepness,
                                      double heatContent,//h
                                      double S_T, double S_e,
                                      double rho_p,
                                      bool useWindLimit,
                                      UnitsType units,
                                      bool debug)
{
	SpreadCalcs calcs;

	calcs = SpreadCalcsRothermelAlbini_Homo(SAV, w_o, fuelBedDepth, M_x, M_f, U, slopeSteepness,
	                                        heatContent, S_T, S_e, rho_p, useWindLimit, units,
	                                        debug);

	return calcs.R;
}

/*This is a wrapper for SpreadRateRothermelAlbini_Homo() that allows it to be called from R via .C():
The parameters are the same as SpreadRateRothermelAlbini_Homo() except:
 - The useWindLimit and debug arguments are type int.  This is port of the .C() interface.  Logical
should be used on the R side.
 - units is passed in as an integer with 1 = US and 2 = metric.
 - The spread rate is returned in the additional argument R.*/
extern "C" void SpreadRateRothermelAlbini_HomoR(const double* SAV, const double* w_o,
                                                const double* fuelBedDepth, const double* M_x,
                                                const double* M_f, const double* U,
                                                const double* slopeSteepness, const double* heatContent,//h
                                                const double* S_T, const double* S_e,
                                                const double* rho_p, const int* useWindLimit,
                                                const int* units, const int* debug, double* R)
{
	bool useWindLimitBool;
	UnitsType cUnits;
	bool debugBool;

	//R logicals are passed as integer values.  This should be seamless to the user:
	if (*useWindLimit == 0)
	{
		useWindLimitBool = false;
	}
	else if (*useWindLimit == 1)
	{
		useWindLimitBool = true;
	}
	else//The NA value is possible it NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for useWindLimit.");
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
	
	if (*debug == 0)
	{
		debugBool = false;
	}
	else if (*debug == 1)
	{
		debugBool = true;
	}
	else//The NA value is possible it NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for debugBool.");
	}
	
	*R = SpreadRateRothermelAlbini_Homo(*SAV, *w_o, *fuelBedDepth, *M_x, *M_f, *U, *slopeSteepness,
                                        *heatContent, *S_T, *S_e, *rho_p, useWindLimitBool, cUnits,
                                        debugBool);
}

/** Albini 1976 modified Rothermel spread model for homogeneous fuels, fuel model version:
 *
 * Calculate the steady state spread rate for surface fuels and environmental conditions passed in.
 *
 * @par Input variables / parameters:
 * This version greatly reduces the number of parameters by taking a fuel model object that contains
 * the properties of the fuel bed.
 *
 * @param fuelModel A fuel model object.
 *
 * Environmental:
 * @param M_f Fuel moisture content (fraction: water weight/dry fuel weight).
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 *
 * Settings: (optional)
 * @param useWindLimit Use the wind limit calculation or not.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 */
double SpreadRateRothermelAlbini_Homo(FuelModel fuelModel,
                                      double M_f, double U, double slopeSteepness,
                                      bool useWindLimit,
                                      bool debug)
{
	SpreadCalcs calcs;

	calcs = SpreadCalcsRothermelAlbini_Homo(fuelModel, M_f, U, slopeSteepness, useWindLimit, debug);

	return calcs.R;
}

/** Albini 1976 modified Rothermel spread model for homogeneous fuels:
 *
 * Calculate the steady state spread rate for surface fuels and environmental conditions passed in.
 * This variant returns many of the component calculated values in the spread rate calculation.
 * The function parameters are the same as SpreadRateRothermelAlbini_Homo().
 *
 * @returns A SpreadCalcs object.
 *
 * @note To avoid code repetition this function is called by the spread rate only functions.
 */
SpreadCalcs SpreadCalcsRothermelAlbini_Homo(double SAV, double w_o, double fuelBedDepth, double M_x,
                                            double M_f, double U, double slopeSteepness,
                                            double heatContent,
                                            double S_T, double S_e,
                                            double rho_p,
                                            bool useWindLimit,
                                            UnitsType units,
                                            bool debug)
{
	double rho_b, packingRatio, optPackingRatio;//Bulk density and packing ratio.
	double GammaPrime, w_n, eta_M, eta_s, I_R;//Reaction intensity & intermediates.
	double xi, phi_s, phi_w, epsilon, Q_ig;//Other intermediate terms.
	double heatSource, heatSink;
	double R;//Spread rate.
	SpreadCalcs calcs;//Return value.

	//Up front calculations:
	//The bulk density is needed to calculate the packing ratio and therefore is used in the numerator
	//and denominator.
	rho_b = BulkDensity(w_o, fuelBedDepth);

	packingRatio = PackingRatio(rho_b, rho_p);
	optPackingRatio = OptimumPackingRatio(SAV, units);

	//The heat source term (numerator) represents the heat flux from the fire front to the fuel in
	//front of it (AKA propagating flux) in BTU/min/ft^2 | kW/m^2:
	//Numerator of Rothermel 1972 equation 52:
	//IR𝜉(1 + 𝜙w + 𝜙s)

	//Reaction intensity I_R:
	GammaPrime = OptimumReactionVelocity(packingRatio, SAV, units);
	w_n = NetFuelLoad_Homo(w_o, S_T);
	eta_M = MoistureDampingCoefficient_Homo(M_f, M_x);
	eta_s = MineralDampingCoefficient_Homo(S_e);
	I_R = ReactionIntensity_Homo(GammaPrime, w_n, heatContent, eta_M, eta_s);

	//Other terms:
	xi = PropagatingFluxRatio(packingRatio, SAV, units);
	phi_s = SlopeFactor(packingRatio, slopeSteepness);

	//Apply wind limit check:
	if (useWindLimit)
	{
		U = WindLimit(U, I_R, units);
	}

	phi_w = WindFactor(SAV, packingRatio, optPackingRatio, U, units);

	heatSource = I_R * xi * (1 + phi_w + phi_s);

	//The heat sink term (denominator) represents the energy required to ignite the fuel in Btu/ft^3 |
	//kJ/m^3:
	//Denominator of Rothermel 1972 equation 52:
	//ρbεQig

	epsilon = EffectiveHeatingNumber(SAV, units);
	Q_ig = HeatOfPreignition(M_f, units);

	heatSink = rho_b * epsilon * Q_ig;

	//Full spread calculation for homogeneous fuels:
	//Rothermel 1972 equation 52:
	//Rate of spread = heat source / heat sink
	//R = I_R𝝃(1 + 𝜙w + 𝜙s) / ρbεQig
	//R = (I_R * xi * (1 + phi_w + phi_s)) / (rho_b * epsilon * Q_ig);
	R = heatSource / heatSink;

	calcs.units = units;
	calcs.homogeneous = true;
	calcs.R = R;
	//calcs.weights can remain empty.  We could make a single element set of weights but that seems unnecessary.
	calcs.GammaPrime = GammaPrime;
	calcs.w_n_x = std::vector <double> {w_n};
	calcs.h_x = std::vector <double> {heatContent};
	calcs.eta_M_x = std::vector <double> {eta_M};
	calcs.eta_s_x = std::vector <double> {eta_s};
	calcs.I_R = I_R;
	calcs.xi = xi;
	calcs.phi_s = phi_s;
	calcs.phi_w = phi_w;
	calcs.heatSource = heatSource;
	calcs.rho_b_x = rho_b;//The mean of a single value is itself.
	calcs.epsilon = epsilon;
	calcs.Q_ig_x = std::vector <double> {Q_ig};
	calcs.heatSink = heatSink;
	calcs.cSAV = SAV;
	calcs.packingRatio = packingRatio;//The mean of a single value is itself.
	calcs.optimumPR = optPackingRatio;

	//For debugging:  Note: Could be replaced with a SpreadCalcs print function.
	if (debug)
	{
		Msg.Log("Homogeneous Spread Calc components:");
		Msg.Log("GammaPrime =", GammaPrime);
		Msg.Log("w_n =", w_n);
		Msg.Log("h =", heatContent);
		Msg.Log("eta_M =", eta_M);
		Msg.Log("eta_s =", eta_s);
		Msg.Log("I_R =", I_R);
		Msg.Log("xi =", xi);
		Msg.Log("phi_s =", phi_s);
		Msg.Log("phi_w =", phi_w);
		Msg.Log("Heat source =", heatSource);
		Msg.Log("rho_b =", rho_b);
		Msg.Log("epsilon =", epsilon);
		Msg.Log("Qig =", Q_ig);
		Msg.Log("Heat sink = ", heatSink);
	}

	return calcs;
}

/** Albini 1976 modified Rothermel spread model for homogeneous fuels, fuel model version:
 *
 * Calculate the steady state spread rate for surface fuels and environmental conditions passed in.
 * This variant returns many of the component calculated values in the spread rate calculation.
 * The function parameters are the same as fuel model version of SpreadRateRothermelAlbini_Homo().
 *
 * @par Input variables / parameters:
 * This version greatly reduces the number of parameters by taking a fuel model object that contains
 * the properties of the fuel bed.
 *
 * @param fuelModel A fuel model object.
 *
 * Environmental:
 * @param M_f Fuel moisture content (fraction: water weight/dry fuel weight).
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 *
 * Settings: (optional)
 * @param useWindLimit Use the wind limit calculation or not.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 */
SpreadCalcs SpreadCalcsRothermelAlbini_Homo(FuelModel fuelModel,
                                            double M_f, double U, double slopeSteepness,
                                            bool useWindLimit,
                                            bool debug)
{
	SpreadCalcs calcs;

	//Check the fuel model units.  Conversions could be applied:
	if (fuelModel.units == US)
	{
		if (fuelModel.w_o_Units == tonPerAc)
		{
			Stop("Loading units are in tons per acre and need to be in lb per ft^2.");
		}
		
		if (fuelModel.M_x_Units == Percent)
		{
			Stop("M_x units are in percent and need to be a fraction.");
		}
	}

	//Check that the fuel model is in fact homogeneous:
	if (fuelModel.numClasses != 1)
	{
		Warning("The fuel model is not homogeneous as expected. Ignoring all but the first fuel class.");
	}

	calcs = SpreadCalcsRothermelAlbini_Homo(fuelModel.SAV_ij[0],
	                                        fuelModel.w_o_ij[0],
	                                        fuelModel.delta,//fuelBedDepth
	                                        fuelModel.M_x_1,//M_x FW_NOTE: Restore!!!!!
	                                        M_f, U, slopeSteepness,
	                                        fuelModel.h_ij[0],
	                                        fuelModel.S_T_ij[0],
	                                        fuelModel.S_e_ij[0],
	                                        fuelModel.rho_p_ij[0],
	                                        useWindLimit,
	                                        fuelModel.units,
	                                        debug);
	return calcs;
}

/* Albini 1976 modified Rothermel spread model for heterogeneous fuels:
 *
 * @par Input variables / parameters:
 * Some of the input variables differ from the homogeneous fuels form in that they are vectors
 * rather than scalars.
 *
 * Fuel array:
 * @param SAV_ij Characteristic surface-area-to-volume ratios for each fuel type (ft^2/ft^3 | cm^2/cm^3).
 * @param w_o_ij An array of oven dry fuel load for each fuel type (lb/ft^2 | kg/m^2).
 * @param fuelBedDepth Fuel bed depth, AKA delta (ft | m).
 * @param M_x_1 Dead fuel moisture of extinction (fraction: water weight/dry fuel weight).
 * @param liveDead An array indicating if each index in each of the other input variables represents a
 * 	dead or live fuel category.  Note: This is placed later in the argument list to allow for a
 * 	default value.  If omitted it is assumed a standard fuel model with 5 fuel types is in use.
 *
 * Environmental:
 * @param M_f_ij Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 *
 * Fuel particle properties: 
 * @param h_ij     Heat content of the fuel types (Btu/lb | kJ/kg).
 *                 All the 53 standard fuel models use 8,000 Btu/lb.
 * @param S_T_ij   An array of total mineral content for each fuel type (unitless fraction).
 *                 For all standard fuel models this is 5.55% (0.0555).
 * @param S_e_ij   Effective mineral content for each fuel type (unitless fraction:
 *                 (mineral mass – mass silica) / total dry mass).  For all standard fuel models this is 1% (0.01).
 * @param rho_p_ij Fuel particle density for each fuel type (lb/ft^3 | kg/m^3).
 *                 All the 53 standard fuel models use 32 lb/ft^3.
 *
 * Settings:
 * @param useWindLimit Use the wind limit calculation or not.  Recent suggestion are that it not be used.
 * @param units Specify the class of units for the inputs.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 *
			 * This function takes a lot of arguments.  These parameters could be combined into fuel model
			 * and environment objects.  Maintaining this generic interface will still need to be retained for
			 * full flexibility of use. !!!!!
 */
double SpreadRateRothermelAlbini_Het(std::vector <double> SAV_ij,
                                     std::vector <double> w_o_ij,
                                     double fuelBedDepth,
                                     double M_x_1,
                                     std::vector <double> M_f_ij,
                                     double U, double slopeSteepness,
                                     std::vector <double> h_ij,
                                     std::vector <double> S_T_ij,
                                     std::vector <double> S_e_ij,
                                     std::vector <double> rho_p_ij,
                                     std::vector <int> liveDead,
                                     bool useWindLimit,
                                     UnitsType units,
                                     bool debug)
{
	SpreadCalcs calcs;

	calcs = SpreadCalcsRothermelAlbini_Het(SAV_ij, w_o_ij, fuelBedDepth, M_x_1, M_f_ij, U,
	                                       slopeSteepness, h_ij, S_T_ij, S_e_ij, rho_p_ij, liveDead,
	                                       useWindLimit, units, debug);

	return calcs.R;
}

//Allow the fuel particle properties to be set with single values or omitted.
double SpreadRateRothermelAlbini_Het(std::vector <double> SAV_ij,
                                     std::vector <double> w_o_ij,
                                     double fuelBedDepth,
                                     double M_x_1,
                                     std::vector <double> M_f_ij,
                                     double U, double slopeSteepness,
                                     double h,
                                     double S_T,
                                     double S_e,
                                     double rho_p,
                                     std::vector <int> liveDead,
                                     bool useWindLimit,
                                     UnitsType units,
                                     bool debug)
{
	int numFuelTypes = SAV_ij.size();
	std::vector <double> h_ij(numFuelTypes, h);
	std::vector <double> S_T_ij(numFuelTypes, S_T);
	std::vector <double> S_e_ij(numFuelTypes, S_e);
	std::vector <double> rho_p_ij(numFuelTypes, rho_p);
	double R;//Return value.

	R = SpreadRateRothermelAlbini_Het(SAV_ij, w_o_ij, fuelBedDepth,  M_x_1,
	                                  M_f_ij, U, slopeSteepness,
	                                  h_ij, S_T_ij, S_e_ij, rho_p_ij,
	                                  liveDead, useWindLimit, units, debug);

	return R;
}

/*This is a wrapper for SpreadRateRothermelAlbini_Het() that allows it to be called from R via .C():
We are only providing an interface for the the explicit arguments version for now.
The number of fuel types must be provided as the size of the input data is not carried from R.*/
extern "C" void SpreadRateRothermelAlbini_HetR(const double* SAV_ij, const double* w_o_ij,
                                               const double* fuelBedDepth,const double* M_x_1,
                                               const double*  M_f_ij, double* U,
                                               double* slopeSteepness, double* h_ij,
                                               double* S_T_ij, double* S_e_ij, double* rho_p_ij,
                                               int* liveDead, const int* useWindLimit,
                                               const int* units, const int* debug,
                                               const int* numFuelTypes, double* R)
{
	//Convert input arrays to vectors:
	std::vector<double> SAV_ijVec(SAV_ij, SAV_ij + *numFuelTypes);
	std::vector<double> w_o_ijVec(w_o_ij, w_o_ij + *numFuelTypes);

	std::vector<double> M_f_ijVec(M_f_ij, M_f_ij + *numFuelTypes);
	std::vector<double> h_ijVec(h_ij, h_ij + *numFuelTypes);
	std::vector<double> S_T_ijVec(S_T_ij, S_T_ij + *numFuelTypes);
	std::vector<double> S_e_ijVec(S_e_ij, S_e_ij + *numFuelTypes);

	std::vector<double> rho_p_ijVec(rho_p_ij, rho_p_ij + *numFuelTypes);
	std::vector<int> liveDeadVec(liveDead, liveDead + *numFuelTypes);

	bool useWindLimitBool;
	UnitsType cUnits;
	bool debugBool;

	//Code repeated from CalcWeightingsR() and SpreadRateRothermelAlbini_HomoR()!!!!!

	//The values used for live and dead differ from R to C++ (due to indexes).  Convert them:
	for (int k = 0; k < *numFuelTypes; k++)
	{
		liveDeadVec[k] -= 1;
	}

	//R logicals are passed as integer values.  This should be seamless to the user:
	if (*useWindLimit == 0)
	{
		useWindLimitBool = false;
	}
	else if (*useWindLimit == 1)
	{
		useWindLimitBool = true;
	}
	else//The NA value is possible it NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for useWindLimit.");
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

	if (*debug == 0)
	{
		debugBool = false;
	}
	else if (*debug == 1)
	{
		debugBool = true;
	}
	else//The NA value is possible it NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for debugBool.");
	}

	*R = SpreadRateRothermelAlbini_Het(SAV_ijVec, w_o_ijVec, *fuelBedDepth, *M_x_1, M_f_ijVec, *U,
	                                   *slopeSteepness, h_ijVec, S_T_ijVec, S_e_ijVec, rho_p_ijVec,
	                                   liveDeadVec, useWindLimitBool, cUnits, debugBool);
}

/* Albini 1976 modified Rothermel spread model for heterogeneous fuels, fuel model version:
 *
 * @par Input variables / parameters:
 * This version greatly reduces the number of parameters by taking a fuel model object that contains
 * the properties of the fuel bed.
 *
 * @param fuelModel A fuel model object.
 *
 * Environmental:
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 * @param M_f_ij Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
 *               This is not required if the fuel moisture has already been added to the fuel model.
 *
 * Settings:
 * @param useWindLimit Use the wind limit calculation or not.  Recent suggestion are that it not be used.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 *
 * @note Note since M_f_ij is optional the argument order is differnt than the equivalent version of
 * SpreadRateRothermelAlbini_Homo().
 */
double SpreadRateRothermelAlbini_Het(FuelModel fuelModel,
                                     double U, double slopeSteepness,
                                     std::vector <double> M_f_ij,
                                     bool useWindLimit,
                                     bool debug)
{
	SpreadCalcs calcs;

	calcs = SpreadCalcsRothermelAlbini_Het(fuelModel, U, slopeSteepness, M_f_ij,
	                                       useWindLimit, debug);

	return calcs.R;
}

//This variant returns many of the component calculated values in the spread rate calculation.
//The function parameters are the same as SpreadRateRothermelAlbini_Het().
//
//Returns: A SpreadCalcs object.
//Note: To avoid code repetition this function is called by the spread rate only functions.
SpreadCalcs SpreadCalcsRothermelAlbini_Het(std::vector <double> SAV_ij,
                                           std::vector <double> w_o_ij,
                                           double fuelBedDepth,
                                           double M_x_1,
                                           std::vector <double> M_f_ij,
                                           double U, double slopeSteepness,
                                           std::vector <double> h_ij,
                                           std::vector <double> S_T_ij,
                                           std::vector <double> S_e_ij,
                                           std::vector <double> rho_p_ij,
                                           std::vector <int> liveDead,
                                           bool useWindLimit,
                                           UnitsType units,
                                           bool debug)
{
	int numFuelTypes;
	FuelWeights weights;
	double meanPackingRatio, fuelBedSAV, optPackingRatio;
	double GammaPrime, IR, xi, phi_s, phi_w, rho_b_bar;//Scalar intermediates.
	double I_R, heatSource, heatSink;
	std::vector <double> w_n_i(2, 0), h_i(2, 0), M_x_i(2, 0), eta_M_i(2, 0), eta_s_i(2, 0), heatSink_i(2, 0);//Live/dead intermediates.
	std::vector <double> Q_ig_ij(w_o_ij.size(), 0);
	double R;//Spread rate.
	SpreadCalcs calcs;//Return value.

	ModelUnits = units;//Temporary!!!!! We need to decide if ModelUnits should be retained or not.

	//Parameter checking and processing:
	if (!SameLengths(SAV_ij, w_o_ij, M_f_ij, h_ij, S_T_ij, S_e_ij, rho_p_ij, liveDead))
	{
		//Stop("SpreadRateRothermelAlbini_Het() expects fuel array arguments to be of the same length.");
		Stop("SpreadXXXXRothermelAlbini_Het() expects fuel array arguments to be of the same length.");
		//Could pass in the function name?
	}

	numFuelTypes = SAV_ij.size();

	//Terms used in numerator and denominator:

	//Calculate the weights:
	weights = CalcWeightings(SAV_ij, w_o_ij, rho_p_ij, liveDead, units);

	//The heat source term (numerator) represents the heat flux from the fire front to the fuel in
	//front of it:
	//Numerator of Rothermel 1972 equation 75:
	//IR𝜉(1 + 𝜙w + 𝜙s)

	//Note: The bulk density is not used to calculate the packing ratio in the heterogeneous form:
	meanPackingRatio = MeanPackingRatio(w_o_ij, rho_p_ij, fuelBedDepth);//AKA beta_bar

	//For heterogeneous fuels we need to calculate the fuel bed level SAV:
	fuelBedSAV = FuelBedSAV(SAV_ij, weights.f_ij, weights.f_i, liveDead);

	optPackingRatio = OptimumPackingRatio(fuelBedSAV);

	//Reaction intensity:
	GammaPrime = OptimumReactionVelocity(meanPackingRatio, fuelBedSAV);
	w_n_i = NetFuelLoad_Het(w_o_ij, S_T_ij, weights.g_ij, liveDead);

	//Heat content by live/dead fuel category:
	h_i = LiveDeadHeatContent(h_ij, weights.f_ij, liveDead);

	//The live fuel moisture of extinction must be calculated:
	M_x_i[Dead] = M_x_1;
	M_x_i[Live] = LiveFuelMoistureOfExtinction(M_f_ij, M_x_1, w_o_ij, SAV_ij, liveDead);

	//Damping coefficients:
	eta_M_i = MoistureDampingCoefficient_Het(M_f_ij, M_x_i, weights.f_ij, liveDead);
	eta_s_i = MineralDampingCoefficient_Het(S_e_ij, weights.f_ij, liveDead);

	I_R = ReactionIntensity_Het(GammaPrime, w_n_i, h_i, eta_M_i, eta_s_i);

	//Other numerator terms:
	xi = PropagatingFluxRatio(meanPackingRatio, fuelBedSAV);
	phi_s = SlopeFactor(meanPackingRatio, slopeSteepness);

	//Apply wind limit check:
	if (useWindLimit)
	{
		U = WindLimit(U, I_R);
	}

	phi_w = WindFactor(fuelBedSAV, meanPackingRatio, optPackingRatio, U, units);

	heatSource = I_R * xi * (1 + phi_w + phi_s);

	//The heat sink term (denominator) represents the energy required to ignite the fuel in Btu/ft^3 |
	//kJ/m^3:
	//The heat sink term is calculated using weights without calculating epsilon explicitly.
	//Rothermel equation 77:
	//ρbεQig = ρb Σi fi Σj fij[exp(-138/σij)](Qig)ij

	rho_b_bar = MeanBulkDensity(w_o_ij, fuelBedDepth);
	Q_ig_ij = HeatOfPreignition(M_f_ij);

	//We'll do it in two steps:
	//Weight and size class:
	for (int k = 0; k < numFuelTypes; k++)
	{
		double savConst;

		if (units == US)
		{
			savConst= -138;
		}
		else
		{
			savConst = -4.527559;
		}

		heatSink_i[liveDead[k]] = heatSink_i[liveDead[k]] +
			weights.f_ij[k] * exp(savConst / SAV_ij[k]) * Q_ig_ij[k];
	}

	//Weigh and sum by live/dead category:
	heatSink = rho_b_bar * ((weights.f_i[Dead] * heatSink_i[Dead]) + (weights.f_i[Live] * heatSink_i[Live]));

	//Full spread calculation for heterogeneous fuels (same as homogeneous in this form):
	//Rothermel 1972 equation 75:
	//Rate of spread = heat source / heat sink
	//R = I_Rξ(1 + φw + φs) / ρbεQig
	//R = (I_R * xi * (1 + phi_w + phi_s)) / heatSink;
	R = heatSource / heatSink;

	//Pack the calculations into the struct:
	calcs.units = units;
	calcs.homogeneous = false;
	calcs.R = R;
	calcs.weights = weights;
	calcs.GammaPrime = GammaPrime;
	calcs.w_n_x = w_n_i;
	calcs.h_x = h_i;
	calcs.eta_M_x = eta_M_i;
	calcs.eta_s_x = eta_s_i;
	calcs.I_R = I_R;
	calcs.xi = xi;
	calcs.phi_s = phi_s;
	calcs.phi_w = phi_w;
	calcs.heatSource = I_R * xi * (1 + phi_s + phi_w);
	calcs.rho_b_x = rho_b_bar;
	calcs.epsilon = -1;//Not calculated for heterogeneous fuels.
	calcs.Q_ig_x = Q_ig_ij;
	calcs.heatSink = heatSink;
	calcs.cSAV = fuelBedSAV;
	calcs.packingRatio = meanPackingRatio;
	calcs.optimumPR = optPackingRatio;
	
	//For debugging:  Note: Could be replaced with a SpreadCalcs print function.
	if (debug)
	{
		Msg.Log("Heterogeneous Spread Calc components:");
		Msg.Log("Weights f_ij =", weights.f_ij);
		Msg.Log("Weights f_i =", weights.f_i);
		Msg.Log("Weights g_ij =", weights.g_ij);
		Msg.Log("GammaPrime =", GammaPrime);
		Msg.Log("w_n_i =", w_n_i);
		Msg.Log("h_i =", h_i);
		Msg.Log("eta_M_i =", eta_M_i);
		Msg.Log("eta_s_i =", eta_s_i);
		Msg.Log("I_R =", I_R);
		Msg.Log("Heat source =", heatSource);
		Msg.Log("rho_b_bar =", rho_b_bar);
		Msg.Log("Q_ig_ij =", Q_ig_ij);
		Msg.Log("Heat sink =", heatSink);
	}

	return calcs;
}

/*This is a wrapper for SpreadCalcsRothermelAlbini_Het() that allows it to be called from R via .C():
The number of fuel types must be provided as the size of the input data is not carried from R.

Note: This interface is incomplete, only returning some of the internal calculations.  It is meant
primarily for verification purposes.
The code is almost completely duplicated from SpreadRateRothermelAlbini_HetR().*/
extern "C" void SpreadCalcsRothermelAlbini_HetR(const double* SAV_ij, const double* w_o_ij,
                                                const double* fuelBedDepth,const double* M_x_1,
                                                const double*  M_f_ij, double* U,
                                                double* slopeSteepness, double* h_ij,
                                                double* S_T_ij, double* S_e_ij, double* rho_p_ij,
                                                int* liveDead, const int* useWindLimit,
                                                const int* units, const int* debug,
                                                const int* numFuelTypes, double* R, double* I_R,
                                                double* heatSource, double* heatSink)
{
	//Convert input arrays to vectors:
	std::vector<double> SAV_ijVec(SAV_ij, SAV_ij + *numFuelTypes);
	std::vector<double> w_o_ijVec(w_o_ij, w_o_ij + *numFuelTypes);

	std::vector<double> M_f_ijVec(M_f_ij, M_f_ij + *numFuelTypes);
	std::vector<double> h_ijVec(h_ij, h_ij + *numFuelTypes);
	std::vector<double> S_T_ijVec(S_T_ij, S_T_ij + *numFuelTypes);
	std::vector<double> S_e_ijVec(S_e_ij, S_e_ij + *numFuelTypes);

	std::vector<double> rho_p_ijVec(rho_p_ij, rho_p_ij + *numFuelTypes);
	std::vector<int> liveDeadVec(liveDead, liveDead + *numFuelTypes);

	bool useWindLimitBool;
	UnitsType cUnits;
	bool debugBool;

	SpreadCalcs calcs;

	//Code repeated from CalcWeightingsR() and SpreadRateRothermelAlbini_HomoR()!!!!!

	//The values used for live and dead differ from R to C++ (due to indexes).  Convert them:
	for (int k = 0; k < *numFuelTypes; k++)
	{
		liveDeadVec[k] -= 1;
	}

	//R logicals are passed as integer values.  This should be seamless to the user:
	if (*useWindLimit == 0)
	{
		useWindLimitBool = false;
	}
	else if (*useWindLimit == 1)
	{
		useWindLimitBool = true;
	}
	else//The NA value is possible it NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for useWindLimit.");
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

	if (*debug == 0)
	{
		debugBool = false;
	}
	else if (*debug == 1)
	{
		debugBool = true;
	}
	else//The NA value is possible if NAOK = TRUE in .C().
	{
		Stop("Invalid value passed for debugBool.");
	}

	calcs = SpreadCalcsRothermelAlbini_Het(SAV_ijVec, w_o_ijVec, *fuelBedDepth, *M_x_1, M_f_ijVec, *U,
	                                   *slopeSteepness, h_ijVec, S_T_ijVec, S_e_ijVec, rho_p_ijVec,
	                                   liveDeadVec, useWindLimitBool, cUnits, debugBool);
	
	//Return values:
	*R = calcs.R;
	*I_R = calcs.I_R;
	*heatSource = calcs.heatSource;
	*heatSink = calcs.heatSink;
}

/* Albini 1976 modified Rothermel spread model for heterogeneous fuels, fuel model version:
 *
 * @par Input variables / parameters:
 * This version greatly reduces the number of parameters by taking a fuel model object that contains
 * the properties of the fuel bed.
 *
 * @param fuelModel A fuel model object.
 *
 * Environmental:
 * @param U Wind speed at midflame height (ft/min | m/min).
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal distance).
 * @param M_f_ij Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).
 *               This is not required if the fuel moisture has already been added to the fuel model.
 *
 * Settings:
 * @param useWindLimit Use the wind limit calculation or not.  Recent suggestion are that it not be used.
 * @param debug Print calculation component values.  This may be removed in the future.
 *
 * @returns R rate of spread in ft/min | m/min.
 *
 * @note Note since M_f_ij is optional the argument order is differnt than the equivalent version of
 * SpreadRateRothermelAlbini_Homo().
 */
SpreadCalcs SpreadCalcsRothermelAlbini_Het(FuelModel fuelModel,
                                           double U, double slopeSteepness,
                                           std::vector <double> M_f_ij,
                                           bool useWindLimit,
                                           bool debug)
{
	std::vector <double> theM_f_ij, fmM_f_ij;
	SpreadCalcs calcs;

	//The fuel moisture may have already been added to the fuel model and if not should be supplied.
	//If both are provided that could be an error so warn and use the argument value as if intended to
	//override the fuel model:
	fmM_f_ij = fuelModel.GetM_f_ij();
	
	if (!M_f_ij.empty())
	{
		theM_f_ij = M_f_ij;

		if (!fmM_f_ij.empty())
		{
			Warning("M_f_ij argument provided but fuel model also has fuel moistures.");
		}
	}
	else if (!fmM_f_ij.empty())
	{
		theM_f_ij = fmM_f_ij;
	}
	else
	{
		Stop("Either the fuel model must have fuel moistures or M_f_ij argument must be provided.");
	}

	//Check the fuel model units.  Conversions could be applied:
	if (fuelModel.units == US)
	{
		if (fuelModel.w_o_Units == tonPerAc)
		{
			Stop("Loading units are in tons per acre and need to be in lb per ft^2.");
		}
		
		if (fuelModel.M_x_Units == Percent)
		{
			Stop("M_x units are in percent and need to be a fraction.");
		}
	}

	calcs = SpreadCalcsRothermelAlbini_Het(fuelModel.SAV_ij,
	                                       fuelModel.w_o_ij,
	                                       fuelModel.delta,//fuelBedDepth
	                                       fuelModel.M_x_1,
	                                       theM_f_ij,
	                                       U, slopeSteepness,
	                                       fuelModel.h_ij,
	                                       fuelModel.S_T_ij,
	                                       fuelModel.S_e_ij,
	                                       fuelModel.rho_p_ij,
	                                       fuelModel.liveDead,
	                                       useWindLimit,
	                                       fuelModel.units,
	                                       debug);

	return calcs;
}

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

/* Calculate the sum of a variable array of the form X_ij by the specified live/dead category:
 *
 * @param x_ij A numerical fuel model variable of the form X_ij.
 * @param liveDead An array indicating if each index in the input variables represents a dead or
 *                 live fuel category.
 * @param liveDeadCat The Live / Dead catagory value to get.
 *
 * This is a working draft.  It could return the sum (an array) for each class rather than specifying one.
 * This function is currently C++ only.
 */
//double SumByClass(std::vector<double> x_ij, std::vector<int> liveDead, FuelCategory liveDeadCat)
double SumByFuelCat(std::vector<double> x_ij, std::vector<int> liveDead, int liveDeadCat)
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

/*--------------------------------------------------------------------------------------------------
Related Fire Property Equations:
	These equations are not part of the Rothermel & Albini spread model per se but can be used with it
to calculate additional fire front properties.  See Andrews 2018 section 4 for more information.
--------------------------------------------------------------------------------------------------*/

//Effective Wind Speed:
//The effective wind speed combines the effect of wind and slope.  It was developed in:
// Albini, Frank A.
// Estimating wildfire behavior and effects.
// Gen. Tech. Rep. INT-GTR-30. Ogden, UT: U.S. Department of Agriculture, Forest Service,
// Intermountain Forest and Range Experiment Station. 92 p. 1976.
//for use in nomograhs but is used by some contexts to calculate the wind limit.
//
//Input variables / parameters:
//U = Wind speed at midflame height (ft/min | m/min).
//phi_w = The wind factor (dimensionless).
//phi_s = The slope factor (dimensionless).
//meanPackingRatio = Mean packing ratio (beta_bar) (dimensionless ratio).
//optPackingRatio = Optimum packing ratio (dimensionless).
//	 Note: optimum packing ratio is a function of SAV.
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//
//Output: effective wind speed at midflame height (ft/min | m/min)
//The functions called handle the units so no conversions need to be done here.
//
//Note: This is included for completeness.  It is not currently used by any other functions.
double EffectiveWindSpeed(double U, double phi_w, double phi_s, double meanPackingRatio,
                          double optPackingRatio, double cSAV, UnitsType units)
{
	double C, B, E;//Wind factor components.
	double phi_E;//Effective wind factor.
	double U_E;//Return value.

	C = WindFactorC(cSAV, units);
	B = WindFactorB(cSAV, units);
	E = WindFactorE(cSAV, units);

	//Effective wind factor:
	//Albini 1976(b) pg. 90:
	phi_E = phi_w * phi_s;

	//The effective wind factor calculated as in WindFactor():
	//ϕE = C * U_E^B (β/βop)^–E

	//Solve for the effective wind speed:
	//U_E = [ϕE (β/βop)E/C]^–B
	U_E = pow((phi_E * pow((meanPackingRatio / optPackingRatio), E) / C ), -B);

	return U_E;
}

//Residence Time:
//The residence time is the how long it takes the fire flame front to pass over a point on the
//ground.  This only considers the flame front and not the residual burning and smoldering. The
//calculation is from Anderson 1969.  This relationship is based on strong observational correlation.
//This can be used to help calculate energy transfer to soil.
//
//Input variables / parameters:
//cSAV = Characteristic surface-area-to-volume ratio (ft^2/ft^3 | cm^2/cm^3).
//
//Output units: minutes
double ResidenceTime(double cSAV, UnitsType units)//ResidenceTimeAnderson
{
	double t_r;//Returen value.

	if (units == US)
	{
		//The original equation predicts the residence time as 8 times the fuel diameter in inches.
		//We use the Rothermel relationship between diameter and SAV, d = 48/SAV:
		t_r = 384 / cSAV;
	}
	else
	{
		t_r = 12.59843 / cSAV;//384 / cmPerFt = 12.59843
		//Wilson 1980 & Andrews 2018 use 12.6 in the alternative formulation of Byram's fireline 
		//intensity.  See ByramsFirelineIntensity().
	}

	return t_r;
}

//Heat per Unit Area:
// This is the total energy released by the flame front as it passed per unit area.
//
//Input variables / parameters:
//I_R = Reaction intensity (Btu/ft^2/min | kJ/m^2/min).
//t_r = Residence time (min).
//
//Output units: Btu/ft^2 | kJ/m^2
//Inputs carry units.  No conversion necessary.
double HeatPerUnitArea(double I_R, double t_r)
{
	double H_A;//Return value.

	//Andrews 2018 section 4.3:
	H_A = I_R * t_r;
	//This can also be calculated as H_A = 384 * I_R/SAV.  See ResidenceTime().

	return H_A;
}

//Byram’s Fireline Intensity:
//	Byram 1959 defines fireline intensity as (equation 3.3):
//I_B = HwR
//Where:
//	H = heat content of the fuel (Btu/lb | kJ/kg)
//	w = weight of "available" fuel (lb/ft^2 | kg/m^2)
//	R = fire rate of spread (ft/s | m/s)
//Albini uses H_A as an approximation of H x W (Andrews 2018).  (Note: I can't find this in the text
//of Albini 1976.  It may be in the code.)
//
//Input variables / parameters:
//H_A = heat per unit area from the flame front (Btu/ft^2 | kJ/m^2)
//R = fire front rate of spread (ft/min | m/min)
//
//Output units: Btu/ft/s | kW/m
//The inputs carry the units.  No unit conversion is required.
//
//Andrews 2018 also presents the alternate formulation:
//I_B = (384/σ)I_R * R
//This is the same as ByramsFirelineIntensity(HeatPerUnitArea(I_R, ResidenceTime(SAV)), R).
double ByramsFirelineIntensity(double H_A, double R)
{
	double I_B;//Return value.

	I_B = H_A * R / 60;//Seconds / minute

	return I_B;
}

//Flame Length:
//	Calculate the flame length (not height) of the flame front.
//The equation is from Byram 1959 (equation 3.4) with notation from Brown and Davis 1973, page 175
//(per Andrews 2018, as I haven't been able to access the book yet.)
//
//Input variables / parameters:
//I_B = Byram's fireline intensity (Btu/ft/s | kW/m).
//
//Output units: ft | m
double ByramsFlameLength(double I_B, UnitsType units)
{
	double F_B;//Return value.

	if (units == US)
	{
		F_B = 0.45 * pow(I_B, 0.46);
	}
	else
	{
		F_B = 0.07749992 * pow(I_B, 0.46);
		//The errata of Wilson 1980 and Andrews 2018 give 0.0775 * I_B_m^0.46, with Wilson using L_f
		//instead of F_B.
	}

	return F_B;//L_f in Wilson 1980.
}

/*--------------------------------------------------------------------------------------------------
SpreadCalcs Functions:
	Should these be moved to new file
--------------------------------------------------------------------------------------------------*/

/** Print the contents of a SpreadCalcs object to an output stream.
 *
 * @param output The output stream to print to.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& SpreadCalcs::Print(std::ostream& output) const
{
	output << "Spread rate calculations:" << std::endl;

	if (units == US)
	{
		output << "Units: US" << std::endl;
	}
	else
	{
		output << "Units: Metric" << std::endl;
	}

	//output << "Homogeneous fuel model: " << homogeneous << std::endl;
	if (homogeneous)
	{
		output << "Homogeneous fuel model." << std::endl;
	}
	else
	{
		output << "Heterogeneous fuel model." << std::endl;
	}

	output << "R: " << R << std::endl;

	if (!homogeneous)//Weights: Only calculated for heterogeneous fuels.
	{
		output << "Weights:" << std::endl;
		//output << "f_ij: ";
		//PrintVector(output, weights.f_ij);
		output << "f_ij: " << VectorToStr(weights.f_ij) << std::endl;
		output << "f_i: " << VectorToStr(weights.f_i) << std::endl;
		output << "g_ij: " << VectorToStr(weights.g_ij) << std::endl;
	}

	//Heat source components:
	output << "GammaPrime: " << GammaPrime << std::endl;
	
	if (homogeneous)
	{
		output << "w_n: " << VectorToStr(w_n_x) << std::endl;
		output << "h: " << VectorToStr(h_x) << std::endl;
		output << "eta_M: " << VectorToStr(eta_M_x) << std::endl;
		output << "eta_s: " << VectorToStr(eta_s_x) << std::endl;
	}
	else
	{
		output << "w_n_i: " << VectorToStr(w_n_x) << std::endl;
		output << "h_i: " << VectorToStr(h_x) << std::endl;
		output << "eta_M_i: " << VectorToStr(eta_M_x) << std::endl;
		output << "eta_s_i: " << VectorToStr(eta_s_x) << std::endl;
	}

	output << "I_R: " << I_R << std::endl;
	output << "xi: " << xi << std::endl;
	output << "phi_s: " << phi_s << std::endl;
	output << "phi_w: " << phi_w << std::endl;
	output << "heatSource: " << heatSource << std::endl;

	//Heat sink components:
	if (homogeneous)
	{
		output << "rho_b: " << rho_b_x << std::endl;
		output << "epsilon: " << epsilon << std::endl;
	}
	else
	{
		output << "rho_b_bar: " << rho_b_x << std::endl;
		output << "epsilon: " << epsilon << " (Not used in heterogeneous fuel models)" << std::endl;
	}

	if (homogeneous)
	{
		output << "Q_ig: ";
	}
	else
	{
		output << "Q_ig_ij: ";
	}
	PrintVector(output, Q_ig_x);

	output << "heatSink: " << heatSink << std::endl;

	//Other components that can be informative:
	output << "cSAV: " << cSAV << std::endl;
	output << "packingRatio: " << packingRatio << std::endl;
	output << "optimumPR: " << optimumPR << std::endl;

	return output;
}

/* Overloaded stream print operator for SpreadCalcs.
 *
 */
std::ostream& operator<<(std::ostream& output, const SpreadCalcs& spreadCalcs)
{
	spreadCalcs.Print(output);
	return output;
}
