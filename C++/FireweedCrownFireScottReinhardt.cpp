/***************************************************************************************************
FireweedCrownFireScottReinhardt.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2026
Reference: Proj. 11 Exp. 26

Description:---------------------------------------------------------------------------------------
  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
crown fire equations of Scott & Reinhardt 2001.

References:-----------------------------------------------------------------------------------------

Van Wagner, C. E.
Conditions for the start and spread of crown fire.
Canadian Journal of Forest Research 7(1): 23-34, 1977. https://doi.org/10.1139/x77-004

Rothermel, R. C.
Predicting behavior and size of crown fires in the Northern Rocky Mountains.
Research Paper INT-RP-438. Ogden, UT: U.S. Department of Agriculture, Forest Service, Intermountain
Research Station. 46 pages, 1991. https://doi.org/10.2737/INT-RP-438

Scott, Joe H. and Reinhardt, Elizabeth D.
Assessing crown fire potential by linking models of surface and crown fire behavior.
Research Paper RMRS-RP-29. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky
Mountain Research Station. 59 pages, 2001.  https://doi.org/10.2737/RMRS-RP-29

Wind speed:-----------------------------------------------------------------------------------------
  The Scott & Reinhardt 2001 model is highly dependent on wind speed.  Wind speed controls the
transition from surface to crown fire.  The model takes open wind speeds (O) in km/hr as one of it's
main inputs a converts it to midflame wind speeds (U) in order to perform component calculations
from the Rothermel & Albini spread rate model.  The conversion factor 54.683 is used in several
places (WindConversionFactor below).  The value is a bit mysterious as it appears to be the
conversion from km/hr to the ft/min used in the original version of the Rothermel & Albini model
(1000 / 60 / mPerFt = 54.68066).  However, in my testing the calculations only give the expected
results when the Rothermel & Albini spread calculations are done in metric units.  Perhaps I am
misunderstanding this or there is a counteracting conversion in the other constants.  In any case
we retain the value as originally published but use kmPerHrToMPerMin for O to U wind speed
conversions elsewhere.

  Our functions all take O by default but since the open wind speed may not always be available we
have provided the ability to input U for the top level calculations.

***************************************************************************************************/

#include <utility>

#include "FireweedCrownFireScottReinhardt.h"
#include "FireweedMessaging.h"
#include "FireweedRAFireSpread.h"

//Constants:
const double WindConversionFactor = 54.683;
const double kmPerHrToMPerMin = 1000.0 / 60.0;//1000 m/km / 60 min/hr = m/min

//Code:---------------------------------------------------------------------------------------------

//Fuel Model Management:

/** Convert a fuel model to the physical properties of fuel model 10.
 * 
 * This is a helper function to deal with the fact that Rothermel 1991 performs calculations with
 * fuel model 10 regardless of the actual fuel model of the system being simulated.
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.
 * @param fuelModel10 Fuel model 10 with default values.  Only needed if fuelModel is not currently
 * fuel model 10.
 * 
 * @returns The converted fuel model.
 */
FuelModel ConvertToFuelModel10(const FuelModel& fuelModel, FuelModel fuelModel10)
{
	if (fuelModel.number == 10)
	{
		return fuelModel;
	}
	else
	{
		//Make sure the units are consistent:
		//We may want to force this to always be metric.
		if (fuelModel.units != fuelModel10.units)
		{
			fuelModel10.ConvertUnits(fuelModel.units);
		}
		
		if (fuelModel.cured || fuelModel.numClasses != 5)
		{
			Stop("Can't convert a fuel model with curing applied or more than the standard 5 classes.");
		}

		//Copy loadings:
		fuelModel10.w_o_ij = fuelModel.w_o_ij;
		
		std::vector <double> fmM_f_ij = fuelModel.GetM_f_ij();
		if (fmM_f_ij.empty())
		{
			Stop("M_f_ij must be provided in fuel model.");
		}
		//fuelModel10.m_f_ij = fmM_f_ij;
		fuelModel10.SetFuelMoisture(fmM_f_ij);

		return fuelModel10;
	}
}

/** Check that a fuel model meets the requirements for the Scott & Reinhardt 2001 crown fire
 * equations and modify it if needed.  This version doesn't enforce that it is fuel model 10.
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 *                  fuel model.
 * 
 * @returns Nothing, but the fuel model passed in may be updated.
 */
//FuelModel CheckFuelModel(const FuelModel& fuelModel, const bool convert = FALSE, FuelModel fuelModel10)
void CheckFuelModel(FuelModel& fuelModel)
{
	//Check for M_f_ij in the incoming fuel model.
	std::vector <double> fmM_f_ij = fuelModel.GetM_f_ij();
	if (fmM_f_ij.empty())
	{
		Stop("M_f_ij must be provided in fuel model.");
	}

	//Check units:
	if (fuelModel.units != Metric)
	{
		fuelModel.ConvertUnits(Metric);
	}
}

/** Check that a fuel model meets the requirements for the Scott & Reinhardt 2001 crown fire
 * equations and modify it if needed.  This version enforces that it is fuel model 10.
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.
 * @param fuelModel10 Fuel model 10 with default values.
 * 
 * @returns Nothing, but the fuel model passed in may be updated.
 */
void CheckFuelModel(FuelModel& fuelModel, FuelModel fuelModel10)
{
	CheckFuelModel(fuelModel);

	//Check model type:
	if (fuelModel.number != 10)
	{
		fuelModel = ConvertToFuelModel10(fuelModel, fuelModel10);
	}
}

//Wind Speed:---------------------------------------------------------------------------------------

/** Check the wind speed passed in for validity and convert to give both the open wind speed O and
 * the midflame wind speed U.
 * 
 * This is a utility to reduce code repetition in top level functions where we allow wind speed to
 * be passed in as either O or U.
 *
 * @param windSpeed The wind speed, by default O the open wind speed at 6.1 m (km/hr), otherwise,
 *                  U the wind speed at midflame height (m/min).  Type set by windType.
 * @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
 *            Needed even if U is supplied.
 * @param windType The wind speed type (see windSpeed).
 *
 * @returns The open wind speed at 6.1 m O and the midflame wind speed U as a pair.
 */
std::pair<double, double> CheckConvertWindSpeed(const double windSpeed, const double WRF,
                                                const char windType)
{
	double O = 0.0;
	double U = 0.0;

	//Validity checking:
	if (windSpeed < 0.0)
	{
		Stop("Invalid wind speed.");
	}
	//Add check for unexpectedly high wind speeds?

	//Check the wind inputs:
	if (windType == 'O')
	{
		O = windSpeed;
		U = O * WRF * kmPerHrToMPerMin;//If O is passed in calculate U.
	}
	else if (windType == 'U')
	{
		U = windSpeed;
		O = U / WRF;//If U is passed in we need to calculate O.
	}
	else
	{
		Stop("Must provide wind speed as O or U.");
	}

	std::pair<double, double> windSpeeds = {O, U};
	return windSpeeds;
}

//Crown Fire Spread Rate:---------------------------------------------------------------------------

/** Calculate an estimate of the active crown fire spread rate based the method of Rothermel 1991.
 * 
 * Rothermel 1991 calculates the crown fire spread rate as simple multiple of the surface fire
 * spread rate but the surface fire spread rate must be calculated with the physical properties of
 * fuel model 10 and a fixed 40% wind reduction factor.  This function handle the needed
 * conversions and returns the resulting rate.
 * 
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
 * those of fuel model 10.
 * @param O Open wind speed at 6.1 m (km/hr)
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
 * @param fuelModel10 Fuel model 10 with default values.  Only needed if fuelModel is not currently
 * fuel model 10.
 *
 * @returns The crown fire rate of spread (m/min).
 */
double SpreadRateCrownRothermel(FuelModel fuelModel, const double O, const double slopeSteepness,
                                FuelModel fuelModel10)
{
	CheckFuelModel(fuelModel, fuelModel10);
	
	double U = O * kmPerHrToMPerMin * 0.4;//Use fixed 40% WRF from Rothermel 1991.
	double R_surface = SpreadRateRothermelAlbini_Het(fuelModel, U, slopeSteepness);
	double R_active = 3.34 * R_surface;
	return R_active;
}

/** Compute the spread rate of a fire regardless of stage (surface, passive crown, or active crown).
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
 * those of fuel model 10.
 * @param windSpeed The wind speed, by default O the open wind speed at 6.1 m (km/hr), otherwise,
 *                  U the wind speed at midflame height (m/min).  Type set by windType.
 * @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
 *            Needed even if U is supplied.
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
 *                       distance).
 * @param CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
 *            fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
 * @param CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
 *            original Van Wagner notation.
 * @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
 * @param fuelModel10 Fuel model 10 with default values.  Only needed if fuelModel is not currently
 * fuel model 10.
 * @param windType The wind speed type (see windSpeed).
 *
 * @returns The fire spread rate (m/min).
 */
double SpreadRateCrownSR(FuelModel fuelModel, const double windSpeed, const double WRF,
                         const double slopeSteepness, const double CBD, const double CBH,
                         const double FMC, FuelModel fuelModel10, const char windType)
{
	CheckFuelModel(fuelModel);//Check the fuel model but don't convert the model number yet.

	std::pair <double, double> windSpeeds = CheckConvertWindSpeed(windSpeed, WRF, windType);
	double O = windSpeeds.first;
	double U = windSpeeds.second;

	//Calculate the component spread rates:
	//I think the surface spread rate should be calculated with the original fuel model but the text is
	//not explicit.  Confirm!
	double R_surface = SpreadRateRothermelAlbini_Het(fuelModel, U, slopeSteepness);
	double R_active = SpreadRateCrownRothermel(fuelModel, O, slopeSteepness, fuelModel10);

	//Get the crown fraction burned:
	double CFB = CrownFractionBurned(fuelModel, O, WRF, slopeSteepness, CBD, CBH, FMC, fuelModel10);

	//Compute the final rate using Scott & Reinhardt 2001 equation 21, pg. 19:
	double R_final = R_surface + CFB * (R_active - R_surface);
	return R_final;
}

//Crown Fire Initiation:----------------------------------------------------------------------------

/** Return the critical surface fire intensity needed to initiate crowning per Van Wagner 1977.
 *
 * @param CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
 *            original Van Wagner notation.
 * @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100).
 *
 * @returns The critical surface (fireline) intensity for crowning (kW/m, I_0 in Van Wagner notation,
 *          I'_initiation in Scott & Reinhardt notation).
 */
double CriticalCrowningIntensityVanWagner(const double CBH, const double FMC)
{
	//Scott & Reinhardt 2001 equation 11, pg. 13:
	//This is a combination of Van Wagner's equations 3, 4, and a value of C = 0.0010 from the text.
	double IPrime_initiation = pow(((CBH * (460.0 + 25.9 * FMC)) / 100.0), (3.0/2.0));
	return IPrime_initiation;
}

/** Return the critical surface fire rate of spread needed to initiate crowning using Van Wagner 1977.
 *
 * @param IPrime_initiation (kW/m, I_0 in Van Wagner notation, I'_initiation in Scott & Reinhardt
 *                           notation).
 * @param HPA Surface fire heat per area (kw/m^2).
 *
 * @returns The critical surface fire rate of spread for crowning (m/min).
 */
double CriticalCrowningROSVanWagner(const double IPrime_initiation, const double HPA)
{
	//Scott & Reinhardt 2001 equation 12, pg. 13:
	double R_initiation = (60.0 * IPrime_initiation) / HPA;
	return R_initiation;
}

/** Calculate the (minimum) critical crown fire rate of spread for active crowning.
 *
 * @param CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
 *            fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
 *
 * @returns The (minimum) critical crown fire rate of spread for active crowning (m/min).
 */
double CriticalActiveROSVanWagner(const double CBD)
{
	//Scott & Reinhardt 2001 equation 14, pg. 14:
	double Rprime_active = 3.0 / CBD;
	return Rprime_active;
}

//Crown Fire Development:----------------------------------------------------------------------------

/** Calculate the torching index, the open wind speed at which torching starts,
 *
 * @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
 * family of functions.
 * @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
 * @param CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
 *            original Van Wagner notation.
 * @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
 *
 * @returns TI, the torching index (open wind speed at 6.1 m, km/hr).
 */
double TorchingIndex(const SpreadCalcs spreadCalcs, const double WRF, const double CBH,
                     const double FMC)
{
	//Scott & Reinhardt 2001 equation 18, pg. 18:
	//TI = (1 / (54.683 * WRF)) * ((((60 * Iprime_initiation * rho_b * epsilon * Q_ig)/(HPA * xi * IR))
	//     - phi_s - 1) / (C * beta/beta_op)))^1/B
	//We break this equation up into the follow parts and compute them in steps:
	//TI = factor_L * (numerator_R / denominator_R)^1/B

	//We need the three wind factors:
	double C = WindFactorC(spreadCalcs.cSAV, Metric);
	double B = WindFactorB(spreadCalcs.cSAV, Metric);
	double E = WindFactorE(spreadCalcs.cSAV, Metric);

	//The leftmost factor:
	double factor_L = 1.0 / (WindConversionFactor * WRF);

	//The numerator on the right is Scott & Reinhardt 2001 equation 16 pg. 18:
	//phiPrime_w_intiation = ((60 * IPrime_initiation * rho_b * epsilon * Q_ig)/(HPA * xi * I_R))
	//                       - phi_s - 1
	//Where:
	//rho_b = bulk density
	//epsilon = the effective heating number
	//Q_ig = heat of preignition
	//HPA = heat per area
	//xi = the propagating flux ratio
	//I_R = Reaction intensity
	//phi_s = slope factor
	
	//The numerator of equation 16 contains the full heat sink (rho_b * epsilon * Q_ig).
	//We substitute for this since we do not calculate epsilon directly for heterogeneous fuels:
	//phiPrime_w_intiation = ((60 * IPrime_initiation * HeatSink)/(HPA * xi * I_R)) - phi_s - 1

	double IPrime_initiation = CriticalCrowningIntensityVanWagner(CBH, FMC);
	double HPA = HeatPerUnitArea(spreadCalcs.I_R, ResidenceTime(spreadCalcs.cSAV, Metric));
	double phiPrime_w_initiation = ((60.0 * IPrime_initiation * spreadCalcs.heatSink) /
	                                (HPA * spreadCalcs.xi * spreadCalcs.I_R)) - spreadCalcs.phi_s - 1;

	//The denominator on the right consists of most of the Rothermel & Albini spread model wind factor
	//calculation.  It has been rearranged to remove U:
	//C * (beta / beta_opt)^-E
	double denominator_R = C * pow((spreadCalcs.packingRatio / spreadCalcs.optimumPR), -E);

	//The complete calculation:
	double TI = factor_L * pow((phiPrime_w_initiation / denominator_R), (1/B));

	return TI;
}

/** Calculate the crowning index, the open wind speed at which active crowning occurs.
 *
 * @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
 * family of functions.
 * @param CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
 *            fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
 *
 * @returns CI, the crowning index (open wind speed at 6.1 m, km/hr).
 */
double CrowningIndex(const SpreadCalcs spreadCalcs, const double CBD)
{
	//Scott & Reinhardt 2001 equation 19, pg. 19 gives the midflame wind speed for active crowning.
	//UPrime_active = (1/54.683) * ((((3.0 / CBD) * rho_b * epsilon * Q_ig) / (3.34 * I_R * xi)
	//                               - phi_s - 1) / (C * (beta/beta_op)^-E))^1/B
	//To get the open wind speed for active crowning (O'_active = CI) using the assumptions of
	//Rothermel 1991 we can divide by a 40% wind reduction factor.
	//We break the resulting equation up into the follow parts and compute them in steps:
	//CI = factor_L * (numerator_R / denominator_R)^1/B

	//We need the three wind factors:
	double C = WindFactorC(spreadCalcs.cSAV, Metric);
	double B = WindFactorB(spreadCalcs.cSAV, Metric);
	double E = WindFactorE(spreadCalcs.cSAV, Metric);

	//The leftmost factor:
	double factor_L = 1.0 / WindConversionFactor / 0.4;

	//The numerator on the right is:
	//(((3.0 / CBD) * rho_b * epsilon * Q_ig) / (3.34 * I_R * xi)) - phi_s - 1
	//Which is the same as:
	//(((3.0 / CBD) * HeatSink) / 3.34 * I_R * xi) - phi_s - 1
	double numerator_R = (((3.0 / CBD) * spreadCalcs.heatSink) /
	                       (3.34 * spreadCalcs.I_R * spreadCalcs.xi)) - spreadCalcs.phi_s - 1.0;

	//The denominator on the right consists of most of the Rothermel & Albini spread model wind factor
	//calculation.  It has been rearranged to remove U  (It is the same as in the torching index.):
	//C * (beta / beta_opt)^-E
	double denominator_R = C * pow((spreadCalcs.packingRatio / spreadCalcs.optimumPR), -E);

	//The complete calculation:
	double CI = factor_L * pow((numerator_R / denominator_R), (1/B));

	//Note: Scott & Reinhardt 2001 uses equation 20, which is a simplified version of the above that
	//assumes parameters of fuel //model 10:
	//CI = 0.0457 * ((((3.0 / 3.34 / xi * HeatSink) / (I_R * CBD)) - phi_s - 1) / 0.001612)^0.7
	//We opt not to use the simplification to keep the code more parallel to the torching index
	//calculation and because rounding and small differences in conversion factors lead to slight
	//differences in the values we get with the explicit calculation.

	return CI;
}

/** Calculate the crown fraction burned.
 * 
 * The representation in Scott & Reinhardt 2001 does not purport to actually predict the crown
 * fraction burned for its own sake.  Instead it is used as a transition function between surface
 * and active crown fire spread rate calculations.  However, this has not stopped people from using
 * the output as an approximation of actual CFB.
 *
 * The underlying Rothermel 1991 crown fire equations expect fuel model 10, timber litter and
 * understory.  We accept any fuel model, which will be converted internally if needed.  This
 * currently requires  fuel model 10 to be passed in, which is inelegant.  If the path to the
 * library's input files could be resolved this argument could be eliminated.
 *
 * The original model uses open wind speed at 6.1 m.  We add the option to use the midflame wind
 * speed as is in the Rothermel & Albini spread model.  In this case the open wind speed will be
 * calculated based on the wind reduction factor.  Only O or U should be provided and note that the
 * units differ for the two.
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
 * those of fuel model 10.
 * @param windSpeed The wind speed, by default O the open wind speed at 6.1 m (km/hr), otherwise,
 *                  U the wind speed at midflame height (m/min).  Type set by windType.
 * @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
 *            Needed even in U is supplied.
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
 *                       distance).
 * @param CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
 *            fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
 * @param CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
 *            original Van Wagner notation.
 * @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100).
 * @param fuelModel10 Fuel model 10 with default values.  Only needed if fuelModel is not currently
 * fuel model 10.
 * @param windType The wind speed type (see windSpeed).
 *
 * @returns CFB, the crown fraction burned (fraction).
 */
double CrownFractionBurned(FuelModel fuelModel, const double windSpeed, const double WRF, 
                           const double slopeSteepness, const double CBD, const double CBH,
                           const double FMC, FuelModel fuelModel10, const char windType)
{
	double CFB;//Return value.
	CheckFuelModel(fuelModel, fuelModel10);

	std::pair <double, double> windSpeeds = CheckConvertWindSpeed(windSpeed, WRF, windType);
	double O = windSpeeds.first;
	double U = windSpeeds.second;
	//Note: We don't actually need an accurate value of U for computing CFB as the spread rate
	//sub-calculations used here are independent of the wind speed.  CI only depends on O.

	//Scott & Reinhardt 2001 never mentions the wind limit and we assume it is not used:
	SpreadCalcs spreadCalcs = SpreadCalcsRothermelAlbini_Het(fuelModel, U, slopeSteepness);

	double TI = TorchingIndex(spreadCalcs, WRF, CBH, FMC);
	double CI = CrowningIndex(spreadCalcs, CBD);

	if (O < TI)//Surface fire:
	{
		CFB = 0.0;
	}
	else if (O >= CI)//Fully active crown fire:
	{
		CFB = 1.0;
	}
	else//Passive crown fire, torching to full crown involvement:
	{
		//Calculate the transition as a linear function of wind speed:
		CFB = (1.0 / (CI - TI)) * (O - TI);
	}

	return CFB;
}

//Crown Fire Intensity:

/** Compute Byram’s fireline intensity of a fire regardless of stage (surface, passive crown, or
 * active crown).
 *
 * @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
 * fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
 * those of fuel model 10.
 * @param windSpeed The wind speed, by default O the open wind speed at 6.1 m (km/hr), otherwise,
 *                  U the wind speed at midflame height (m/min).  Type set by windType.
 * @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
 *            Needed even in U is supplied.
 * @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
 *                       distance).
 * @param CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
 *            fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
 * @param CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
 *            original Van Wagner notation.
 * @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
 * @param W_canopy Canopy fuel load (kg/m^2), AKA CFL.
 * @param fuelModel10 Fuel model 10 with default values.  Only needed if fuelModel is not currently
 *                    fuel model 10.
 * @param windType The wind speed type (see windSpeed).
 * @param H_canopy Heat yield of canopy fuel, heat content - heat of drying (kJ/kg, default from 
 *                 FARSITE).
 *
 * @returns Byram's fireline intensity (kW/m).
 * 
 * @note Requires verification!
 */
double CrownFireIntensity(FuelModel fuelModel, const double windSpeed, const double WRF,
                          const double slopeSteepness, const double CBD, const double CBH,
                          const double FMC, const double W_canopy, FuelModel fuelModel10,
                          const char windType, const double H_canopy)
{
	CheckFuelModel(fuelModel, fuelModel10);

	std::pair <double, double> windSpeeds = CheckConvertWindSpeed(windSpeed, WRF, windType);
	double O = windSpeeds.first;
	double U = windSpeeds.second;

	//Scott & Reinhardt 2001 equation 22, pg. 21:
	//I_final = ((HPA_surface + (W_canopy * H_canopy * CFB)) * R_final) / 60

	//Surface fire heat per area (kw/m^2):
	SpreadCalcs spreadCalcs = SpreadCalcsRothermelAlbini_Het(fuelModel, U, slopeSteepness);
	double HPA_surface = HeatPerUnitArea(spreadCalcs.I_R, ResidenceTime(spreadCalcs.cSAV, Metric));

	double CFB = CrownFractionBurned(fuelModel, O, WRF, slopeSteepness, CBD, CBH, FMC, fuelModel10);
	double R_final = SpreadRateCrownSR(fuelModel, O, WRF, slopeSteepness, CBD, CBH, FMC,
	                                   fuelModel10);

	double I_final = ((HPA_surface + (W_canopy * H_canopy * CFB)) * R_final) / 60.0;
	return I_final;
}
