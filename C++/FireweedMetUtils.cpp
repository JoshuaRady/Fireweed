/***************************************************************************************************
FireweedMetUtils.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:---------------------------------------------------------------------------------------
	This is part of the Fireweed wildfire code library.
	This file defines utility functions for calculating and converting meteorological properties.

	This is not intended to be a comprehensive library.  Rather it provides a few calculations that
are helpful when using and connecting other components of the library.

References:----------------------------------------------------------------------------------------

F. W. Murray.
On the Computation of Saturation Vapor Pressure.
Journal of Applied Meteorology and Climatology 6(1): 203-204, 1967. DOI: https://doi.org/10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2
	This paper compares the calculations of Goff and Gratch (1946) and Tetens (1930).  I have not
been able to find the original text of Tetens 1930, which is in German.

Arden L. Buck.
New Equations for Computing Vapor Pressure and Enhancement Factor.
Journal of Applied Meteorology and Climatology 20(12): 527-1532, 1981. https://doi-org.library.proxy.mbl.edu/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2
	This is the original source of the Buck equation.

Arden L. Buck.
Model CR-1A hygrometer with autofill operating manual.
Buck Research Instruments LLC: Aurora, CO, USA. 2012 (1996).
Obtained from https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf.
	This manual contains an appendix in which Buck provides a set of equations converting different
humidity related quantities, including new parameters for the Buck equations and a enhancement
factor equation, which can replace the table values from the original paper.
	A manual is an unusual reference but it cited by several sources.  The manual was revised in
2012 but this appendix has was last revised in 1996 and sources cite it as Buck 1996.

***************************************************************************************************/

#include <cmath>//pow()

#include "FireweedMessaging.h"
#include "FireweedMetUtils.h"

//Saturation Vapor Pressure:------------------------------------------------------------------------

/** Calculate the saturation vapor pressure of water vapor using the Tetens equation:
 * 
 * The Tetens equation is the Magnus equation form with different parameters.  I have not be able to
 * obtain the original German language manuscript (Otto Tetans 1930).  This function is based on that
 * presented in Murray 1967.
 * 
 * Note that this calculation is for pure water vapor.  The difference from moist air is small but not
 * negligible.
 * 
 * @param tempC The air temperature (degrees Celsius).
 * 
 * @returns: The saturation vapor pressure of water vapor, P_s, P_sat, vp_sat, or e_s (hPa = millibar).
 */
double SaturationVaporPressureTetens(double tempC)
{
	double e_s;//Return value.

	//The formulation in Murray 1967 is:
	//e_s = 6.1078 * exp(a x (T - 273.16) / (T - b))
	//Where T is temperature in Kelvin and a and b are parameters that vary for e_s over ice and water.
	//The equations below are rearranged to take temperature in degrees Celsius.
	if (tempC > 0)//Over water:
	{
		e_s = 6.1078 * exp((17.2693882 * tempC) / (tempC + 237.3));
	}
	else//Over ice:
	{
		e_s = 6.1078 * exp((21.8745584 * tempC) / (tempC + 265.5));
	}

	return e_s;
}

/** Calculate the saturation vapor pressure of water in moist air using the Arden Buck equation:
 * 
 * The Arden Buck equation (Buck 1981) gives more accurate estimates of saturation vapor pressure than
 * many earlier equations and includes an "enhancement factor" that adjusts the output to that of
 * moist air.  This function incorporates improvements from Buck 1996.
 * 
 * @param tempC The air temperature (degrees Celsius).
 * @param p_hPa (Atmospheric) pressure (hPa).  Defaults to typical air pressure at sea level.
 *              This is used for the enhancement factor calculation.
 * 
 * @returns The saturation vapor pressure of water in moist air, P_s or e' sub s (hPa = millibars).
 * 
 * p_hPa could be made optional.  If p_hPa was null the enhancement factor would not be applied and
 * saturation vapor pressure of water vapor would be returned.
 */
double SaturationVaporPressureBuck(double tempC, double p_hPa)
{
	double P_s;//Return value.

	//Check temp:
	//The upper bound of the range is unclear with these parameters, it may be as high as 100C.
	if (tempC < -80 || tempC > 50)
	{
		Warning("Temperature" + std::to_string(tempC) +
				"C may be outside the accuracy range of this Buck equation parameter set.");
	}

	//Buck 1981 actually includes a set related equations with three or four parameters fit for
	//performance optimized for several different temperature ranges.  This used the four parameter
	//version with values as follows.  These values are shared with some other conversions equations
	//in Buck 1996 and could be made global constants.

	//Buck coefficients for over water:
	//These are the values from Buck 1996.  They are most similar to Buck 1981 parameter set e_w4.
	double buck_a_w = 6.1121;//hPa
	double buck_b_w = 18.678;
	double buck_c_w = 257.14;//Degrees C
	double buck_d_w = 234.5;//Degrees C
	
	//Buck coefficients for over ice:
	//These are the same for Buck 1981 parameter set e_i3 and Buck 1996.
	double buck_a_i = 6.1115;//hPa
	double buck_b_i = 23.036;
	double buck_c_i = 279.82;//Degrees C
	double buck_d_i = 333.7;//Degrees C
	
	//The original Buck 1981 notation is (equation 4a):
	//e_w or e_i = a x exp[(b - T/d)T / (T + c)]
	//Where e_w is e_s over water and e_i is e_s over ice.
	//This is rearranged below to:
	//e_w or e_i = a x exp((b - T/d) x (T / (T + c)))
	//This is multiplied by the enhancement factor, ƒ_x in Buck 1981, EF_x in Buck 1996, to produce the
	//moist air value.  This is denoted e' sub w and e' sub i respectively.

	if (tempC > 0)//Over water:
	{
		//Enhancement factor over water:
		//This EF version as function of pressure and temp is from Buck 1996.
		double EF_w = 1 + pow(10, -4) * (7.2 + p_hPa * (0.0320 + 5.9 * pow(10, -6) * pow(tempC, 2)));

		P_s = EF_w * buck_a_w * exp((buck_b_w - (tempC / buck_d_w)) * (tempC / (buck_c_w + tempC)));
	}
	else//Over ice:
	{
		//Enhancement factor over ice:
		//This EF version as a function of pressure and temp is from Buck 1996.
		double EF_i = 1 + pow(10, -4) * (2.2 + p_hPa * (0.0383 + 6.4 * pow(10, -6) * pow(tempC, 2)));
	
		P_s = EF_i * buck_a_i * exp((buck_b_i - (tempC / buck_d_i)) * (tempC / (buck_c_i + tempC)));
	}

	return P_s;
}

//Relative Humidity:--------------------------------------------------------------------------------

/* Calculate relative humidity from vapor pressures.
 * 
 * This function is not that useful by itself because having both vapor pressures as inputs is not
 * that common.  Rather this function is used by other functions with more useful inputs.  This also
 * allows error checks to be centralized.
 * 
 * @aram P The partial pressure of water vapor in air (hPa or other).
 * @aram P_s The saturation vapor pressure of water in air (hPa or other).
 * 
 * @returns Relative humidity (%).
 * 
 * Note: As long as the units of both pressures are the same the the result will be valid.
 */
double RHfromVP(double P, double P_s)
{
	double rhPct;//Return value.

	//Negative pressures are not physically possible:
	if (P < 0)
	{
		Stop("The partial pressure of water can't be negative.");
	}
	if (P_s < 0)
	{
		Stop("The saturation vapor pressure of water can't be negative.");
	}

	rhPct = P / P_s * 100;

	//Check that the calulated value is valid:
	if (rhPct > 100)
	{
		//Because of small numerical inaccuracies in the pressures the calculated RH might exceed 100%.
		//The error should be small though, and the user may want to know:
		//I need to determine a reasonable upper cutoff.
		Warning("Calculated relative humidity it above 100: " + std::to_string(rhPct) + "%. Setting to 100.");
		rhPct = 100;
	}

	return rhPct;
}

/* Convert dew point at a given temperature to relative humidity using the Buck equation.
 * 
 * @aram tempC The air temperature (degrees Celsius).
 * @aram T_d Dew point temperature (degrees Celsius).
 * @aram p_hPa (Atmospheric) pressure (hPa).  Defaults to typical air pressure at sea level.
 * 
 * @returns Relative humidity (%).
 */
double RHfromDewPointBuck(double tempC, double T_d, double p_hPa)
{
	double rhPct;//Return value.

	double P_s = SaturationVaporPressureBuck(tempC, p_hPa);
	double P = SaturationVaporPressureBuck(T_d, p_hPa);

	rhPct = RHfromVP(P, P_s);
	return rhPct;
}

//Vapor Pressure Deficit:---------------------------------------------------------------------------
//Vapor pressure deficit seems to be most commonly reported in kPs.  We used hPa here to be
//consistent across equations.

/** Calculate the vapor pressure deficit (VPD) from relative humidity, temperature, and the
 * saturation vapor pressure of water.
 * 
 * This version allows the user to choose how they calculate P_s.
 * 
 * @aram tempC The air temperature (degrees Celsius).
 * @aram rhPct Relative humidity (%).
 * @aram P_s The saturation vapor pressure of water in moist air (hPa).
 * 
 * @returns Vapor pressure deficit (hPa).
 */
double VPDfromRH(double tempC, double rhPct, double P_s)
{
	double vpd_hPa = P_s * (1 - (rhPct / 100));

	//Check value for validity:
	if (vpd_hPa < 0)
	{
		Stop("Invalid VPD:" + std::to_string(vpd_hPa));
	}

	return vpd_hPa;
}

/** Calculate the vapor pressure deficit (VPD) from relative humidity and temperature and using the
 * Buck equation for the saturation vapor pressure of water in moist air:
 * 
 * @aram tempC The air temperature (degrees Celsius).
 * @aram rhPct Relactive humidity (%).
 * @aram p_hPa (Atmospheric) pressure (hPa).  Defaults to typical air pressure at sea level.
 * 
 * @returns Vapor pressure deficit (hPa).
 */
double VPDfromRHBuck(double tempC, double rhPct, double p_hPa)
{
	double P_s = SaturationVaporPressureBuck(tempC, p_hPa);
	double vpd_hPa = VPDfromRH(tempC, rhPct, P_s);
	return vpd_hPa;
}
