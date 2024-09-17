#FireweedMetUtils.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 20
#
#Description:---------------------------------------------------------------------------------------
#  This is part of the Fireweed wildfire code library.
#  This file contains utility function for meteorological properties.
#
#References:----------------------------------------------------------------------------------------
#
#F. W. Murray.
#On the Computation of Saturation Vapor Pressure.
#Journal of Applied Meteorology and Climatology 6(1): 203-204, 1967. DOI: https://doi.org/10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2
#  This paper compares the calculations of Goff and Gratch (1946) and Tetens (1930).  I have not
#been able to find the original text of Tetens 1930, which is in German.
#
#Arden L. Buck.
#New Equations for Computing Vapor Pressure and Enhancement Factor.
#Journal of Applied Meteorology and Climatology 20(12): 527-1532, 1981. https://doi-org.library.proxy.mbl.edu/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2
#  This is the original source of the Buck equation.
#
#Arden L. Buck.
#Model CR-1A hygrometer with autofill operating manual.
#Buck Research Instruments LLC: Aurora, CO, USA. 2012 (1996).
#Obtained from https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf.
#  This manual contains an appendix in which Buck provides a set of equations converting different
#humidity related quantities, including new parameters for the Buck equations and a enhancement
#factor equation, which can replace the table values from the original paper.
#  A manual is an unusual reference but it cited by several sources.  The manual was revised in
#2012 but this appendix has was last revised in 1996 and sources cite it as Buck 1996.
#
#___________________________________________________________________________________________________

# source("FireweedUnits.r")
# source("FireweedUtils.r")

#Code:----------------------------------------------------------------------------------------------

#Calculate the saturation vapor pressure of water vapor using the Tetens equation:
#
#The Tetens equation is the Magnus equation form with different parameters.  I have not be able to
#obtain the original German language manuscript (Otto Tetans 1930).  This function is based on that
#presented in Murray 1967.
#
#Note that this calculation is for pure water vapor.  The difference from moist air is small but not
#negligible.
#
#Parameters:
#tempC = The air temperature (degrees Celsius).
#
#Returns: The saturation vapor pressure of water vapor, P_s or e_s (hPa = millibar).
SaturationVaporPressureTetens <- function(tempC)
{
  #The formulation in Murray 1967 is:
  #e_s = 6.1078 * exp(a x (T - 273.16) / (T - b))
  #Where T is temperature in Kelvin and a and b are parameters that vary for e_s over ice and water.
  #The equations below are rearranged to take temperature in degrees Celsius.
  if (tempC > 0)#Over water:
  {
    e_s = 6.1078 * exp((17.2693882 * tempC) / (tempC + 237.3))
  }
  else#Over ice:
  {
    e_s = 6.1078 * exp((21.8745584 * tempC) / (tempC + 265.5))
  }
  
  return(e_s)
}

#Calculate the saturation vapor pressure of water in moist air using the Arden Buck equation:
#
#The Arden Buck equation (Buck 1981) gives more accurate estimates of saturation vapor pressure than
#many earlier equations and includes an "enhancement factor" that adjusts the output to that of
#moist air.  This function incorporates improvements from Buck 1996.
#
#Parameters:
#tempC = The air temperature (degrees Celsius).
#p_hPa = (Atmospheric) pressure (hPa).  Defaults to typical air pressure at sea level.
#  This is used for the enhangment factor calculation,
#
#Returns: The saturation vapor pressure of water in moist air, P_s or e' sub s (hPa = millibars).
#
#p_hPa could be made optional.  If p_hPa was null the enhancement factor would not be applied and
#saturation vapor pressure of water vapor would be returned.
SaturationVaporPressureBuck <- function(tempC, p_hPa = 1013)
{
  #Check temp:
  #The upper bound of the range is unclear with these parameters, it may be as high as 100C.
  if (tempC < -80 || tempC > 50)
  {
    warning(paste("Temperature", tempC, "C is outside accuracy range of Buck equation parameters."))
  }
  
  #Buck 1981 actually includes a set related equations with three or four parameters fit for
  #performance optimized for several different temperature ranges.  This used the four parameter
  #version with values as follows.  These values are shared with some other conversions equations
  #in Buck 1996 and could be made global constants.
  
  #Buck coefficients for over water:
  #These are the values from Buck 1996.  They are most similar to Buck 1981 parameter set e_w4.
  buck_a_w = 6.1121#hPa
  buck_b_w = 18.678
  buck_c_w = 257.14#Degrees C
  buck_d_w = 234.5#Degrees C
  
  #Buck coefficients for over ice:
  #These are the same for Buck 1981 parameter set e_i3 and Buck 1996.
  buck_a_i = 6.1115#hPa
  buck_b_i = 23.036
  buck_c_i = 279.82#Degrees C
  buck_d_i = 333.7#Degrees C
  
  #The original Buck 1981 notation is (equation 4a):
  #e_w or e_i = a x exp[(b - T/d)T / (T + c)]
  #Where e_w is e_s over water and e_i is e_s over ice.
  #This is rearranged below to:
  #e_w or e_i = a x exp((b - T/d) x (T / (T + c)))
  #This is multiplied by the enhancement factor, ƒ_x in Buck 1981, EF_x in Buck 1996, to produce the
  #moist air value.  This is denoted e' sub w and e' sub i respectively.
  
  if (tempC > 0)#Over water:
  {
    #Enhancement factor over water:
    #This EF version as function of pressure and temp is from Buck 1996.
    EF_w = 1 + 10^-4 * (7.2 + p_hPa * (0.0320 + 5.9 * 10^-6 * tempC^2))
    
    P_s = EF_w * buck_a_w * exp((buck_b_w - (tempC / buck_d_w)) * (tempC / (buck_c_w + tempC)))
  }
  else#Over ice:
  {
    #Enhancement factor over ice:
    #This EF version as a function of pressure and temp is from Buck 1996.
    EF_i = 1 + 10^-4 * (2.2 + p_hPa * (0.0383 + 6.4 * 10^-6 * tempC^2))
    
    P_s = EF_i * buck_a_i * exp((buck_b_i - (tempC / buck_d_i)) * (tempC / (buck_c_i + tempC)))
  }
  
  return(P_s)
}
