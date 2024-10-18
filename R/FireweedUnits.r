#FireweedUnits.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 19
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed wildfire code library.
#  This file contains unit constants, conversion factors, and functions used in the library.
#
#___________________________________________________________________________________________________

source("FireweedUtils.r")

#The units class to be used or in use is specified by a string.  The options are 'US' for United
#States customary units or 'Metric' for metric.

#Symbols for identifying specific units are: lbPer_ft2, tonPerAc, kgPer_m2, Percent, Fraction.

#Unit Conversion Factors:---------------------------------------------------------------------------

#Length: (exact per international yard and pound act)
cmPerIn = 2.54
cmPerFt = 30.48
mPerFt = 0.3048
ftPerM = 3.28084#1 / mPerFt
ftPerMi = 5280#* for conversion of windspeed (U, MPH * ftPerMi / 60 = ft/min)

#SAV is in ft^2/ft^3 = 1/ft or cm^2/cm^3 = 1/cm
#Therefore units convert: ft^2/ft^3 * cmPerFt^2/cmPerFt^2 = 1/ft * 1/cmPerFt = 1/cm
#So: SAVft * 1/cmPerFt = SAVft / cmPerFt = SAVcm

#Area:
ft2PerAcre = 43560

#Mass:
kgPerLb = 0.453592
lbsPerTon = 2000

#Density:
lbPerFtCuToKgPerMCu = 16.0185#kgPerLb * (ftPerM)^3, 16.01846337396

#JPerBtu = 1055.06 or 1,054.35
#The definition of a BTU can vary resulting in several different conversion factors.  Wilson 1980
#seems to have used a value close to the themochemical value of 1.05435 J/BTU, based on his heat of
#preignition conversion.  We will use that to be consistent with his converted constant values.
#The IT value of 1.05506 would be a reasonable alternative.
kJPerBtu = 1.05435

#tons/ac -> lb/ft^2: (See fuel loading note in FireweedFuelModels.r.)
tonsPerAcToLbPerSqFt = lbsPerTon / ft2PerAcre#*

#Code:----------------------------------------------------------------------------------------------

#Convert a temperature in degrees Celsius to degrees Fahrenheit:
#
#Parameters:
#degreesC = The temperature to convert (degrees Celsius).
#
#Returns: Temperature in degrees Fahrenheit.
CtoF <- function(degreesC)
{
  return((degreesC * 9/5) + 32)
}

#Convert a temperature in degrees Fahrenheit to degrees Celsius:
#
#Parameters:
#degreesF = The temperature to convert (degrees Fahrenheit).
#
#Returns: Temperature in degrees Celsius.
FtoC <- function(degreesF)
{
  return((degreesF - 32) * 5/9)
}

#Convert from percent slope to slope steepness as used by the Rothermel Albini spread model.
#
#Parameters:
#slopePct = The slope or grade etc. of the landscape in percent.
#
#Returns: Slope steepness (unitless fraction: vertical rise / horizontal distance), AKA tan ϕ.
#Slope values are forced to postive values.
SlopePctToSteepness <- function(slopePct)
{
  #The percent slope is simply the rise / run x 100%:
  slopeSteepness = abs(slopePct) / 100
  
  CheckSlope(slopeSteepness)#Check that the value is reasonable.
  
  return(slopeSteepness)
}

#Convert from slope in degrees to slope steepness as used by the Rothermel Albini spread model.
#
#Parameters:
#slopeDegrees = The slope or grade etc. of the landscape in degrees.
#
#Returns: Slope steepness (unitless fraction: vertical rise / horizontal distance), AKA tan ϕ.
#Slope values are forced to postive values.
SlopeDegreesToSteepness <- function(slopeDegrees)
{
  #While any degree value may be mathematically valid values beyond 90 are not really used to
  #express slopes.  Large values may indicate the input value was miscomputed.  Negative slopes
  #could be valid too but are questionable.  Allow but warn:
  if (!InRange(slopeDegrees, 0, 90))
  {
    warning("Slope is outside the expected range of 0 - 90 degrees.")
  }
  
  #tan() takes degrees:
  slopeSteepness = tan(abs(slopeDegrees) * (pi/180))
  
  CheckSlope(slopeSteepness)#Check that the value is reasonable.
  
  return(slopeSteepness)
}

#Check if a slope is reasonable:
#
#The problem with fractional or percent slope is that is becomes huge as it approaches 90 degrees.
#However, realistically very high slopes should be very rare on the landscape.
#This check is minimal and might be improved.
#
#Returns: Nothing.  Posts a warning.  Could return a status instead.  ValidSlope()?
CheckSlope <- function(slopeSteepness)
{
  if (!InRange(slopePct, 0, 11.43005))#0 - ~85 degrees.
  {
    warning("Questionable slope value.")
  }
}
