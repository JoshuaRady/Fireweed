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
