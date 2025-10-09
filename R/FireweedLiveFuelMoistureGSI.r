#FireweedLiveFuelMoistureGSI.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 20
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
#GSI live fuel moisture model used in the National Danger Rating System (NFDRS) 2016.
#
#  Growing Season Index (GSI) is a index that can be used to predict phenological stage, such as
#greenup and greendown globally.  It is calculated from daily weather conditions and during the
#growing season its value can be viewed as a predictor of moisture content in plants that respond
#quickly to environmental conditions.  NFDRS uses GSI as a basis for prediction of live fuel
#moisture.
#  GSI, as well as herbaceous, and woody live fuel moisture are calculated daily but these values
#should not be used directly.  GSI values can change quickly with daily weather and plants should
#not be expected to respond as quickly.  Therefore the standard procedure is to a 21 day running
#average.  The authors also suggest that using monthly mean values as inputs may also provide an
#appropriate smoothing.
#
#References:----------------------------------------------------------------------------------------
#
#A Generalized, Bioclimatic Index to Predict Foliar Phenology in Response to Climate.
#William M. Jolly, Ramakrishna Nemani, and Steven W. Running.
#Global Change Biology 11(4): 619 - 632, 2005.  DOI:10.1111/j.1365-2486.2005.00930.x
#  This paper presents the Growing Season Index (GSI), an index for predicting phenological status
#driven by minimum daily temperature, VPD, and photoperiod / day length.
#
#Overview of NFDRS2016.
#W. Matt Jolly, USFS, RMRS, Missoula Fire Sceinces Laboratory.
#National NFDRS 2016 Rollout Workshop, 4/28/2018.
#https://gacc.nifc.gov/eacc/predictive_services/fuels_fire-danger/documents/Overview%20of%20NFDRS2016%20and%20Implementation%20and%20Evaluation.pdf
#  This presentation contains the equations used to transform GSI values to live fuel moisture
#predictions at is is done in the National Danger Rating System 2016.  This is not an ideal
#reference but it is the only documentation I have been able to find thus far.
#
#___________________________________________________________________________________________________

#Code:----------------------------------------------------------------------------------------------

#Calculate the Growing Season Index for a location based on daily weather conditions:
#
#Parameters:
#tempCMin = Minimum daily temperature (degrees Celsius).
#vpdPa = Mean daily vapor pressure (VPD) deficit (Pa).  The documentation mentions maximum daily VPD
#  (VPD_max, not to be confused with the VPD_max parameter in the calculation) can also be used but
#  it is unclear if there are any changes in the calculation in this case.
#dayLength = Daylight length or photoperiod (seconds).
#
#Returns: The Growing Season Index (GSI) ranging from 0 (inactive) - 1 (unconstrained) (unitless).
GrowingSeasonIndex<- function(tempCMin, vpdPa, dayLength)
{
  #Model parameters defining the range over with plants go from phenologially inactive to
  #unconstrained. 
  #These could be made into function default parameters.  While these are the standard values used
  #in all published work they could in theory be used to turn the model for specific conditions.
  T_MMin = -2#C
  T_MMax = 5#C
  VPD_Min = 900#Pa
  VPD_Max = 4100#Pa
  Photo_Min = 36000#Seconds = 10 hours
  Photo_Max = 39600#Seconds = 11 hours
  
  #Validity checking:
  #For temperatures any value above 0K is valid, but this doesn't seem like a very useful check.
  #We could check that the value is in a likely range.
  #VPD reaches ~10,000 Pa at 46C and 0% RH at sea level.  That is extremely high so this may not
  #catch much, especially since we use Pa.  Using the wrong units will likely decrease the value
  #(e.g. using kPa or hPa):
  if (vpdPa < 0)
  {
    stop(paste("VPD must be positive. vpdPa =", vpdPa))
  }
  else if (vpdPa > 10000)
  {
    warning(paste("Very unlikely VPD value:", vpdPa))
  }
  #Day length can reach extreems in the far north.  Since we are using seconds we are unlikely to
  #catch the wrong units (e.g. hours or minutes):
  if (dayLength < 0 || dayLength > 86400)#24 hours
  {
    stop(paste("Invalid day length:", dayLength))
  }
  
  #Minimum temperature index:
  #Jolly, Nemani, & Running 2005 equation 1.
  #tempCMin = T_min in the original notation.
  if (tempCMin <= T_MMin)
  {
    iT_Min = 0
  }
  else if (tempCMin < T_MMax)#(tempCMin > T_MMin) is implied.
  {
    iT_Min = (tempCMin - T_MMin) / (T_MMax - T_MMin)
  }
  else#(tempCMin >= T_MMax) is implied.
  {
    iT_Min = 1
  }
  
  #VPD index:
  #Jolly, Nemani, & Running 2005 equation 2.
  #vpdPa = VPD in the original notation.
  if (vpdPa >= VPD_Max)
  {
    iVPD = 0
  }
  else if (vpdPa > VPD_Min)#(VPD_Max > vpdPa) is implied.
  {
    iVPD = 1 - ((vpdPa - VPD_Min) / (VPD_Max - VPD_Min))
  }
  else#(vpdPa <= VPD_Min) is implied.
  {
    iVPD = 1
  }
  
  #Photoperiod index:
  #Jolly, Nemani, & Running 2005 equation 3.
  #dayLength = Photo in the original notation.
  if (dayLength <= Photo_Min)
  {
    iPhoto = 0
  }
  else if (dayLength < Photo_Max)#(dayLength > Photo_Min) is implied.
  {
    iPhoto = (dayLength - Photo_Min) / (Photo_Max - Photo_Min)
  }
  else#(dayLength >= Photo_Max) is implied.
  {
    iPhoto = 1
  }
  
  gsi = iT_Min * iVPD * iPhoto
  return(gsi)
}

#Convert the Growing Season Index to live fuel moisture:
#GSI is linearly scaled to a range of live fuel moisture estimates for a live fuel type with the
#moisture range passed in.
#This generic function should not normally be called directly.  Use HerbaceousLiveFuelMoisture()
#and WoodyLiveFuelMoisture() for standard NFDRS 2016 behavior.
#
#Parameters:
#gsi = The Growing Season Index value for the location.
#lfmMin = Minimum live fuel moisture parameter.
#lfmMax = Maximum live fuel moisture parameter.
#gu = GSI greenup threshold, defaults to 0.5 (unitless index value).
#
#Returns: Percent live fuel moisture (water weight/dry fuel weight * 100%).
GSI_LiveFuelMoisture <- function(gsi, lfmMin, lfmMax, gu = 0.5)
{
  #Validity checking:
  if (gsi < 0.0 || gsi > 1.0)
  {
    stop(paste("Invalid GSI value:", gsi))
  }
  if (gu < 0.0 || gu > 1.0)#Values at the edges are unlikley too.
  {
    stop(paste("Invalid greenup threshold value:", gu))
  }
  
  m = (lfmMax - lfmMin) / (1.0 - gu)#Slope
  b = lfmMax - m#Intercept
  
  if (gsi < gu)
  {
    lfm = lfmMin
  }
  else#(gsi >= gu)
  {
    lfm = m * gsi + b
  }
  
  return(lfm)
}

#Live herbaceous fuel moisture:
#
#Parameters:
#gsi = The Growing Season Index value for the location.
#
#Returns: Percent live fuel moisture (water weight/dry fuel weight * 100%).
HerbaceousLiveFuelMoisture <- function(gsi)
{
  Min_H = 30#Minimum live fuel moisture parameter.
  Max_H = 250#Maximum live fuel moisture parameter.
  gu_H = 0.5#Default GSI greenup threshold.

  return(GSI_LiveFuelMoisture(gsi, Min_H, Max_H, gu_H))
}

#Live woody fuel moisture:
#
#Parameters:
#gsi = The Growing Season Index value for the location.
#
#Returns: Percent live fuel moisture (water weight/dry fuel weight * 100%).
WoodyLiveFuelMoisture <- function(gsi)
{
  Min_W = 60#Minimum live fuel moisture parameter.
  Max_W = 200#Maximum live fuel moisture parameter.
  gu_W = 0.5#Default GSI greenup threshold.
  
  return(GSI_LiveFuelMoisture(gsi, Min_W, Max_W, gu_W))
}
