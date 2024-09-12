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
#
#Modified slight from Proj_11_Exp_20_Analysis.r GrowingSeasonIndexD2().
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
  
  GSI = iT_Min * iVPD * iPhoto
  return(GSI)
}
