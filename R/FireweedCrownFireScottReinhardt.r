#FireweedCrownFireScottReinhardt.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2026
#Reference: Proj. 11 Exp. 26
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
#crown fire equations of Scott & Reinhardt 2001.
#
#  This code is in development and may be incomplete for a while.
#
#References:----------------------------------------------------------------------------------------
#
#Van Wagner 1977...
#Rothermel 1991 ...
#
#Assessing crown fire potential by linking models of surface and crown fire behavior.
#Scott, Joe H. and Reinhardt, Elizabeth D.
#Res. Pap. RMRS-RP-29. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky
#Mountain Research Station. 59 p. 2001.  https://doi.org/10.2737/RMRS-RP-29
#
#___________________________________________________________________________________________________

source("FireweedRAFireSpread.r")

#Constants:
#The Scott & Reinhardt 2001 model takes open wind speeds (O) in km/hr and converts them to midflame
#wind speeds (U) in order to perform component calculations from the Rothermel & Albini spread
#model.  They assume these calculations are done in the model's original English units.  Therefore
#they use a constant of 54.683 to convert wind speeds from km/hr to ft/min.  Fireweed can perform
#the spread calculations in United States customary units or metric units.  We choose to do
#everything in metric units so we replace the original constant with a conversion from km/hr to m/min:
#kmPerHrToMPerMin = 1000 / 60#1000 m/km / 60 min/hr
kmPerHrToMPerMin = 54.683
tempUnits = "Metric"#"US"

#Code:----------------------------------------------------------------------------------------------

#' Return an estimate of the active fire spread rate based on surface spread rate using the method
#' of Rothermel 1991.
#'
#' @param R_surface Surface fire rate of spread (m/min, for FM10).
#'
#' @returns The crown fire rate of spread (m/min).
#' 
#' @note The correlation between surface and crown fire ROS in Rothermel 1991 used fuel model 10
#' specifically.  ...
SpreadRateCrownRothermel <- function(R_surface)
{
  R_active = 3.34 * R_surface
}

#' Return an estimate of the active fire spread rate based on surface spread rate using the method
#' of Rothermel 1991.
#' 
#' The calculations per Rothermel 1991 used fuel model 10 for surface fire calculation specifically.
#' Is it useful to have a function devoted to doing the rate of spread calculations all in one step?
#' Since we have to pass in fuel model 10 and the other inputs.  I'm not sure it helps anything.
#' We could have a function to check it? 
CrownSpreadRateRothermel <- function()
{
  
}

#Crown Fire Initiation:-----------------------------------------------------------------------------

#' Return the critical surface intensity needed to initiate crowning using Van Wagner 1977.
#'
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (% water weight/dry fuel weight x 100).
#'
#' @returns The critical surface (fireline) intensity for crowning (kW/m, I_0 in Van Wagner notation,
#'          I'_initiation in Scott & Reinhardt notation).
CrownFireInitiationVanWagner <- function(CBH, FMC)#CriticalCrowningIntensityVanWagner
{
  #Scott & Reinhardt 2001 equation 11, pg. 13:
  #This is a combination of Van Wagner's equations 3, 4, and value of C = 0.0010 from the text.
  IPrime_initiation = ((CBH * (460 + 25.9 * FMC)) / 100)^(3/2)#(2/3)
  return(IPrime_initiation)
}

#' Return the critical surface fire rate of spread needed to initiate crowning using Van Wagner 1977.
#'
#' @param I_initiation (kW/m, I_0 in Van Wagner notation, I'_initiation in Scott & Reinhardt notation).
#' @param HPA Surface fire heat per area (kw/m^2).
#'
#' @returns The critical surface fire rate of spread for crowning (m/min).
CriticalCrowningROSVanWagner <- function(I_initiation, HPA)
{
  #Scott & Reinhardt 2001 equation 12, pg. 13:
  R_initiation = (60 * I_initiation) / HPA
  return(R_initiation)
}

#Active Crown Fire:

#' Calculate the (minimum) critical crown fire rate of spread for active crowning.
#'
#' @param CBD Crown bulk density (kg/m^3).  This is the dry mass of foliage and fine branches per
#'            volume in the tree canopy.
#'
#' @returns The (minimum) critical crown fire rate of spread for active crowning (m/min).
CriticalActiveROSVanWagner <- function(CBD)
{
  #Scott & Reinhardt 2001 equation 14, pg. 14:
  Rprime_active = 3.0 / CBD
  return(Rprime_active)
}

#' Calculate the torching index, the open wind speed?? at which torching starts
#'
#' @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
#' family of functions.
#' @param WRF Wind reduction factor.  Ratio used to convert from mid-flame to open wind speed.
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
#'
#' @returns TI, the torching index (open wind speed at 6.1 m, km/hr).
TorchingIndex <- function(spreadCalcs, WRF, CBH, FMC)
{
  #Scott & Reinhardt 2001 equation 18, pg. 18:
  #TI = (1 / (54.683 * WRF)) * ((((60 * Iprime_initiation * rho_b * epsilon * Q_ig)/(HPA * xi * IR))
  #     - phi_s - 1) / (C * beta/beta_op)))^1/B
  #We break this equation up into the follow parts for readability:
  #TI = factor_L * (numerator_R / denominator_R)^1/B
  
  #We need the three wind factors:
  C = WindFactorC(spreadCalcs$cSAV, tempUnits)
  B = WindFactorB(spreadCalcs$cSAV, tempUnits)
  E = WindFactorE(spreadCalcs$cSAV, tempUnits)
  #WRF = wind reduction factor, the ratio of midflame wind speed to open wind speed...
  
  #The leftmost factor:
  #factor_L = 1 / (54.683 * WRF)
  factor_L = 1 / (kmPerHrToMPerMin * WRF)
  
  #The numerator on the right is Scott & Reinhardt 2001 equation 16 pg. 18:
  #phiPrime_w_intiation = ((60 * IPrime_initiation * HeatSink)/(HPA * xi * I_R)) - phi_s - 1
  #[Bulk density]?????
  
  #The numerator of equation 16 contains the full heat sink (rho_b * epsilon * Q_ig) where
  #rho_b = Bulk density
  #...
  #and:
  #HPA = heat per area
  #xi = Effective heating number
  #I_R = Reaction intensity
  #phi_s = slope factor
  
  IPrime_initiation = CrownFireInitiationVanWagner(CBH, FMC)
  HPA = HeatPerUnitArea(spreadCalcs$I_R, ResidenceTime(spreadCalcs$cSAV, tempUnits))
  
  phiPrime_w_initiation = ((60 * IPrime_initiation * spreadCalcs$HeatSink) /
                             (HPA * spreadCalcs$xi * spreadCalcs$I_R)) - spreadCalcs$phi_s - 1
  
  #The denominator on the right consists of most of the terms from the wind factor calculation.  It
  #has been rearranged to remove U.
  #C * (beta / beta_opt)^-E
  #beta = mean packing ratio
  #beta_op = optimum packing ratio
  #beta/beta_op = relative packing ratio
  denominator_R = C * (spreadCalcs$MeanPackingRatio / spreadCalcs$OptimumPR)^-E
  
  #The complete calculation:
  TI = factor_L * (phiPrime_w_initiation / denominator_R)^(1/B)#1/B
  
  return(TI)
}

#' Calculate the crowning index, the open wind speed?? at which crowning starts.
#'
#' @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
#' family of functions.
#' @param CBD Crown bulk density, the foliage (needles) and fine branches (kg/m^3).
#'
#' @returns CI, the crowning index (open wind speed at 6.1 m, km/hr).
CrowningIndex <- function(spreadCalcs, CBD)
{
  #This is Scott & Reinhardt 2001 equation 19, pg. 19:
  #...
  #We break this equation up into the follow parts for readability:
  #CI = factor_L * (numerator_R / denominator_R)^1/B
  #Note: Scott & Reinhardt 2001 equation 20 is a simplifed version that assumes parameters of fuel
  #model 10.  We use equation 19 for flexibility.  It allows the user to break this assumption.
  
  #We need the three wind factors:
  #tempUnits = "US"
  C = WindFactorC(spreadCalcs$cSAV, tempUnits)
  B = WindFactorB(spreadCalcs$cSAV, tempUnits)
  E = WindFactorE(spreadCalcs$cSAV, tempUnits)
  
  #The leftmost factor:
  #factor_L = 1 / 54.683
  #factor_L = 1 / kmPerHrToMPerMin
  factor_L = 1 / 54.683 / 0.4
  
  #The numerator on the right is:
  #(((3.0 / CBD) * rho_b * epsilon * Q_ig) / 3.34 * I_R * xi) - phi_s - 1
  #Which is the same as:
  #(((3.0 / CBD) * HeatSink) / 3.34 * I_R * xi) - phi_s - 1
  # numerator_R = (((3.0 / CBD) * spreadCalcs$HeatSink) /
  #                  3.34 * spreadCalcs$I_R * spreadCalcs$xi) - spreadCalcs$phi_s - 1
  numerator_R = (((3.0 / CBD) * spreadCalcs$HeatSink) /
                   (3.34 * spreadCalcs$I_R * spreadCalcs$xi)) - spreadCalcs$phi_s - 1
  
  #The denominator on the right consists of most of the terms from the wind factor calculation.  It
  #has been rearranged to remove U.  (It is the same as in the torching index.)
  #C * (beta / beta_opt)^-E
  #beta = mean packing ratio
  #beta_op = optimum packing ratio
  #beta/beta_op = relative packing ratio
  denominator_R = C * (spreadCalcs$MeanPackingRatio / spreadCalcs$OptimumPR)^-E
  
  #The complete calculation:
  CI = factor_L * (numerator_R / denominator_R)^(1/B)#1/B
  #CI = (1/54.683 / 0.4) * (numerator_R / denominator_R)^(1/B)#1/B
  #CI = 0.0457 * (numerator_R / denominator_R)^(1/B)#1/B
  
  #Equation 20:
  # CI2 = 0.0457 * ((((3.0 / 3.34 / spreadCalcs$xi * spreadCalcs$HeatSink) /
  #                     (spreadCalcs$I_R * CBD))- spreadCalcs$phi_s - 1) / 0.001612)^0.7
  
  return(CI)
}

#' Calculate the crown fraction burned.
#' 
#' The representation in Scott & Reinhardt 2001 does not purport to actually predict the crown
#' fraction burned for its own sake.  Instead it is used as a transition function between surface
#' and active crown fire spread rate calculations.  However, this has not stopped people from using
#' the output as an approximation of actual CFB and we do so well.
#'
#' @param fm10 Fuel model 10 with fuel loadings. ... 
#' M_f_ij must be included in the fuel model.
#' Expect FM 10 and change it if not?  We shoud have an option to not change it.  No, we still need
#' the data from fuel model 10 so it needs to be passed in.
#' 
#' We need to decide on how to deal with the wind.  We only need U or O and WRF.  Decide on units too.
#' #I guess O = U when WRF = 1????
#' @param U Wind speed at midflame height (m/min).
#' @param O Open wind speed at 6.1 m (km/hr)
#' @param WRF Wind reduction factor.  Ratio used to convert from mid-flame to open wind speed.
#' 
#' @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
#'                       distance).
#' 
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (% water weight/dry fuel weight x 100)
#' @param CBD Crown bulk density, the foliage (needles) and fine branches (kg/m^3).
#'
#' @returns CFB, the crown fraction burned (fraction).
CrownFractionBurned <- function(fm10, #U,
                                O, WRF,
                                slopeSteepness,
                                #spreadCalcs,
                                CBH, FMC, CBD)
                                #useWindLimit = FALSE)#This probably has an assumed value?????
{
  #Check for M_f_ij in the incoming fuel model.
  
  #Check the wind speed:
  if (O < 0)#or U
  {
    Stop("Invalid wind speed.")
  }
  #Add check for unexpectedly high wind speeds?
  
  #If O is passed in we need to calculate U...
  #U = O * WRF
  U = O * WRF * kmPerHrToMPerMin
  #U = O * 0.4#See pg. 11...
  
  spreadCalcs = SpreadRateRothermelAlbini_HetFM(fm10, U, slopeSteepness, components = TRUE)
  
  TI = TorchingIndex(spreadCalcs, WRF, CBH, FMC)
  CI = CrowningIndex(spreadCalcs, CBD)
  
  if (O < TI)#Surface fire
  {
    CFB = 0.0
  }
  else if (O >= CI)#Fully active crown fire
  {
    CFB = 1.0
  }
  else#Torching to full crown involvement:
  {
    #Calculate the transition as a linear function of wind speed:
    CFB = (1 / (CI - TI)) * (O - TI)
  }
  
  return(CFB)
}
