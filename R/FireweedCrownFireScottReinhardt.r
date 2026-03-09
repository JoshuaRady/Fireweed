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
#Van Wagner, C. E.
#Conditions for the start and spread of crown fire.
#Canadian Journal of Forest Research 7(1): 23-34, 1977. https://doi.org/10.1139/x77-004
#
#Rothermel, R. C.
#Predicting behavior and size of crown fires in the Northern Rocky Mountains.
#Research Paper INT-RP-438. Ogden, UT: U.S. Department of Agriculture, Forest Service, Intermountain
#Research Station. 46 pages, 1991. https://doi.org/10.2737/INT-RP-438
#
#Scott, Joe H. and Reinhardt, Elizabeth D.
#Assessing crown fire potential by linking models of surface and crown fire behavior.
#Research Paper RMRS-RP-29. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky
#Mountain Research Station. 59 pages, 2001.  https://doi.org/10.2737/RMRS-RP-29
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
kmPerHrToMPerMin = 1000 / 60#1000 m/km / 60 min/hr = m/min
WindConversionFactor = 54.683

#Code:----------------------------------------------------------------------------------------------

#' Convert a fuel model to the physical properties of fuel model 10.
#' 
#' This is a helper function to deal with the fact that Rothermel 1991 performs calculations with
#' fuel model 10 regardless of the actual fuel model of the system being simulated.
#'
#' @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
#' fuel model.
#' @param fuelModel10 Fuel model 10 with default values.  Only needed if fm is not fuel model 10.
#' @return The converted fuel model.
ConvertToFuelModel10 <- function(fuelModel, fuelModel10)
{
  if (fuelModel$Number == 10)
  {
    return(fuelModel)
  }
  else
  {
    #Make sure the units are consistent:
    #We may want to force this to always be metric.
    if (fuelModel$Units != fuelModel10$Units)
    {
      FuelModelConvertUnits(fuelModel10, fuelModel$Units)
    }
    
    if (fuelModel$cured || fuelModel$NumClasses != 5)
    {
      Stop("Can't convert a fuel model with curing applied or more than the standard 5 classes.")
    }
    
    #Copy loadings:
    fuelModel10$w_o_ij = fuelModel$w_o_ij
    
    if (!"M_f_ij" %in% names(fuelModel))
    {
      Stop("M_f_ij must be provided in fuel model.")
    }
    fuelModel10$m_f_ij = fuelModel$m_f_ij
    
    return(fuelModel10)
  }
}

#Crown Fire Spread Rate:----------------------------------------------------------------------------

#' Calculate an estimate of the active crown fire spread rate based the method of Rothermel 1991.
#' 
#' Rothermel 1991 calculates the crown fire spread rate as simple multiple of the surface fire
#' spread rate but the surface fire spread rate must be calculated with the physical properties of
#' fuel model 10 and a fixed 40% wind reduction factor.  This function handle the needed
#' conversions and returns the resulting rate.
#' 
#' @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
#' fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
#' those of fuel model 10.
#' @param O Open wind speed at 6.1 m (km/hr)
#' @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
#' @param fuelModel10 Fuel model 10 with default values.  Only needed if fm is not fuel model 10.
#'
#' @returns The crown fire rate of spread (m/min).
SpreadRateCrownRothermel <- function(fuelModel, O, slopeSteepness, fuelModel10 = NULL)
{
  #Checks repeated from CrownFractionBurned():
  
  #Check for M_f_ij in the incoming fuel model.
  if (!"M_f_ij" %in% names(fuelModel))
  {
    Stop("M_f_ij must be provided in fuel model.")
  }
  
  #Check units:
  if (fuelModel$Units == "US")
  {
    fuelModel = FuelModelConvertUnits(fuelModel, "Metric")
  }
  
  #Check model type:
  if (fuelModel$Number != 10)
  {
    if (!is.null(fuelModel10))
    {
      Stop("Fuel model 10 needed.")
    }
    fuelModel = ConvertToFuelModel10(fuelModel, fuelModel10)
  }
  
  U = O * kmPerHrToMPerMin * 0.4#Use fixed 40% WRF from Rothermel 1991.
  R_surface = SpreadRateRothermelAlbini_HetFM(fuelModel, U, slopeSteepness)
  R_active = 3.34 * R_surface
  return(R_active)
}

#Crown Fire Initiation:-----------------------------------------------------------------------------

#' Return the critical surface fire intensity needed to initiate crowning per Van Wagner 1977.
#'
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100).
#'
#' @returns The critical surface (fireline) intensity for crowning (kW/m, I_0 in Van Wagner notation,
#'          I'_initiation in Scott & Reinhardt notation).
CriticalCrowningIntensityVanWagner <- function(CBH, FMC)
{
  #Scott & Reinhardt 2001 equation 11, pg. 13:
  #This is a combination of Van Wagner's equations 3, 4, and a value of C = 0.0010 from the text.
  IPrime_initiation = ((CBH * (460 + 25.9 * FMC)) / 100)^(3/2)
  return(IPrime_initiation)
}

#' Return the critical surface fire rate of spread needed to initiate crowning using Van Wagner 1977.
#'
#' @param IPrime_initiation (kW/m, I_0 in Van Wagner notation, I'_initiation in Scott & Reinhardt
#'                           notation).
#' @param HPA Surface fire heat per area (kw/m^2).
#'
#' @returns The critical surface fire rate of spread for crowning (m/min).
CriticalCrowningROSVanWagner <- function(IPrime_initiation, HPA)
{
  #Scott & Reinhardt 2001 equation 12, pg. 13:
  R_initiation = (60 * IPrime_initiation) / HPA
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

#' Calculate the torching index, the open wind speed at which torching starts,
#'
#' @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
#' family of functions.
#' @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
#'
#' @returns TI, the torching index (open wind speed at 6.1 m, km/hr).
TorchingIndex <- function(spreadCalcs, WRF, CBH, FMC)
{
  #Scott & Reinhardt 2001 equation 18, pg. 18:
  #TI = (1 / (54.683 * WRF)) * ((((60 * Iprime_initiation * rho_b * epsilon * Q_ig)/(HPA * xi * IR))
  #     - phi_s - 1) / (C * beta/beta_op)))^1/B
  #We break this equation up into the follow parts and compute them in steps:
  #TI = factor_L * (numerator_R / denominator_R)^1/B
  
  #We need the three wind factors:
  C = WindFactorC(spreadCalcs$cSAV, "Metric")
  B = WindFactorB(spreadCalcs$cSAV, "Metric")
  E = WindFactorE(spreadCalcs$cSAV, "Metric")
  
  #The leftmost factor:
  factor_L = 1 / (WindConversionFactor * WRF)
  
  #The numerator on the right is Scott & Reinhardt 2001 equation 16 pg. 18:
  #phiPrime_w_intiation = ((60 * IPrime_initiation * rho_b * epsilon * Q_ig)/(HPA * xi * I_R))
  #                       - phi_s - 1
  #Where;
  #rho_b = bulk density
  #epsilon = the effective heating number
  #Q_ig = heat of preignition
  #HPA = heat per area
  #xi = the propagating flux ratio
  #I_R = Reaction intensity
  #phi_s = slope factor
  
  #The numerator of equation 16 contains the full heat sink (rho_b * epsilon * Q_ig).
  #We substitute for this since we do not calculate epsilon directly for heterogeneous fuels:
  #phiPrime_w_intiation = ((60 * IPrime_initiation * HeatSink)/(HPA * xi * I_R)) - phi_s - 1
  
  IPrime_initiation = CriticalCrowningIntensityVanWagner(CBH, FMC)
  HPA = HeatPerUnitArea(spreadCalcs$I_R, ResidenceTime(spreadCalcs$cSAV, "Metric"))
  phiPrime_w_initiation = ((60 * IPrime_initiation * spreadCalcs$HeatSink) /
                             (HPA * spreadCalcs$xi * spreadCalcs$I_R)) - spreadCalcs$phi_s - 1
  
  #The denominator on the right consists of most of the Rothermel & Albini spread model wind factor
  #calculation.  It has been rearranged to remove U:
  #C * (beta / beta_opt)^-E
  denominator_R = C * (spreadCalcs$MeanPackingRatio / spreadCalcs$OptimumPR)^-E
  
  #The complete calculation:
  TI = factor_L * (phiPrime_w_initiation / denominator_R)^(1/B)
  
  return(TI)
}

#' Calculate the crowning index, the open wind speed at which active crowning occurs.
#'
#' @param spreadCalcs The spread calculations (list) output from the SpreadRateRothermelAlbini_*()
#' family of functions.
#' @param CBD Crown bulk density, the foliage (needles) and fine branches (kg/m^3).
#'
#' @returns CI, the crowning index (open wind speed at 6.1 m, km/hr).
CrowningIndex <- function(spreadCalcs, CBD)
{
  #Scott & Reinhardt 2001 equation 19, pg. 19 gives the midflame wind speed for active crowning.
  #UPrime_active = (1/54.683) * ((((3.0 / CBD) * rho_b * epsilon * Q_ig) / (3.34 * I_R * xi)
  #                               - spreadCalcs$phi_s - 1) / (C * (beta/beta_op)^-E))^1/B
  #To get the open wind speed for active crowning (O'_active = CI) using the assumptions of
  #Rothermel 1991 we can divide by a 40% wind reduction factor.
  #We break the resulting equation up into the follow parts and compute them in steps:
  #CI = factor_L * (numerator_R / denominator_R)^1/B
  
  #We need the three wind factors:
  C = WindFactorC(spreadCalcs$cSAV, "Metric")
  B = WindFactorB(spreadCalcs$cSAV, "Metric")
  E = WindFactorE(spreadCalcs$cSAV, "Metric")
  
  #The leftmost factor:
  factor_L = 1 / WindConversionFactor / 0.4
  
  #The numerator on the right is:
  #(((3.0 / CBD) * rho_b * epsilon * Q_ig) / (3.34 * I_R * xi)) - phi_s - 1
  #Which is the same as:
  #(((3.0 / CBD) * HeatSink) / 3.34 * I_R * xi) - phi_s - 1
  numerator_R = (((3.0 / CBD) * spreadCalcs$HeatSink) /
                   (3.34 * spreadCalcs$I_R * spreadCalcs$xi)) - spreadCalcs$phi_s - 1
  
  #The denominator on the right consists of most of the Rothermel & Albini spread model wind factor
  #calculation.  It has been rearranged to remove U  (It is the same as in the torching index.):
  #C * (beta / beta_opt)^-E
  denominator_R = C * (spreadCalcs$MeanPackingRatio / spreadCalcs$OptimumPR)^-E
  
  #The complete calculation:
  CI = factor_L * (numerator_R / denominator_R)^(1/B)
  
  #Note: Scott & Reinhardt 2001 uses equation 20, which is a simplified version of the above that
  #assumes parameters of fuel #model 10:
  #CI = 0.0457 * ((((3.0 / 3.34 / spreadCalcs$xi * spreadCalcs$HeatSink) /
  #                   (spreadCalcs$I_R * CBD))- spreadCalcs$phi_s - 1) / 0.001612)^0.7
  #We opt not to use the simplification to keep the code more parallel to the torching index
  #calculation and because rounding and small differences in conversion factors lead to slight
  #differences in the values we get with the explicit calculation.
  
  return(CI)
}

#' Calculate the crown fraction burned.
#' 
#' The representation in Scott & Reinhardt 2001 does not purport to actually predict the crown
#' fraction burned for its own sake.  Instead it is used as a transition function between surface
#' and active crown fire spread rate calculations.  However, this has not stopped people from using
#' the output as an approximation of actual CFB.
#'
#' The underlying Rothermel 1991 crown fire equations expect fuel model 10, timber litter and
#' understory.  We accept any fuel model, which will be converted internally if needed.  This
#' currently requires  fuel model 10 to be passed in, which is inelegant.  If the path to the
#' library's input files could be resolved this argument could be eliminated.
#'
#' The original model uses open wind speed at 6.1 m.  We add the option to use the midflame wind
#' speed as is in the Rothermel & Albini spread model.  In this case the open wind speed will be
#' calculated based on the wind reduction factor.  Only O or U should be provided and note that the
#' units differ for the two.
#'
#' @param fuelModel The fuel model representing the surface fuelbed.  M_f_ij must be included in the
#' fuel model.  If the fuel model it not fuel model 10 its physical properties will be converted to
#' those of fuel model 10.
#' @param O Open wind speed at 6.1 m (km/hr)
#' @param WRF Wind reduction factor.  Ratio to convert from open (6.1 m) to mid-flame wind speed.
#'            Needed even in U is supplied.
#' @param U Wind speed at midflame height (m/min).  Not needed if O is supplied.
#' @param slopeSteepness Slope steepness, maximum (unitless fraction: vertical rise / horizontal
#'                       distance).
#' @param CBD Crown bulk density, the foliage (needles) and fine branches (kg/m^3).
#' @param CBH Crown base height (m, z in original Van Wagner notation).
#' @param FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100)
#' @param fuelModel10 Fuel model 10 with default values.  Only needed if fm is not fuel model 10.
#'
#' @returns CFB, the crown fraction burned (fraction).
CrownFractionBurned <- function(fuelModel, O = NULL, WRF, U = NULL, slopeSteepness, CBD, CBH, FMC,
                                fuelModel10 = NULL)
{
  #Check for M_f_ij in the incoming fuel model.
  if (!"M_f_ij" %in% names(fuelModel))
  {
    Stop("M_f_ij must be provided in fuel model.")
  }
  
  #Check units:
  if (fuelModel$Units == "US")
  {
    fuelModel = FuelModelConvertUnits(fuelModel, "Metric")
  }
  
  #Check model type:
  if (fuelModel$Number != 10)
  {
    if (!is.null(fuelModel10))
    {
      Stop("Fuel model 10 needed.")
    }
    fuelModel = ConvertToFuelModel10(fuelModel, fuelModel10)
  }
  
  #Check the wind inputs:
  if (is.null(O) && is.null(U))
  {
    stop("Must provide wind speed as O or U.")
  }
  else if (!is.null(O))
  {
    #If O is passed in calculate U:
    U = O * WRF
    #Note: We don't actually need and accurate value for U as the spread calculations we use in
    #function are independent of the wind speed.  CI only depends on O.
  }
  else
  {
    #If U is passed in we need to calculate O:
    O = U / WRF
  }
  
  #Check the wind speed:
  if (O < 0)#or U
  {
    Stop("Invalid wind speed.")
  }
  #Add check for unexpectedly high wind speeds?
  
  #Scott & Reinhardt 2001 never mentions the wind limit and we assume it is not used:
  spreadCalcs = SpreadRateRothermelAlbini_HetFM(fuelModel, U, slopeSteepness, components = TRUE)
  
  TI = TorchingIndex(spreadCalcs, WRF, CBH, FMC)
  CI = CrowningIndex(spreadCalcs, CBD)
  
  if (O < TI)#Surface fire:
  {
    CFB = 0.0
  }
  else if (O >= CI)#Fully active crown fire:
  {
    CFB = 1.0
  }
  else#Passive crown fire, torching to full crown involvement:
  {
    #Calculate the transition as a linear function of wind speed:
    CFB = (1 / (CI - TI)) * (O - TI)
  }
  
  return(CFB)
}
