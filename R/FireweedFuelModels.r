#FireweedFuelModels.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 7, 19
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed fire code library.  It implements in R a fire behavior fuel
#model representation for use with the Rothermel Albini (Rothermel 1972, Albini 1976) fire spread
#model and related models.
#
#  The fule model object is a implemented as a R list with an fixed set of elements.  See
#GetFuelModelFromDF() for the list structure.
#___________________________________________________________________________________________________

source("FireweedUnits.r")
source("FireweedUtils.r")

#Constants:-----------------------------------------------------------------------------------------

#Fuel Category Constants:
#The values of the live dead categories are forced to match the matching array indexes so they may
#be used to access arrays of the form X_i (values are language specific):
Dead = 1
Live = 2

#Code:----------------------------------------------------------------------------------------------

#Find a fuel model in the specified file and return it as a fuel model object (list).
#
#Parameters:
#fuelModelFilePath = The path to the tab delimited file containing the table of fuel models.
#modelID = The standard fuel model number, alphanumeric code, or index of the model requested.  If
#  a number is passed and does not match a known model number it is interpreted as an index, that is
#  the position in the table of fuel models.  For 'the 13' the number, code, and index are the same.
#spreadModelUnits = If true then convert units used in the file that differ from those used in the
#  Rothermel & Albini spread model.
GetFuelModelFromTabDelimited <- function(fuelModelFilePath, modelID, spreadModelUnits = TRUE)
{
  fuelModelDF = read.delim(fuelModelFilePath, skip = 3)#The file has three lines of header.
  GetFuelModelFromDF(fuelModelDF, modelID, spreadModelUnits)
}

#Find a fuel model in the specified file and return it as a fuel model object (list).
#
#Parameters:
#fuelModelFilePath = The path to the CSV file containing the table of fuel models.
#modelID = The standard fuel model number, alphanumeric code, or index of the model requested.  If
#  a number is passed and does not match a known model number it is interpreted as an index, that is
#  the position in the table of fuel models.  For 'the 13' the number, code, and index are the same.
#spreadModelUnits = If true then convert units used in the file that differ from those used in the
#  Rothermel & Albini spread model.
GetFuelModelFromCSV <- function(fuelModelFilePath, modelID, spreadModelUnits = TRUE)
{
  fuelModelDF = read.csv(fuelModelFilePath, skip = 3)#The file has three lines of header.
  GetFuelModelFromDF(fuelModelDF, modelID, spreadModelUnits)
}

#Find a fuel model in the data frame passed return it as a fuel model object (list).
#
#Parameters:
#fuelModelDF = A data frame containing the table of fuel models.
#modelID = The standard fuel model number, alphanumeric code, or index of the model requested.  If
#  a number is passed and does not match a known model number it is interpreted as an index, that is
#  the position in the table of fuel models.  For 'the 13' the number, code, and index are the same.
#spreadModelUnits= If true then convert units used in the file that differ from those used in the
#  Rothermel & Albini spread model.
#
#ToDo:
# - The function assumes the input data is in its original English units if spreadModelUnits is true
#  but there handling for when this is false or is missing.  We could add an inputUnits parameter
#  to allow the input file to be in metric or alter the parameter behavior.
# - There is a question of whether to add M_f / M_f_ij to the data structure.
#
#Note: This expects draft 3 (D3) of the standard fuel models spreadsheet (loaded to a data frame).
#Note: This code currently assumes the units of the file are in United States customary units with
#loadings in ton/acre and moisture of extinction in percent.
GetFuelModelFromDF <- function(fuelModelDF, modelID, spreadModelUnits = TRUE)
{
  #Find and extract the row representing the fuel model requested:
  if (is.numeric(modelID))
  {
    #Check number for validity:
    if (length(modelID) != 1)
    {
      stop("GetFuelModel4() expects a single model ID:")
    }
    
    #Assume that a model index / row number is being requested if the number is low and not a
    #standard model number:
    if (modelID > 13 && modelID < 54)
    {
      theModelRow = fuelModelDF[modelID,]
    }
    else#Otherwise it must be a model number:
    {
      #Custom fuel models are possible so don't treat it as an error if it doesn't match one of
      #the 53 standard models:
      if (!(modelID %in% c(1:13, 101:109, 121:124, 141:149, 161:165, 181:189, 201:204)))
      {
        warning("Non-standard fuel model passed to GetFuelModel4().")
      }
      
      theModelRow = fuelModelDF[fuelModelDF$Number == modelID,]
    }
  }
  else if (is.character(modelID))
  {
    theModelRow = fuelModelDF[fuelModelDF$Code == modelID,]
  }
  else
  {
    stop("Expecting a fuel model code, number, or index.")
  }
  
  #Check that a matching row was found:
  if (nrow(theModelRow) == 0)
  {
    stop(paste("Model ID", modelID, "not found."))
  }
  
  #These parameters need to be rearranged into vectors:
  loadingCols = c("w_o_11", "w_o_12", "w_o_13", "w_o_21", "w_o_22")
  savCols = c("SAV_11", "SAV_12", "SAV_13", "SAV_21", "SAV_22")
  
  #Turn rows that do not need to be processed directly into members:
  #Unfortunately this doesn't work because of the column names being in a variable:
  #fuelModel = as.list(subset(theModelRow, select = -c(loadingCols, savCols)))
  #This ugly thing works:
  fuelModel = theModelRow[, !names(theModelRow) %in% c(loadingCols, savCols)]
  fuelModel = as.list(fuelModel)
  
  #Restructure fuel loadings and SAV into vectors:
  fuelModel$w_o_ij = as.numeric(subset(theModelRow, select = loadingCols))
  
  SAV_ij = as.numeric(subset(theModelRow, select = savCols))
  SAV_ij[is.na(SAV_ij)] = 0
  fuelModel$SAV_ij = SAV_ij
  
  #Typically the units of some parameters are different in the published tables than in the
  #model equations:
  #It would be an improvement to detect the units used in the file.  That information is currently
  #in the header information but not in a form that would be ideal to parse.
  if (spreadModelUnits)
  {
    fuelModel$w_o_ij = fuelModel$w_o_ij * lbsPerTon / ft2PerAcre#ton/acre to lb/ft^2
    fuelModel$M_x = fuelModel$M_x / 100#% to fraction
    
    #Record the units used:
    fuelModel$w_o_Units = "lbPer_ft2"
    fuelModel$M_x_Units = "Fraction"
  }
  else
  {
    fuelModel$w_o_Units = "tonPerAc"
    fuelModel$M_x_Units = "Percent"
  }
  
  #Expand parameters with fixed values across all fuel classes:
  #The standard values for these variables are whole numbers so they created as integers.  This
  #doesn't cause issues in R but is a problem when they are passed into the C++ spread equations.
  fuelModel$h = as.numeric(fuelModel$h)
  fuelModel$rho_p = as.numeric(fuelModel$rho_p)
  
  #It may be better to put zeros for classes not present:
  fuelModel$h_ij = rep(fuelModel$h, times = 5)
  fuelModel$S_T_ij = rep(fuelModel$S_T, times = 5)
  fuelModel$S_e_ij = rep(fuelModel$S_e, times = 5)
  fuelModel$rho_p_ij = rep(fuelModel$rho_p, times = 5)
  
  fuelModel$M_x_1 = fuelModel$M_x
  
  #Other:
  fuelModel$liveDead = as.integer(c(1,1,1,2,2))#Standard fuel models. Change to LiveDead?
  fuelModel$NumClasses = 5#length(fuelModel$liveDead)
  fuelModel$Units = "US"
  fuelModel$Cured = FALSE#Only relevant to dynamic fuel models.
  
  #Reorder members:
  numMembers = length(fuelModel)
  fuelModel = fuelModel[c("Number", "Code", "Name",#Fuel model identifiers:
                          #Model model properties:
                          "Type",#Static vs. Dynamic
                          "Units",#The model units type.
                          "Cured",
                          "NumClasses",
                          "w_o_Units",
                          "M_x_Units",
                          #Fuel model parameters / data members:
                          "SAV_ij",
                          "w_o_ij",
                          "delta",
                          "liveDead",
                          #For convenience homogeneous notation aliases are provided for the
                          #moisture of extinction and fuel particle properties. (Scalars could go?)
                          "M_x", "M_x_1",
                          "h", "h_ij",
                          "S_T", "S_T_ij",
                          "S_e", "S_e_ij",
                          "rho_p", "rho_p_ij",
                          #Precalculated columns (could be removed):
                          "CharacteristicSAV", "BulkDensity", "RelativePackingRatio")]
  
  #Make sure we didn't forget anything when reordering the members:
  if (length(fuelModel) != numMembers)
  {
    stop("Lost some members.")
  }
  
  return(fuelModel)
}

#Fuel Type Index Functions:-------------------------------------------------------------------------
#' @par Fuel Type Index Functions:
#' For standard fuel models the number and order of fuel types is fixed but there are complicating
#' factors in practice.  First, many of the fuel models only represent either one of the herbaceous
#' or woody fuel types.  Both index positions will be present but the underrepresented class should
#' not be used without making changes to the model.  Adding fuel loadings for these 'missing' size
#' classes should not prevent the Rothermel Albini spread model from running but the fuel will
#' effectively be ignored and no have impact on fire behavior output, though this is not well
#' tested.  Second, when dynamic fuel models are used the number of fuel types and indexes may
#' change.
#'
#' The public liveDead member can be used to infer the fuel categories and the FuelClassIndex()
#' function can be used to get indexes by category more conveniently.  Still additional assumptions
#' must be used to know exactly what each index represents.  Additionally we could change the way we
#' represent fuel types in the future.  For example we may add a live moss category in the future.
#'
#' To futureproof the code and to reduce the need for calling code to understand the way the fuel
#' types are organized the following accessor and information functions are provided.
#' 
#' @note These functions have not been fully integrated into the existing code.  A review is needed.

#' How many size classes are there in the dead fuel category?
#' 
#' @returns The number of dead fuel types for the model.
NumDeadClasses <- function(fm)
{
  return(sum(fm$liveDead == Dead))
}

#' How many size classes are there in the live fuel category?
#' 
#' @returns The number of live fuel types for the model.
NumLiveClasses <- function(fm)
{
  return(sum(fm$liveDead == Live))
}

#' Return the index (k) of the live herbaceous fuel type.
#' 
#' @param fm A fuel model object.
#' 
#' @returns The index (k) of the live herbaceous fuel type for the model.
LiveHerbaceousIndex <- function(fm)
{
  return(FuelClassIndex(fm$liveDead, Live, 1))#In the future the 1 index may not be guaranteed.
}

#' Return the index (k) of the live woody fuel type.
#' 
#' @param fm A fuel model object.
#' 
#' @returns The index (k) of the live woody fuel type for the model.
LiveWoodyIndex <- function(fm)
{
  return(FuelClassIndex(fm$liveDead, Live, 2))#In the future the 2 index may not be guaranteed.
}

#' Return the index (k) of the dead herbaceous fuel type, if present.  If not present returns -1.
#' Check cured first to see if it is present.
#'
#' @param fm A fuel model object.
#'
#' @returns The index (k) of the dead herbaceous fuel type, if present, -1 if not.
DeadHerbaceousIndex <- function(fm)
{
  if (fm$cured)
  {
    #We put the dead herbaceous fuel in the second dead position.  See CalculateDynamicFuelCuring().
    #This could change in the future.
    return(FuelClassIndex(fm$liveDead, Dead, 2))
  }
  else
  {
    return(-1)
  }
}

#' Is the live herbaceous fuel type active in this fuel model?
#'
#' @param fm A fuel model object.
#' 
#' @returns Is the live herbaceous fuel real of just a placeholder in the model?
LiveHerbaceousPresent <- function(fm)
{
  #We indicate that a live fuel is not present with a SAV of 0, which is unique to our implementation. 
  return(fm$SAV_ij[LiveHerbaceousIndex()] != 0)
}

#' Is the live woody fuel type active in this fuel model?
#'
#' @param fm A fuel model object.
#' 
#' @returns Is the live woody fuel type real of just a placeholder in the model?
LiveWoodyPresent <- function(fm)
{
  #We indicate that a live fuel is not present with a SAV of 0, which is unique to our implementation. 
  return(fm$SAV_ij[LiveWoodyIndex()] != 0)
}

#Fuel Moisture Functions:---------------------------------------------------------------------------

#Record the fuel moisture content values.
#
#Parameters:
#fm = The fuel model to apply curing to.
#M_f_ij = Fuel moisture content for for each fuel type (fraction: water weight/dry fuel weight).
#
#This function checks M_f_ij and stores it if valid, without apply curing.  If curing has been
#previously applied the  value can not be overwritten and an error is generated.  If the value has
#been previously set but curing has not been applied we allow the value to be overwritten.  The
#utility of overwriting is uncertain so we currently warn when this happens.
SetFuelMoisture <- function(fm, M_f_ij)
{
  functionName = match.call()[[1]]
  
  #Check the type:
  if (!is.numeric(M_f_ij))
  {
    stop(paste(functionName, "(): M_f_ij not of proper type", sep = ''))
  }
  
  #Check length is appropriate:
  if (length(M_f_ij) != fm$NumClasses)
  {
    stop(paste(functionName, "(): M_f_ij has length ", length(M_f_ij), " while numClasses = ",
               fm$NumClasses, sep = ''))
  }
  
  #Check values:
  for (i in 1:length(M_f_ij))
  {
    if (M_f_ij[i] < 0)#Clearly invalid:
    {
      stop(paste("M_f_ij contains invalid values:", paste(M_f_ij, collapse = ", ")))
    }
    else
    {
      #The maximum fuel moisture varies by fuel type but also by moisture model:
      #These limits could be researched further but making this a warning lowers the stakes.
      maxDeadFM = 1.5#100% may be fine but go higher.
      maxHerbFM = 2.5#GSI maxes out at 250%.
      maxWoodyFM = 2.0#GSI maxes out at 200%.
      
      if (fm$liveDead[i] == Dead && M_f_ij[i] > maxDeadFM)
      {
        warning("Dead fuel moisture value seems high: " + paste(M_f_ij, collapse = ", "))
      }
      else if (i == LiveHerbaceousIndex() && M_f_ij[i] > maxHerbFM)
      {
        warning("Herbaceous fuel moisture value seems high: " + paste(M_f_ij, collapse = ", "))
      }
      else if (i == LiveWoodyIndex() && M_f_ij[i] > maxWoodyFM)
      {
        warning("Woody fuel moisture value seems high: " + paste(M_f_ij, collapse = ", "))
      }
    }
  }
  
  #Checking if curing has already been run for this fuel model:
  if (cured == TRUE)
  {
    stop("Fuel model has already had curing applied.")
  }
  
  if (!this->M_f_ij.empty)
  {
    warning("M_f_ij is being overwritten.")
  }
  
  fm$M_f_ij = M_f_ij
  
  return(fm)
}

#Take a dynamic fuel model and calculate and apply the curing of herbaceous fuels based on the
#herbaceous fuel moisture (per Scott & Burgan 2005).  For dynamic fuels curing moves some live
#herbaceous fuel to a new dead herbaceous fuel class. As a result the number of fuel classes may
#increase with this call.  The function takes a fuel model object and returns a modified version of
#it reflecting the curing process.
#  Trying to apply curing to a static fuel model has no effect except to store the moisture values.
#By default the code posts a warning when this occurs.
#
#Parameters:
#fm = The fuel model to apply curing to.
#M_f_ij = Fuel moisture content for for each fuel type (fraction: water weight/dry fuel weight).
#  If omitted curing must be specified.
#curing = The percent herbaceous fuel curing to apply.  This is provided as an alternative to
#  specifying moisture content.  Curing is normally calculated based on the live herbaceous fuel
#  moisture but it can be useful to specify the curing directly, especially for testing.
#warn = If true warn about attempts to apply curing to static models.
#
#Notes:
#  On return Curing and possibly M_f_ij elements will be added to the fuel model (fm).  It may be
#better to add these on initialization of the model.  However in that case both will be undefined
#initially.  M_f_ij will be remain undefined for static fuel models or until it is set for dynamic
#fuel models.  With this implementation Curing and M_f_ij don't exist until this function is called.
#
#  This implementation results in different weights for missing fuels that those published in
#Andrews 2018, but this does not have material effects.
#  If created the dead herbaceous fuel is inserted at the second position (1,2).  For the standard
#fuel model this results in dead fuels being in descending SAV order with the exception of SH9.
CalculateDynamicFuelCuring <- function(fm, M_f_ij = NULL, curing = NULL, warn = TRUE)
{
  functionName = match.call()[[1]]
  
  fm = SetFuelMoisture(fm, M_f_ij)#Check validity of M_f_ij and save.
  
  if (fm$Type == "Dynamic")
  {
    #The live herbaceous is the first dead fuel.  The index should be 4 for standard fuel models:
    liveHerbIndex = match(Live, fm$liveDead)
    
    #Either M_f_ij or curing must be provided but not both:
    if (!is.null(M_f_ij))
    {
      #Curing is a function of live herbaceous fuel moisture:
      M_f_21 = M_f_ij[liveHerbIndex]
      
      #Calculate the transfer of herbaceous fuel loading from live to dead:
      #Curing is 0 at 120% fuel moisture and 1 at 30%:
      #Equation from Andrews 2018 pg. 36:
      #T = -1.11(M_f)_21 + 1.33, 0<=T<=1.0
      #This implements this exactly but is not numerically exact at the ends of the curing range:
      #cureFrac = -1.11 * M_f_21 + 1.33
      #Note: We can't use T, since T = TRUE in R.
      #This is exact at 30 and 120%:
      cureFrac = 1.0 - ((M_f_21 - 0.3) / 0.9)#1.2 - 0.3 = 0.9
      
      if (cureFrac < 0)
      {
        cureFrac = 0
        #Could break out here as no curing needs to be applied.
      }
      else if (cureFrac > 1)
      {
        cureFrac = 1
      }
    }
    else if (!is.null(curing))
    {
      #Curing is generally presented as a percentage and that is what we expect:
      if (curing < 0 || curing > 100)
      {
        stop(paste(functionName, "() expects fuel curing as a percentage.", sep = ''))
      }
      
      cureFrac = curing / 100
    }
    else
    {
      stop(paste(functionName, "() requires either M_f_ij or curing to be provided.", sep = ''))
    }
    
    #Checking if this function has already been run for this fuel model:
    if (fm$Cured == TRUE)
    {
      stop(paste(functionName, "(): Fuel model has already had curing applied.", sep = ''))
    }
    
    #Expand the number of fuel classes, inserting the cured herbaceous at the second dead position,
    #(Dead, 2) in 1 based space:
    fm$w_o_ij = c(fm$w_o_ij[1], 0, fm$w_o_ij[2:fm$NumClasses])#Initial loading = 0.
    #Curing doesn't change SAV.  Inherit from the live herbaceous:
    fm$SAV_ij = c(fm$SAV_ij[1], fm$SAV_ij[liveHerbIndex], fm$SAV_ij[2:fm$NumClasses])
    
    #For the standard fuel models all values for these parameters should be the same but don't
    #assume that for robustness.  Inherit from the live class:
    fm$h_ij = c(fm$h_ij[1], fm$h_ij[liveHerbIndex], fm$h_ij[2:fm$NumClasses])
    fm$S_T_ij = c(fm$S_T_ij[1], fm$S_T_ij[liveHerbIndex], fm$S_T_ij[2:fm$NumClasses])
    fm$S_e_ij = c(fm$S_e_ij[1], fm$S_e_ij[liveHerbIndex], fm$S_e_ij[2:fm$NumClasses])
    fm$rho_p_ij = c(fm$rho_p_ij[1], fm$rho_p_ij[liveHerbIndex], fm$rho_p_ij[2:fm$NumClasses])
    
    #Expand liveDead:
    fm$liveDead = c(fm$liveDead[1], Dead, fm$liveDead[2:fm$NumClasses])
    
    #Update the moisture content vector if provided:
    if (!is.null(M_f_ij))
    {
      fm$M_f_ij = c(M_f_ij[1], M_f_ij[1], M_f_ij[2:fm$NumClasses])#Inherit from 1-hr dead moisture.
    }
    # else
    # {
    #   fm$M_f_ij = NA or NULL
    # }
    
    fm$NumClasses = fm$NumClasses + 1
    liveHerbIndex = liveHerbIndex + 1#Update after all data members are restructured.
    
    #Transfer the loading from live to dead:
    fm$w_o_ij[2] = cureFrac * fm$w_o_ij[liveHerbIndex]
    #fm$w_o_ij[liveHerbIndex] = fm$w_o_ij[liveHerbIndex] - cureFrac * fm$w_o_ij[liveHerbIndex]
    fm$w_o_ij[liveHerbIndex] = fm$w_o_ij[liveHerbIndex] - fm$w_o_ij[2]#Same as above.
    
    fm$Cured = TRUE#Record that curing has been applied.
    fm$Curing = cureFrac * 100#Record the curing percentage.  Is this useful?
  }
  else if (warn)
  {
    warning("Fuel model is static. No curing applied.")
  }
  
  return(fm)
}

#Take a fuel model object and convert its units.
#
#Parameters:
#fm = The fuel model to convert.
#newUnits = The units to convert to.  If the same as the current units do nothing.
FuelModelConvertUnits <- function(fm, newUnits)#ConvertFuelModelUnits
{
  if (fm$Units == newUnits)
  {
    warning("Fuel model is already in requested units.")
  }
  else
  {
    if (newUnits == "US")
    {
      #Leave Number, Code, Name, and Type unchanged.
      #Cured and NumClasses don't change.
      
      fm$SAV_ij = fm$SAV_ij * cmPerFt#1/cm -> 1/ft | cm^2/cm^3 -> 3ft^2/ft^3
      
      #The metric units are alway kg/m^2. We don't include an option to convert to ton/Ac:
      fm$w_o_ij = fm$w_o_ij / kgPerLb * mPerFt^2#kg/m^2 -> lb/ft^2
      fm$w_o_Units = "lbPer_ft2"
      
      fm$delta = fm$delta / mPerFt#m -> ft
      
      #Leave liveDead unchanged.
      #M_x and M_x_1 are either unitless fractions or percentages and can be left unchanged.
      
      fm$h = fm$h / kJPerBtu * kgPerLb#kJ/kg -> Btu/lb
      fm$h_ij = fm$h_ij / kJPerBtu * kgPerLb#kJ/kg -> Btu/lb
      
      #S_T, S_T_ij, S_e, and S_e_ij are unitless fractions and can be left unchanged.
      
      fm$rho_p = fm$rho_p / lbPerFtCuToKgPerMCu
      fm$rho_p_ij = fm$rho_p_ij / lbPerFtCuToKgPerMCu
      
      fm$CharacteristicSAV = fm$CharacteristicSAV * cmPerFt#cm^2/cm^3 -> ft^2/ft^3
      
      fm$BulkDensity = fm$BulkDensity / lbPerFtCuToKgPerMCu#kg/m^3 -> lb/ft^3
      
      #RelativePackingRatio is a dimensionless ratio.
      #M_f_ij is a fraction.
      #Curing is a percent.
      
      fm$Units = "US"
    }
    else if (newUnits == "Metric")
    {
      #Leave Number, Code, Name, and Type unchanged.
      #Cured and NumClasses don't change.
      
      fm$SAV_ij = fm$SAV_ij / cmPerFt#1/ft -> 1/cm | ft^2/ft^3 -> cm^2/cm^3
      
      if (fm$w_o_Units == "tonPerAc")
      {
        fm$w_o_ij = fm$w_o_ij * tonsPerAcToLbPerSqFt#Convert loadings to lb/ft^2.
      }
      #We could check for invalid w_o_Units here.
      fm$w_o_ij = fm$w_o_ij * kgPerLb / mPerFt^2#lb/ft^2 -> kg/m^2
      fm$w_o_Units = "kgPer_m2"
      
      fm$delta = fm$delta * mPerFt#ft -> m
      
      #Leave liveDead unchanged.
      #M_x and M_x_1 are either unitless fractions or percentages and can be left unchanged.
      
      fm$h = fm$h * kJPerBtu / kgPerLb#Btu/lb -> kJ/kg
      fm$h_ij = fm$h_ij * kJPerBtu / kgPerLb#Btu/lb -> kJ/kg
      
      #S_T, S_T_ij, S_e, and S_e_ij are unitless fractions and can be left unchanged.
      
      fm$rho_p = fm$rho_p * lbPerFtCuToKgPerMCu
      fm$rho_p_ij = fm$rho_p_ij * lbPerFtCuToKgPerMCu
      
      fm$CharacteristicSAV = fm$CharacteristicSAV / cmPerFt#ft^2/ft^3 -> cm^2/cm^3
      
      fm$BulkDensity = fm$BulkDensity * lbPerFtCuToKgPerMCu#lb/ft^3 -> kg/m^3
      
      #RelativePackingRatio is a dimensionless ratio.
      #M_f_ij is a fraction.
      #Curing is a percent.
      
      fm$Units = "Metric"
    }
    else
    {
      stop("Requested fuel model units not recognized.")
    }
  }
  
  return(fm)
}

#' Return the index (k) in a variable array of the form X_ij given the (live/dead, size) index pair.
#'
#' Representing fuel model X_ij varaibles as vectors has the disadvantage of making mapping to
#' individual classes awkward.  This function makes this simple but remains somewhat inelegant.
#'
#' @param liveDead An array indicating if each index in the input variables represents a dead or
#'                 live fuel category.
#' @param liveDeadCat The Live / Dead catagory value to get.
#' @param sizeIndex The index (1 based) of the size class to get.
#'
#' @returns The index (k) of the requested fuel.
FuelClassIndex <- function(fm, liveDead, liveDeadCat, sizeIndex)
{
  numDead = sum(liveDead == Dead)
  numLive = length(liveDead) - numDead
  #Could add a check for invalid values here.
  
  if (liveDeadCat == Dead)
  {
    if (sizeIndex < 1 || sizeIndex > numDead)
    {
      stop("Invalid dead size index.")
    }
    else
    {
      return(sizeIndex)
    }
  }
  else if (liveDeadCat == Live)
  {
    if (sizeIndex < 1 || sizeIndex > numLive)
    {
      stop("Invalid live size index.")
    }
    else
    {
      return(numDead + sizeIndex)
    }
  }
  else
  {
    stop("Invalid live / dead category.")
  }
  
  return(-1)#We could make the size index errors above warnings and return.
}

