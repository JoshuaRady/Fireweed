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
#___________________________________________________________________________________________________

#Find a fuel model in the specified file and return it as a fuel model object (list).
#
#Parameters:
#modelID = The standard fuel model number, alphanumeric code, or index of the model requested.  If
#  a number is passed and does not match a known model number it is interpreted as an index, that is
#  the position in the table of fuel models.  For 'the 13' the number, code, and index are the same.
#fuelModelPath = The path to the CSV file containing the table of fuel models.
#spreadModelUnits= If true then convert units used in the file that differ from those used in the
#Rothermel & Albini spread model.
#
#ToDo:
# - The function assumes the input data is in its original English units if spreadModelUnits is true
#  but there handling for when this is false or is missing.  We could add an inputUnits parameter
#  to allow the input file to be in metric or alter the parameter behavior.
# - There is a question of whether to add M_f / M_f_ij to the data structure.
#
#Note: This expects draft 3 (D3) of the standard fuel models spreadsheet.
#Note: This code currently assumes the units of the file are in United States customary units with
#loadings in ton/acre and moisture of extinction in percent.
GetFuelModelFromCSV <- function(modelID, fuelModelPath, spreadModelUnits = TRUE)
{
  fuelModelDF = read.delim(fuelModelPath, skip = 3)#The file has three lines of header.
  
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
  #It would be an improvement to detect the units used the file.  That information is currently
  #in the header information but not in a form that would be ideal to parse.
  if (spreadModelUnits)
  {
    fuelModel$w_o_ij = fuelModel$w_o_ij * lbsPerTon / ft2PerAcre#ton/acre to lb/ft^2
    fuelModel$M_x = fuelModel$M_x / 100#% to fraction
    
    #Record the units used:
    w_o_Units = "lbPer_ft2"
    M_x_Units = "Fraction"
  }
  else
  {
    w_o_Units = "tonPerAc"
    M_x_Units = "Percent"
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
  fuelModel = fuelModel[c("Number", "Code", "Name",#Fuel model identifiers:
                          #Model model properties:
                          "Type",#Static vs. Dynamic
                          "Units",#The model units type.
                          "Cured",
                          "NumClasses",
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
  
  return(fuelModel)
}

#Copied from Proj_11_Exp_7_Analysis.r:
##6/2/2023:
#Take a dynamic fuel model and calculate and apply the curing of herbaceous fuels based on the
#herbaceous fuel moisture (per Scott & Burgan 2005).  The function takes a fuel model object and
#returns a modified version of it reflecting the curing process.
#Changes: Add optional curing parameter as second way to specify curing.
ApplyDynamicFuelCuring3 <- function(fm, M_f_ij = NULL, curing = NULL, warn = TRUE)
{
  functionName = match.call()[[1]]
  
  if (fm$Type == "Dynamic")
  {
    #The live herbaceous is the first dead fuel.  The index should be 4 for standard fuel models:
    liveHerbIndex = match(2, fm$liveDead)
    
    #Either M_f_ij or curing must be provided but not both:
    if (!is.null(M_f_ij))
    {
      #Check length is appropriate:
      if (!is.numeric(M_f_ij) || length(M_f_ij != fm$NumClasses))
      {
        stop(paste(functionName, "(): M_f_ij not of proper form.", sep = ''))
      }
      
      #Curing is a function of live herbaceous fuel moisture:
      M_f_21 = M_f_ij[liveHerbIndex]
      
      #Calculate the transfer of herbaceous fuel loading from live to dead:
      #Curing is 0 at 120% fuel moisture and 1 at 30%:
      #Equation from Andrews 2018 pg. 36:
      #T = -1.11(M_f)_21 + 1.33, 0<=T<=1.0
      #This implements this exactly but is not numerically exact:
      #Note: We can't use T, since T = TRUE in R.
      #cureFrac = -1.11 * M_f_21 + 1.33
      #This is exact at 30 and 120%:
      cureFrac = (M_f_21 - 0.3) / 0.9#1.2 - 0.3 = 0.9
      
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
        stop(paste(functionName, "() expects curing as a percentage.", sep = ''))
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
    
    #liveHerbIndex = match(2, fm$liveDead)#Should be 4 for standard fuel model.
    
    #Expand the number of fuel classes, inserting the cured herbaceous at position 1,2:
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
    fm$liveDead = c(fm$liveDead[1], 1, fm$liveDead[2:fm$NumClasses])
    
    #Curing is a function of live herbaceous fuel moisture:
    #M_f_21 = M_f_ij[liveHerbIndex]
    
    #Update the moisture content vector if provided:
    if (!is.null(M_f_ij))
    {
      #I'm not sure about how to return this.  Adding it to the fuel model makes some sense but
      #M_f_ij will be undefined for static fuel models or where it is not passed in until it is
      #added to the fuel model upstream.
      fm$M_f_ij = c(M_f_ij[1], M_f_ij[1], M_f_ij[2:fm$NumClasses])#Inherit from 1-hr dead moisture.
    }
    # else
    # {
    #   fm$M_f_ij = NA or NULL
    # }
    
    fm$NumClasses = fm$NumClasses + 1
    liveHerbIndex = liveHerbIndex + 1#Update after all date members are revised.
    
    #Calculate the transfer of herbaceous fuel loading from live to dead:
    #Curing is 0 at 120% fuel moisture and 1 at 30%:
    #Equation from Andrews 2018 pg. 36:
    #T = -1.11(M_f)_21 + 1.33, 0<=T<=1.0
    #This implements this exactly but is not numerically exact:
    #Note: We can't use T, since T = TRUE in R.
    #transFrac = -1.11 * M_f_21 + 1.33
    #This is exact at 30 and 120%.
    # transFrac = (M_f_21 - 0.3) / 0.9#1.2 - 0.3 = 0.9  transFrac -> cureFrac?
    # 
    # if (transFrac < 0)
    # {
    #   transFrac = 0
    #   #Could break out here as no curing needs to be applied.
    # }
    # else if (transFrac > 1)
    # {
    #   transFrac = 1
    # }
    
    #Transfer the loading from live to dead:
    fm$w_o_ij[2] = cureFrac * fm$w_o_ij[liveHerbIndex]#w_o_ij
    #fm$w_o_ij[liveHerbIndex] = fm$w_o_ij[liveHerbIndex] - cureFrac * fm$w_o_ij[liveHerbIndex]
    fm$w_o_ij[liveHerbIndex] = fm$w_o_ij[liveHerbIndex] - fm$w_o_ij[2]#Same as above.
    
    fm$Cured = TRUE#Record that curing has been applied.
    #It could be better to record the curing fraction, which = cureFrac.
    #The fraction cured could also be taken as an alternative argument to fuel moisture.  There
    #could be conditions when it would be useful to specify that directly.
    fm$Curing = cureFrac * 100#Is this useful?
  }
  else if (warn)
  {
    warning("Fuel model is static. No curing applied.")
  }
  
  return(fm)
}
