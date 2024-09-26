#FireweedDeadFuelMoistureFosberg.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 20
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed wildfire code library.  It contains an implementation of the
#Fosberg dead fuel moisture model, specifically the Rothermel 1981 / NWCG variant.
#
#References:----------------------------------------------------------------------------------------
#
#Fosberg, M.A. and Deeming, J.E.
#Derivation of the 1-and 10-hour timelag fuel moisture calculations for fire-danger rating.
#Research Note RM-RN-207. Rocky Mountain Forest and Range Experiment Station, Fort Collins, CO.
#Forest Service, US Department of Agriculture, 1971.
#  This is the first publication documenting the 'Fosberg' dead fuel moisture calculation.  It gives
#calcuations for 1-hour and 10-hour fuels.  Subsequent publications present methods for 100-hour and
#1000-hour fuels.  The Fosberg method was highly influential and is still used today, with some
#modifications.
#
#Rothermel, Richard C.
#How to predict the spread and intensity of forest and range fires.
#Gen. Tech. Report INT-143. Ogden, UT: U.S. Department of Agriculture, Forest Service, Intermountain
#Forest and Range Experiment Station. 161 p., 1983.
#  This presents a modified version of the tables presented in Fosberg and Deeming 1971 for use in
#field estimation of 1-hour fuel moisture.  These tables remain in use today by the National
#Wildfire Coordinating Group (NWCG) and are included in the NWCG Incident Response Pocket Guide
#(IPTG, as of 2022) carried by wildland firefighters.
#
#___________________________________________________________________________________________________

source("FireweedUnits.r")
source("FireweedUtils.r")

#Code:----------------------------------------------------------------------------------------------

#Return an estimate of the 1-hour dead fuel moisture based on the conditions passed.  The modified
#Fosberg lookup table method presented in Rothermel 1981 and adopted by the NWCG is used.
#
#Parameters:
#tableA_Path - tableD_Path = Paths to tab delimited files holding the lookup table data.
#  These tables are used in Rothermel 1981 and are also the NWCG website and the 2022 IPTG.
#  These paths add to an already large number of parameters.  Their file location can not be safely
#  assumed but there values could be consolidated in a list of placed in globals.
#
#temp = The air temperature at 4.5 feet / 1.37 meters above ground level (degrees F / C).
#rh = Relative humidity (percent).
#monthOfYear = The numeric month of year (1 - 12).
#hourOfDay = The numeric hour of day in 24 hour format (1 - 24).
#slopePct = Percent slope of the location (rise / run x 100%).
#aspect = The aspect of the local slope as a cardinal direction (N, E, S, W passed as a character)
#  or a bearing (degrees from North, 0-360).
#shaded = Does the location have 50% or greater canopy cover or is there full cloud cover.
#  This could be converted to a percentage.
#elevation = A code indicating the slope position for the prediction relative to where the weather
#  conditions were taken.  Values are:
#  B: The fire is 1000 - 2000 feet below the location where weather was recorded,
#  L: The fire is within 1000 feet of elevation from the the location where weather was recorded,
#  A: The fire is 1000 - 2000 feet above the location where weather was recorded,
#  This is for field use only and may be omitted for other applications.
#units: The units to use.  Only relevant to temp.
#
#Returns: 1-hour fuel moisture (fraction: water weight/dry fuel weight).
#
#Based on Proj_11_Exp_20_Analysis.r FosbergNWCG_Table1HrFM2().
FosbergNWCG_1HrFM <- function(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
                              temp, rh, monthOfYear, hourOfDay, slopePct, aspect,
                              shaded = FALSE, elevation = "L", units = "Metric")
{
  #Convert the temperature if needed:
  if (units == "US")
  {
    tempF = temp
  }
  else#units == "Metric"
  {
    tempF = CtoF(temp)
  }
  
  rfm = FosbergNWCG_GetRFM(tableA_Path, tempF, rh)#Look up the reference fuel moisture.
  
  if (monthOfYear %in% 5:7)#May - July
  {
    correctionTablePath = tableB_Path
  }
  else if (monthOfYear %in% c(2:4, 8:10))#Feb - April, Aug - Oct
  {
    correctionTablePath = tableC_Path
  }
  else if (monthOfYear %in% c(1, 11, 12))#Nov - Jan
  {
    correctionTablePath = tableD_Path
  }
  else
  {
    stop("Invalid month value. Please supply the integer month of year.")
  }
  
  #If the aspect is input as degrees check and convert before passing it on:
  aspectErrStr = "Invalid aspect. Aspect must be a cardinal direction (N, E, S, W) or bearing in degrees."
  
  if (is.numeric(aspect))
  {
    #This is imperfect as NE = 45 degrees must be lumped into either N or E, and so on.  We lump
    #counterclockwise:
    if (aspect <= 45)
    {
      aspectCardinal = "N"
    }
    else if (aspect <= 135)
    {
      aspectCardinal = "E"
    }
    else if (aspect <= 225)
    {
      aspectCardinal = "S"
    }
    else if (aspect <= 315)
    {
      aspectCardinal = "W"
    }
    else if (aspect <= 360)
    {
      aspectCardinal = "N"
    }
    else
    {
      stop(aspectErrStr)
    }
  }
  else if (class(aspect) != "character")
  {
    stop(aspectErrStr)
  }
  #Values are checked in FosbergNWCG_GetCorrection().
  
  #Look up the correction factor for the conditions specified:
  correction = FosbergNWCG_GetCorrection(correctionTablePath, hourOfDay, slopePct, aspectCardinal,
                                         shaded, elevation)
  
  return(rfm + correction)#Combine and return.
}

#Look up the Fosberg reference fuel moisture (NWCG variant) given the temperature and humidity:
#
#Parameters:
#tableA_Path = A path to a tab delimited file holding lookup table A from Rothermel 1981 & the NWCG.
#
#tempF = The air temperature at 4.5 feet / 1.37 meters above ground level (degrees F).
#rh = Relative humidity (percent).
#
#Returns: Reference fuel moisture (fraction: water weight/dry fuel weight).
#
#Modified from Proj_11_Exp_20_Analysis.r FosbergNWCG_GetRFM().
FosbergNWCG_GetRFM <- function(tableA_Path, tempF, rh)
{
  #Parameter checking:
  #I need to determine a valid temp range.  The tables have a bottom but not a top.
  if (!InRange(rh, 0, 100))
  {
    stop("Relative humidity must be a valid percentage")
  }
  
  #Read in the header IDs as data to avoid conversion to prepended strings:
  df = read.delim(tableA_Path, header = FALSE)
  
  #Break into IDs and data:
  
  #Relative humidity column IDs:
  rhRangeBottoms = df[1, -1]
  rhRangeBottoms = as.numeric(as.vector(rhRangeBottoms))
  
  #Temperature row IDs:
  tempRangeBottoms = df[-1, 1]
  
  #Extract the lookup table values:
  luDF = df[-1, -1]#A dataframe we will treat as a matrix.
  
  #Search for the matching value:
  numTempBins = length(tempRangeBottoms)
  for (i in 1:numTempBins)
  {
    if (i < numTempBins)
    {
      #The top bound for each bin is the bottom of the next:
      if (tempF >= tempRangeBottoms[i] && tempF < tempRangeBottoms[i + 1])
      {
        tempIndex = i
        break
      }
    }
    else
    {
      #The last bin only has a lower temperature bound:
      if (tempF >= tempRangeBottoms[i])
      {
        tempIndex = i
      }
      else#Should not possible with a valid values.
      {
        stop("Match not found.")
      }
    }
  }
  
  numRHBins = length(rhRangeBottoms)
  for (j in 1:numRHBins)
  {
    if (j < numRHBins)
    {
      #The top bound for each bin is the bottom of the next:
      if (rh >= rhRangeBottoms[j] && rh < rhRangeBottoms[j + 1])
      {
        rhIndex = j
        break
      }
    }
    else
    {
      #The value must match the final bin but we check anyway:
      if (rh == rhRangeBottoms[j])
      {
        rhIndex = j
      }
      # else#Should not possible with a valid values.
      # {
      #   stop("Match not found.")
      # }
    }
  }
  
  return(luDF[tempIndex, rhIndex])
}

#Look up the Fosberg fuel moisture correction (NWCG variant) for the conditions passed:
#Note: This doesn't handle the nighttime values from table C or D.
#
#Parameters:
#tableFilePath = A path to a tab delimited file holding a correction lookup table (B - D) from
#  Rothermel 1981 & the NWCG.
#
#hourOfDay = The numeric hour of day in 24 hour format (1 - 24).
#slopePct = Percent slope of the location (rise / run x 100%).
#aspectCardinal = The slope aspect as a cardinal direction (N, E, S, W).
#  Note: The ability to take the aspect as degrees is enabled in the calling code..
#shaded = Does the location have 50% or greater canopy cover or is there full cloud cover.
#  Note: We could also accept a percentage.
#elevation = A code indicating the slope position for the prediction relative to where the weather
#  conditions were taken.  Values are:
#  B: The fire is 1000 - 2000 feet below the location where weather was recorded,
#  L: The fire is within 1000 feet of elevation from the the location where weather was recorded,
#  A: The fire is 1000 - 2000 feet above the location where weather was recorded,
#  This is for field use only and may be omitted for other applications.
#
#Returns: Fuel moisture correction (fraction: water weight/dry fuel weight).
#
#Modified from Proj_11_Exp_20_Analysis.r FosbergNWCG_GetCorrection2().
FosbergNWCG_GetCorrection <- function(tableFilePath, hourOfDay, slopePct, aspectCardinal, shaded,
                                      elevation = "L")
{
  #Read in the header IDs as data to avoid conversion to prepended strings:
  df = read.delim(tableFilePath, header = FALSE)
  
  #Break into IDs and data.  There are 3 column IDs and 3 row IDs:
  
  #Column IDs:
  idCols = 1:3
  
  #The column IDs give start and end times for each column in military time hours and minutes without
  #a colon.  Convert to hour of day:
  hourStart = as.numeric(as.vector(df[1, -idCols])) / 100
  hourEnd = as.numeric(as.vector(df[2, -idCols])) / 100
  
  #Parameter checking (after extracting hours):
  if (!InRange(hourOfDay, min(hourStart), max(hourEnd)))
  {
    stop("Invalid hour of day value or nighttime value passed.")
  }
  #The problem with percent slope is that is becomes huge as it approaches 90 degrees.  However,
  #realistically very high slopes should be very rare.  This check is minimal.
  if (!InRange(slopePct, 0, 1000))#0 - ~85 degrees.
  {
    stop("Invalid percent slope")
  }
  if (!(aspectCardinal %in% c("N", "S", "E", "W")))
  {
    stop("Invalid aspect.  Must be a cardinal direction (N, E, S, W).")
  }
   if (!elevation %in% c("B", "L", "A"))#Could allow lowercase.
  {
    stop("Invalid relative elevation code.")
  }
  
  elevationCode = as.character(as.vector(df[3, -idCols]))
  
  #Row IDs:
  headerRows = 1:3
  shadeID = df[-headerRows, 1]
  aspects = df[-headerRows, 2]
  slopeClass = df[-headerRows, 3]
  
  #Extract the lookup table values:
  luDF = df[-headerRows, -idCols]
  #Since the the elevation header is character all the columns will be interpreted as character:
  luDF = data.frame(apply(luDF, 2, as.numeric))#Convert to numeric.
  #luMatrix = as.matrix(apply(luDF, 2, as.numeric))#Equivalent.
  
  #Search for the matching values, narrowing rows down by ID:
  if (shaded)
  {
    rowRange = which(shadeID == "Shaded")
  }
  else
  {
    rowRange = which(shadeID == "Unshaded")
    
    if (slopePct <= 30)
    {
      matches = which(slopeClass[rowRange] == "0-30%")
      rowRange = rowRange[matches]
    }
    else
    {
      matches = which(slopeClass[rowRange] == ">30%")
      rowRange = rowRange[matches]
    }
  }
  
  matches = which(aspects[rowRange] == aspectCardinal)
  theRow = rowRange[matches]
  
  #We will now be down to one row:
  if (length(theRow) != 1)
  {
    stop("Search does not narrow to a single row.")
  }
  
  #Find the matching column:
  steps = seq(1, length(hourStart), by = 3)#The hours ranges are repeated by threes (B, L, A).
  for (j in steps)
  {
    if (hourOfDay >= hourStart[j] && hourOfDay < hourEnd[j])
    {
      for (k in j:(j+2))
      {
        if (elevation == elevationCode[k])
        {
          theCol = k
          break
        }
      }
    }
  }
  
  return(luDF[theRow, theCol])
}

#Make an estimate of the 100-hour dead fuel moisture:
#The estimation is based on the 1-hour Fosberg fuel moisture.  Rothermel 1981 suggests adding 1%.
#NWCG suggests adding ~1-2%.  We use the middle of the latter recommendation.
#
#Parameters:
#oneHrFM = Fosberg NWCG prediction of 1-hr fuel moisture (fraction: water weight/dry fuel weight).
#
#Returns: 10-hour fuel moisture (fraction: water weight/dry fuel weight).
NWCG_10hrFM <- function(oneHrFM)
{
  return(oneHrFM * 1.015)#Add 1.5%.
}

#Make an estimate of the 100-hour dead fuel moisture:
#The estimation is based on the 1-hour Fosberg fuel moisture.  Rothermel 1981 suggests adding 2%.
#NWCG suggests adding ~2-4%.  We use the middle of the latter recommendation.
#
#Parameters:
#oneHrFM = Fosberg NWCG prediction of 1-hr fuel moisture (fraction: water weight/dry fuel weight).
#
#Returns: 100-hour fuel moisture (fraction: water weight/dry fuel weight).
NWCG_100hrFM <- function(oneHrFM)
{
  return(oneHrFM * 1.03)#Add 3%.
}

#Unit Tests:----------------------------------------------------------------------------------------

#Test FosbergNWCG_1HrFM() using the set of daytime scenarios outlined on Rothermel 1983 pg. 20:
#This is not an ideal unit test since it requires external file specificaitons.
#
#Parameters:
#tableA_Path - tableD_Path = Paths to tab delimited files holding the lookup table data.
#  These tables are used in Rothermel 1981 and are also the NWCG website and the 2022 IPTG.
#
#Returns: True if the test passes, false otherwise.
#
#Modified from Proj_11_Exp_20_Analysis.r Rothermel1983_FMTest2().
FosbergNWCG_1HrFM_UnitTest <- function(tableA_Path, tableB_Path, tableC_Path, tableD_Path)
{
  pass = TRUE
  #print("Daytime fuel moisture scenarios (Rothermel 1983, pg. 20):")#Verbose?????
  
  #Should be 6%:
  valA = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
                           temp = 80, rh = 20, monthOfYear = 8, hourOfDay = 13, slopePct = 35,
                           aspectCardinal = "N")
  if (valA != 6)
  {
    print(paste("a. Expected value = 6%, calculated: ", valA, "%", sep = ""))
    pass = FALSE
  }
  
  #Should be 4%:
  valB = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
                           temp = 80, rh = 20, monthOfYear = 8, hourOfDay = 13, slopePct = 35,
                           aspectCardinal = "E", elevation = "B")
  if (valB != 4)
  {
    print(paste("b. Expected value = 4%, calculated: ", valB, "%", sep = ""))
    pass = FALSE
  }
  
  #Should be 7%:
  valC = FosbergNWCG_1HrFM(tableA_Path, tableB_Path, tableC_Path, tableD_Path,
                           temp = 80, rh = 20, monthOfYear = 8, hourOfDay = 13, slopePct = 35,
                           aspectCardinal = "N", shaded = TRUE)
  if (valC != 7)
  {
    print(paste("c. Expected value = 7%, calculated: ", valC, "%", sep = ""))
    pass = FALSE
  }
  
  return(pass)
}
