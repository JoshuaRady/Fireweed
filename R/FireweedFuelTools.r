#FireweedFuelTools.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 15
#
#Description:---------------------------------------------------------------------------------------
#  This is part of the Fireweed wildfire code library.
#  This file contains functions for manipulating the fuel amounts and converting them to fuel model
#loadings.
#
#___________________________________________________________________________________________________

source("FireweedUtils.r")

#Code:----------------------------------------------------------------------------------------------

#Take a set of fuel loadings and redistribute the fuel to a second set of size classes.  All the
#fuel for each input class is placed into the nearest output size class.
#
#Parameters:
#inputSizes = A vector of fuel sizes (SAVs, diameters, etc.) for the input fuel loadings.
#loadings = A vector of fuel loading masses for each size class in inputSizes.
#outputSizes = The size classes to redistribute to.  Units much match inputSizes.
#
#Returns: A vector of fuel loadings for each size class specified with outputSizes.  Units will
#match loadings.
RedistributeFuelNearest <- function(inputSizes, loadings, outputSizes)
{
  #Parameter checking:
  if (!SameLengths(inputSizes, loadings))
  {
    stop("distribSAVs and distribWts must have equal lengths.")
  }
  
  #If the size classes are the same we are done:
  if (length(inputSizes) == length(outputSizes) && all(inputSizes == outputSizes))
  {
    w_o_final = loadings#w_o_initial
  }
  else#Otherwise redistribute the mass to the requested size classes:
  {
    w_o_final = vector("numeric", length(outputSizes))
    
    #Find the closest output size class for each input class:
    for (i in 1:length(loadings))
    {
      #Need to convert to loop in C++:
      diffs = abs(inputSizes[i] - outputSizes)
      closest = which.min(diffs)
      w_o_final[closest] = w_o_final[closest] + loadings[i]
    }
  }
  
  return(w_o_final)
}

#Take a set of fuel loadings and redistribute the fuel to a second set of size classes.  When a
#fuel size class from the original falls between two of the output set it's fuel is split
#proportionally based on proximity to its neighbors.
#
#This make sense when there are about one or fewer input classes fall between sets of output
#classes.  If multiple classes fall between two output classes some fuel from the smaller input
#class will be contributed to an output size class larger than it adjacent input class.  This isn't
#logical.  In that case the nearest bin approach will be better.
#
#Parameters:
#inputSizes = A vector of fuel sizes (SAVs, diameters, etc.) for the input fuel loadings.
#loadings = A vector of fuel loading masses for each size class in inputSizes.
#outputSizes = The size classes to redistribute to.  Units much match inputSizes.
#
#Returns: A vector of fuel loadings for each size class specified with outputSizes.  Units will
#match inputSizes.
RedistributeFuelProportional <- function(inputSizes, loadings, outputSizes)
{
  #Parameter checking:
  if (!SameLengths(inputSizes, loadings))
  {
    stop("distribSAVs and distribWts must have equal lengths.")
  }
  
  #The output size classes need to be in order for the algorithm below to work properly.  As a
  #matter of convention we expect SAVs to be in descending order and diameters to be in ascending
  #order. The following code was written primarily for SAV in descending order.  In the case of
  #ascending sizes we reverse the size bins and then reverse the output before returning it below.
  sizesAscending = FALSE
  
  if (!is.unsorted(outputSizes))#Sizes are in ascending order.
  {
    outputSizes = rev(outputSizes)#In place.
    sizesAscending = TRUE
  }
  else if (is.unsorted(rev(outputSizes)))
  {
    stop("outputSizes must be in ascending or descending order.")
  }
  
  #It is not strictly necessary for the input sizes to be ordered to be sorted.  They do need to be
  #in the same order as the output for comparison purposes though:
  if (sizesAscending)
  {
    inputSizesSorted = sort(inputSizes)
  }
  else
  {
    inputSizesSorted = sort(inputSizes, decreasing = TRUE)
  }
  
  #If the size classes are the same we are done:
  if (length(inputSizes) == length(outputSizes) && all(inputSizesSorted == outputSizes))
  {
    w_o_output = loadings
  }
  else#Otherwise redistribute the mass to the requested size classes:
  {
    w_o_output = vector("numeric", length(outputSizes))
    
    #Redistribute the mass for each initial size class to the proper recipient classes:
    for (i in 1:length(loadings))
    {
      for (j in 1:length(outputSizes))
      {
        if (inputSizes[i] > outputSizes[j])#Greater = to the left of.
        {
          if (j == 1)#The SAV is greater than any in the output.  Put it into the the largest class:
          {
            w_o_output[j] = w_o_output[j] + loadings[i]
            break
          }
          else
          {
            #Calculate the distance to the adjacent sizes and use to split between them:
            range = outputSizes[j - 1] - outputSizes[j]
            diffLeft = outputSizes[j - 1] - inputSizes[i]
            diffRight = inputSizes[i] - outputSizes[j]
            
            w_o_output[j - 1] = w_o_output[j - 1] + loadings[i] * (diffRight / range)
            w_o_output[j] = w_o_output[j] + loadings[i] * (diffLeft / range)
            break
          }
        }
        else if (inputSizes[i] == outputSizes[j])#If it is a perfect match put it all in this class.
        {
          w_o_output[j] = w_o_output[j] + loadings[i]
          break
        }
        else#if (inputSizes[i] < outputSizes[j])
        {
          #If the fuel is smaller than any in the output put it into the the smallest size class:
          if (j == length(outputSizes))
          {
            w_o_output[j] = w_o_output[j] + loadings[i]
            break
          }
          #Otherwise continue to the next size bin.
        }
      }
    }
  }
  
  if (sizesAscending)
  {
    w_o_output = rev(w_o_output)
  }
  
  return(w_o_output)#We could name the values with sizes.
}

#Distribute a pool of fuel across a set of size classes.
#
#  This function takes an estimated distribution of fuel sizes for a fuel bed.  The distribution is
#provided as a discreet list of fuel sizes and weighing factors.  The fuel mass input is
#distributed according to this distribution and is then redistributed to the set of requested
#output size classes.
#
#  The purpose of this function is to make it easy to convert a total dead fuel estimate to a fuel
#model loading using any fuel distribution estimate from the literature.  Fuel bed descriptions can
#take different forms and may not match the size bins for the desired fuel model.  The function
#handles this conversion in the process.  The function expects weighting factors (fractions) but
#for convenience fuel loadings may be input and they will be converted to weights.
#
#Parameters:
#distribSizes = A vector of fuel sizes (SAVs, diameters, etc.) describing the fuel distribution.
#distribWts = A vector of weightings (fractions summing to 1) for each size class in distribSizes.
#  Alternatively typical fuel loadings (i.e lb/ft^2 | kg/m^2) for each size class in distribSizes.
#totalLoading = The total fuel mass to distribute. (lb/ft^2 | kg/m^2)
#outputSizes = The SAV size classes to redistribute to.
#method = The method for sorting: "Nearest" or "Proportional".
#
#Returns: A vector of fuel loadings for each size class in outputSizes.  Units will match
#totalLoading.
DistributeFuel <- function(distribSizes, distribWts, totalLoading, outputSizes,
                           method = "Proportional")
{
  #Parameter checking:
  if (!SameLengths(distribSizes, distribWts))
  {
    stop("distribSizes and distribWts must have equal lengths.")
  }
  #Further checking will occur in the called fuctions.
  
  #Make sure the weights are valid and add up to 1 (allowing for floating point error):
  if (!isTRUE(all.equal(sum(distribWts), 1)))
  {
    #If they aren't we convert them:
    warning("Adjusting weights to sum to zero.")
    distribWts = distribWts / sum(distribWts)
  }
  
  #Distribute the litter mass to the provided distribution. Multiply the mass by the weights:
  w_o_initial = totalLoading * distribWts
  
  #Redistribute the mass to the requested size classes:
  if (method == "Nearest")
  {
    w_o_final = RedistributeFuelNearest(distribSizes, w_o_initial, outputSizes)
  }
  else if (method == "Proportional")
  {
    w_o_final = RedistributeFuelProportional(distribSizes, w_o_initial, outputSizes)
  }
  else
  {
    stop("Unknown method.")
  }
  
  return(w_o_final)
}
