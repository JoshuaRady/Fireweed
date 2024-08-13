#FireweedUtils.r
#Programmed by: Joshua M. Rady
#Woodwell Climate Research Center
#2024
#Reference: Proj. 11 Exp. 19
#
#Description:---------------------------------------------------------------------------------------
#  This file is part of the Fireweed wildfire code library.
#  This file contains (general) utilities.
#
#___________________________________________________________________________________________________

#This utility checks that the parameters (vectors) passed have the same length.  Between 2 and 4
#arguments are accepted.
#Adapted from code originally in in CalcWeightings().
#
#Alternatively we could return the length if true and otherwise -1, but this seems likely to create
#more work in practice.  In C positive = TRUE, 0 = FALSE, which could be more useful.
#In the current usage we expect the arguments to be a mix of numeric and logical vectors.  We
#could add checking for this.
SameLengths <- function(arg1, arg2, arg3 = NULL, arg4 = NULL)
{
  #Put the argments in a list removing NULL elements:
  argList = list(arg1, arg2, arg3, arg4)
  argList = argList[!sapply(argList, is.null)]
  #This would work too for omitted arguments but would ignore any zero length vectors passed in:
  #argList = argList[length(argList) != 0]
  
  #Are arguments the same length?
  if (all(sapply(argList, length) == length(arg1)))
  {
    #return(length(arg1))
    return(TRUE)
  }
  else
  {
    #return(-1)
    return(FALSE)
  }
}

#This utility checks that a value falls in a valid range.
#
#R array aware.
InRange <- function(value, low, high)
{
  if (any(value < low) || any(value > high))
  {
    return(FALSE)
  }
  else
  {
    return(TRUE)
  }
}

#Check if a value is from 0 to 1, a valid proportion.
#
#R array aware.
ValidProportion <- function(value)
{
  return(InRange(value, low = 0, high = 1))
}
