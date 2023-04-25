#FwRothermelAlbiniFireSpread.r
#Programed by: Joshua M. Rady
#Woodwell Climate Research Center
#Started: 2/27/2023
#
#Description:---------------------------------------------------------------------------------------
# This file is part of the Fireweed fire code library.  It contains R implementations of the
#Rothermel fire spread model (Rothermel 1972) with the modifications of Albini (Albini 1972).
#
# The Rothermel/Albini equations have been implemented in a form as close to the originals as was
#possible.  The modifications of Albini 1972 are considered official parts of the Rothermel fire
#spread model and all equations include these modifications where they apply.  The Rothermel model
#consists of a nested set of equations. The library presents functions for each equation with useful
#output.  Some functions are very simple and could be inlined but a focus on clarity and modularity
#have been prioritized over compactness.
#
# This code was started in Proj_11_Exp_3_Analysis_D5.r from Project 11 Experiment 3.
#
#References:----------------------------------------------------------------------------------------
#Richard C. Rothermel.
#A mathematical model for predicting fire spread in wildland fuels.
#Res. Pap. INT-115. Ogden, UT: U.S. Department of Agriculture, Intermountain Forest and Range
#Experiment Station. 40 p. 1972
#
#Computer-based models of wildland fire behavior: a user's manual.
#Albini, F. A.
#Intermountain Forest and Range Experiment Station, Forest Service, U.S. Department of Agriculture.
#1976.
#
#Andrews, Patricia L., Miguel G. Cruz, and Richard C. Rothermel.
#Examination of the wind speed limit function in the Rothermel surface fire spread model.
#International Journal of Wildland Fire 22(7): 959-69, 2013. http://dx.doi.org/10.1071/WF12122
# This review summarizes and contextualizes the equations of Rothermel and Albini as well as
#related work and was an important reference for preparing this code.
#
#Notation:------------------------------------------------------------------------------------------
# The equations contain mathematical notation and variables with characters that cannot be directly
#represented in R/C++.  To represent the variables the following translations were used.
#
#Greek letters in variable names:
# Many variables use Greek characters.  In most cases these are translated using their English names
#with the case indicating the case of the character, e.g. β -> Beta and σ -> sigma.  In a few cases
#Greek variable names have been changed to abbreviations or descriptive names.  Greek is used in
#the comments.
#[See table of variables below.]
#
#Subscripts:
# Variables with subscripts are represented with underscores, e.g. Ab (A sub b) -> A_b.  A number of
#variables have two levels of subscript, the second representing fuel type indexes (i and j, see
#below).  These are represented with underscores as well, e.g (Ab)ij ((A sub b) sub ij) -> A_b_ij.
#
#[Add variables table...]
#
#Variable notes:
# - The surface-area-to-volume ratio for fuels is notated as σ (sigma) and abbreviated as sav or
#SAV in different places in the papers.  We use SAV in the code.
# - It is unclear if fuel loading is w0 or wo.  In Rothermel 1972 it is not clear and in reprints it
#varies.  We usse w sub o (w_0).
# - Total mineral content is sometimes notated S sub t and sometimes S sub T.  We use S sub t (S_t).
#
#Units:
# The original equations used English units.  Units are indicated for function inputs and outputs.
#SI units will be added later.
#___________________________________________________________________________________________________

#Globals:-------------------------------------------------------------------------------------------

#Code:----------------------------------------------------------------------------------------------

#Rothermel & Albini Spread Model:-------------------------------------------------------------------

#Propagating Flux Ratio:
#The propagating flux ratio, represented as lower case xi, is the proportion of the reaction
#intensity that heats fuels adjacent to the fire front.
#
#The equation is the same for Rothermel 1972 (eq. 42/76) and Albini 1976 (pg?????):
#ξ = (192 + 0.2595σ)^-1 exp[(0.792 + 0.681σ^0.5)(β + 0.1)]
#
#Input variables / parameters:
#β = (mean) packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
#  For heterogeneous fuels the mean packing ratio is passed in.
#σ (sav) = characteristic surface-area-to-volume ratio (ft^2/ft^3)
#  For heterogeneous fuels the fuel bed level SAV is used.
#
#Output units: Dimensionless proportion
PropagatingFluxRatio <- function(packingRatio, sav)
{
  xi = (192 + 0.2595 * sav)^-1 * exp((0.792 + 0.681 * sav^0.5) * (packingRatio + 0.1))
  
  #Convert units????
  
  return(xi)
}

#Packing Ratio:
# The packing ratio is the fraction of the (surface) fuel bed volume occupied by fuel, aka
#compactness.
#
#Rothermel 1972 equation 31:
#β = ρb/ρp
#
#Input variables / parameters:
#ρb (rho sub b) = fuel array bulk density (density, originally lb/ft^3)
#ρp (rho sub p) = fuel particle density (density, originally lb/ft^3)
#  For standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
#
#Output units: Dimensionless ratio
PackingRatio <- function(fuelArrayBulkDensity, fuelParticleDensity)
{
  return(fuelArrayBulkDensity / fuelParticleDensity)
}

#Mean Packing Ratio:
#For heterogeneous fuelbeds the mean packing ratio must be calculated.
#
#Rothermel 1972 equation 74:
#β = 1/δ Σi Σj (wo)ij/(ρp)ij
#For i = 1 to m fuel categories (live vs. dead) and j = 1 to n fuel size classes.
#The original notation includes from and to sum subscripts and bars over beta, rho, and w.
#
#Input variables / parameters:
#(w sub o [bar]) = An array of fuel loadings for each fuel category (i: live vs. dead) and fuel size
#class (j) (lb/ft^2).
#w_oij?????
#ρp (rho sub p [bar]) = particle densities for each fuel class (lb/ft^3)
#  For standard fuel models particle density is 32 lb/ft^3. (30-46 in some others.)
#δ (delta) = fuel bed depth (ft)
#
#Output units: Dimensionless ratio
MeanPackingRatio <- function(fuelLoadings, fuelBedDepth, particleDensities = 32)#Order for defaults?
{
  #Parameter checking:
  numLoadings = length(fuelLoadings)
  numDensities = length(particleDensities)
  
  #If only one particle density is provided assume that is it the same for all fuel classes:
  if (numDensities == 1)
  {
    particleDensities = array(data = particleDensities, dim = numLoadings)
  }
  else#Otherwise one should be provided for each fuel class.
  {
    if (numDensities != numLoadings)
    {
      stop("The number of fuel loadings and particle densities do not match.")
    }
  }
  
  #Calculate w_o / rho_p for each fuel element:
  x = fuelLoadings / particleDensities
  
  beta_bar = sum(x) / fuelBedDepth#Mean packing ratio
  return(beta_bar)
}

#Optimum Packing Ratio:
# This is the packing ratio (compactness) at which the maximum reaction intensity will occur.
#
#Rothermel 1972 equations 37,69:
#βop = 3.348(σ)^-0.8189
#
#Input variables / parameters:
#σ (sav) = characteristic surface-area-to-volume ratio (ft^2/ft^3)
#  For heterogeneous fuels the SAV of the fuel bed / complex is used.
#
#Output units: Dimensionless ratio
OptimumPackingRatio <- function(sav)
{
  optPackingRatio = 3.348 * (sav)^-0.8189
  return(optPackingRatio)
}

#Weighting Factors:-----
#The spread model uses weights in many of the calculations for heterogeneous fuels.  These need to
#be calculated...
#
#Input variables / parameters:
#σ (SAVs) = An array of characteristic surface-area-to-volume ratios for the fuel classes
#  (ft^2/ft^3).
#w_oij (w sub o ij) = An array of fuel loadings for each fuel category (i: live vs. dead) and fuel
#  size class (j) (lb/ft^2).
#     Rotethermel used w_o while Albini used w_oij.
#ρp (rho sub p) = ...
#liveDead = An array indicating if each index in each of the other input variables represents a
#  dead (1) or live (2) or fuel categories.
#Need to know the live dead indexes!!!!!
#
#Output units: ...
#
#Note: It makes sense to calculate these togeather and they only need to be calculated onces for a
#give spread rate scenario.  However, I'm not sure the best way the handle the outputs.  Would it
#be better to put them in a global?
#Calc_f_ij_Weights
CalcWeightings <- function(SAV_ij, w_os, rho_ps,#fuelParticleDensity = 32,#Change name?
                           liveDead)
{
  #Add error checking...
  
  numFuelTypes = length(SAV_ij)#Types = sum of size classes in both categories.
  
  #Calcualte the (mean) total surface area for each fuel component:
  #Rothermel equation 53?:
  #Aij = (σ)ij (w0)ij ⁄(ρp)ij
  A_ij = SAV_ij * w_os / rho_ps
  
  #Mean total surface area by live / dead fuel categories:
  #Rothermel equation 54?:
  #Ai = ΣjAij
  A_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    A_i[liveDead[k]] = A_i[liveDead[k]] + A_ij[k]
  }
  
  #Mean total surface area of the fuel:
  #Rothermel equation 55?:
  #AT = ΣiAi
  A_T = sum(A_i)#Single scalar value.
  
  #f_ij weighting factor...
  #Rothermel equation 56?:
  #fij = Aij/Ai
  
  #f_ij = array(dim = numFuelTypes)
  f_ij = array(data = 0, dim = numFuelTypes)#A vector seems nicer.
  #f_ij = vector(mode = "numeric", length = numFuelTypes)
  for (l in 1:numFuelTypes)
  {
    if (A_i[liveDead[l]] != 0)
    {
      f_ij[l] = A_ij[l] / A_i[liveDead[l]]
    }
    #A_i can be 0 if there in no fuel in either the live or dead fuel category.  If A_i = 0 then
    #A_ij for this fuel should also be 0.  In this case we avoid a divide by 0 and give an
    #appropriate weight of 0.  Given the initialization we don't have to do this explicitly.
  }
  
  #fi
  #Rothermel equation 57?:
  #fi = Ai/AT
  f_i = A_i / A_T#Array of 2.
  
  #g_ij:
  #The final set of weights was added in Albini 1976 to get arround a logical problem of using f_ij
  #for fuel loading.
  #Albini 1972 pg. 15:
  #...
  #The notation here is a bit confusing in my opinion.
  
  #What size subclass is each fuel type in?
  #A better names is needed.  This is the subclass that each fuel type is in: subclassOfFuelType?
  subclass = array(data = 0, dim = numFuelTypes)
  #subclassWt = array(data = 0, dim = 6)#Didn't end up needed in this since f_ij is used.
  subclassTotal = array(data = 0, dim = 6)
  for (n in 1:numFuelTypes)
  {
    if (SAV_ij[n] >= 1200)
    {
      subclass[n] = 1
      #subclassWt[1] = subclassWt[1] + f_ij[n]
    }
    else if (SAV_ij[n] >= 192)
    {
      subclass[n] = 2
      #subclassWt[2] = subclassWt[2] + f_ij[n]
    }
    else if (SAV_ij[n] >= 192)
    {
      subclass[n] = 3
      #subclassWt[3] = subclassWt[3] + f_ij[n]
    }
    else if (SAV_ij[n] >= 96)
    {
      subclass[n] = 4
      #subclassWt[4] = subclassWt[4] + f_ij[n]
    }
    else if (SAV_ij[n] >= 48)
    {
      subclass[n] = 5
      #subclassWt[5] = subclassWt[5] + f_ij[n]
    }
    else if (SAV_ij[n] >= 16)
    {
      subclass[n] = 6
      #subclassWt[6] = subclassWt[6] + f_ij[n]
    }
    #Add the 
    #subclassWt[subclass[n]] = subclassWt[subclass[n]] + f_ij[n]
    #Bin the fuels type weights by subclass they belong to... 
    subclassTotal[subclass[n]] = subclassTotal[subclass[n]] + f_ij[n]
  }
  
  #Complete the weight calculation:
  #Determine the fraction of surface area in each size subclass...
  #subclassWt = subclassWt / A_T
  #g_ij = subclassTotal / A_T#Not quite...
  
  #subclassWt = subclassTotal / A_T#Length 6.
  #Hold on how did A_T get in here in the first place?
  
  #Assign the subclass weights to each size class.  Some may share the same weight.
  g_ij = vector(mode = "numeric", length = numFuelTypes)#Or array?
  for (k in 1:numFuelTypes)#k reused...
  {
    #g_ij[k] = subclassWt[subclass[k]]
    #g_ij[k] = subclassTotal[subclass[k]]
    
    if (subclass[k] != 0)
    {
      g_ij[k] = subclassTotal[subclass[k]]
    }
    #If a fuel class is not fully specified, i.e. has an invalid SAV of 0, it will not be mapped to
    #a size subclass.  In that case leave g_ij[k] = 0.
    #A value of NA might be more logical but a 0 wieght makes the math simpler.
  }
  
  return(list(f_ij = f_ij, f_i = f_i, g_ij = g_ij))
  #Or save in globals?
}


#Fuel Bed Surface-area-to-volume Ratio:
#  For heterogeneous fuels a SAV for the entire fuel bed must be calculated.
#
#Input variables / parameters:
#σ (SAV_ij) = An array of characteristic surface-area-to-volume ratios for the fuel classes
#(ft^2/ft^3).
#f_ij (f sub ij) = Weighting factors for each fuel type (dimensionless).
#f_i (f sub i) = Weighting factors for each fuel live/dead category (dimensionless).
#liveDead = An array indicating if each index in each of the other input variables represents a
#  dead (1) or live (2) or fuel categories.
#
#Output units: ft^2/ft^3
#FuelBedSAV <- function(SAVs, weights, liveDead)#f_ij, SAVs -> SAV_ij?
FuelBedSAV <- function(SAV_ij, f_ij, f_i, liveDead)#f_ij, SAVs -> SAV_ij?
{
  #Add error checking...
  
  #numFuelTypes = length(SAVs)
  numFuelTypes = length(SAV_ij)#Types = sum of size classes in both categories.
  
  #Characteristic live and dead SAVs:
  #Rothermel 1972 equation 72:
  #σi = Σj fijσij (~ over sigma sub i and bar over sigma sub ij in original)
  
  sav_i = c(0,0)#Or sigma_i?
  for (k in 1:numFuelTypes)
  {
    #sigma_i[liveDead[k]] = sigma_i[liveDead[k]] + SAVs[k]
    #sav_i[liveDead[k]] = sav_i[liveDead[k]] + SAV_ij[k]
    #Left out weights!
    
    sav_i[liveDead[k]] = sav_i[liveDead[k]] + (f_ij[k] * SAV_ij[k])
  }
  
  #Sum the live and dead components to get the final value...
  #Rothermel 1972 equation 71:
  #σ = Σi fiσi (~ over sigma and sigma sub i in original)
  fuelBedSAV = sum(f_i * sav_i)
  
  return(fuelBedSAV)
}

#Optimum (Potential) Reaction Velocity:
#  This is a measure of the optimum (potential) fuel consumption rate (fire efficiency / reaction
#time).  The 'optimum' rate is the ideal rate that would occur for alpha cellulose in the absence
#of minerals and moisture.
#Notation: Γ' (Gamma prime)
#
#  This version includes the Albini modification...
#
#Input variables / parameters:
#β = (mean) packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
#  For heterogeneous fuels the mean packing ratio is passed in.  (see PackingRatio()?????
#σ (sav) = characteristic surface-area-to-volume ratio (ft^2/ft^3)
#  For heterogeneous fuels the SAV of the fuel bed / complex is used.
#
#Output units: min^1
OptimumReactionVelocity <- function(meanPackingRatio, sav)#, liveDead)
{
  #Calculate the maximum reaction velocity (min^-1):
  #This is rate for moisture free fuel with mineral compostion of alpha cellulose.
  #Rothermel 1972 equations 36,68:
  #Γ'max = σ^1.5/(495 + 0.0594σ^1.5)
  GammaPrimeMax = sav^1.5 / (495 + 0.0594 * sav^1.5)
  
  optPackingRatio = OptimumPackingRatio(sav)
  
  #"Arbitrary" variable (no units?????):
  #Albini 1976 pg. 15????
  #A = 133σ^-0.7913
  A = 133 * sav^-0.7913
  
  #These are combined to produce the optimal reaction velocity (min^-1):
  #Rothermel 1972 equation 38:
  #Γ' = Γ'max(β/βop)^A exp[A(1 - β/βop)]
  GammaPrime = GammaPrimeMax * (meanPackingRatio/optPackingRatio)^A  *
    exp(A * (1 - meanPackingRatio/optPackingRatio))
  
  return(GammaPrime)
}

#Net Fuel Load:
# The fuel load (w sub n) is mass per ground area of dry combustible fuel removing the
#noncombustible mineral mass.
#
#Albini 1976 pg. 14?????
#wn = wo(1 - St)
#Or alt notation? w_n = w_0(1 - S_t)
#
#Input variables / parameters:
#w_o (w sub o) = Oven dry fuel load (lb/ft^2).  This includes combustible and mineral fraction.
#S_t (S sub t) = Total mineral content (fuel particle property: mineral mass / total dry mass,
#unitless fraction).  For all standard fuel models this is 5.55% (0.0555).
#
#Output units: lb/ft^2
#For heterogeneous fuel beds this is computed for each fuel catagory?
NetFuelLoad_Albini <- function(w_o, S_t)
{
  w_n = w_o * (1 - S_t)
  
  #Convert units????
  
  return(w_n)
}

#For hetogenous fuels the net fuel load for each fuel catagory (live/dead) is calculated using
#weights.
#
#Input variables / parameters:
#(w0)ij ((w sub o) sub ij) = Array of oven dry fuel load for each fuel class (lb/ft^2).
#(ST)ij ((S sub T) sub ij) = Array of total mineral content for each fuel class (unitless fraction).
#g_ij = Net fuel load weights (Albini 1976).
#
#Returns: w_n_i ((w sub o) sub i) = An 
#Output units: lb/ft^2
NetFuelLoad_Albini_Het <- function(w_o_ij, S_T_ij, g_ij, liveDead)#Name?????
{
  #Argument checking...
  
  numFuelTypes = length(w_o_ij)
  
  #Calculate the net fuel load for each fuel class:
  #Rothermel equation ?????? modified by Albini 1976 pg. 14:
  #(wn)ij = (w0)ij (1 – (ST)ij)
  w_n_ij = w_o_ij * (1 - S_T_ij)
  
  #Sum the net fuel load for each 
  #Albini 1976 pg. 15:
  #(wn)i = Σjgij(wn)ij
  w_n_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    w_n_i[liveDead[k]] = w_n_i[liveDead[k]] + (g_ij[k] * w_n_ij[k])
  }
  
  return(w_n_i)
}

#Damping Coefficients:-----

#Moisture Damping Coefficient:
# This returns the extent to which fuel moisture reduces combustion for one fuel component.
#
#Input variables / parameters:
#Mf (M sub f) = Fuel moisture content (fraction, water weight/dry fuel weight)
#Mx (M sub x) = Moisture of extinction (fraction, water weight/dry fuel weight)
#
#Output units: Dimensionless coefficient
MoistureDampingCoefficient <- function(M_f, M_x)
{
  #Calculate the ratio of fuel moisture content to moisture of extinction:
  #Rothermel 1972 equation 29?,65:
  #rM = Mf/Mx (max = 1.0)
  r_M = M_f/M_x
  
  if (r_M > 1.0)
  {
    r_M = 1.0
  }
  
  #Use the ratio to calculate the damping coefficient:
  #Rothermel 1972 equations 29,64:
  #ηM = 1 - 2.59rM + 5.11(rM)^2 - 3.52(rM)^3
  eta_M = 1 - 2.59 * r_M + 5.11 * (r_M)^2 - 3.52 * (r_M)^3
  
  return(eta_M)
}

#Moisture Damping Coefficient (heterogeneous fuels):
#  For heterogenous fuel beds the moisture damping coefficient is calculated for each fuel category
#(live/dead).
#
#Input variables / parameters:
#(Mf)ij ((M sub f) sub ij) = Fuel moisture content for for each fuel type (fraction, water
#  weight/dry fuel weight)
#(Mx)i ((M sub x) sub i = Moisture of extinction each fuel category (fraction, water weight/dry
#  fuel weight).
#f_i (f sub i) = Weighting factors for each fuel live/dead category.
#liveDead = An array indicating if each index in each of the other input variables represents a
#  dead (1) or live (2) or fuel categories.
#
#Output units: Dimensionless coefficient (array length 2)
MoistureDampingCoefficient_Het <- function(M_f_ij, M_x_i, f_ij, liveDead)
{
  numFuelTypes = length(M_f_ij)
  
  #Calculate the weighted moisture content for each fuel category:
  #Rothermel 1972 equations 66:
  #(Mf)i = Σj fij (Mf)ij
  M_f_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    M_f_i[liveDead[k]] = M_f_i[liveDead[k]] + (f_ij[k] * M_f_ij[k])
  }
  
  #Calculate the moisture damping coefficient for each fuel catagory:
  #Rothermel 1972 equations 64,65:
  #(ηM)i = 1 – 2.59(rM)i + 5.11(rM)i2 – 3.52(rM)i3 (max = 1)
  eta_m_i = c(0,0)
  eta_m_i[1] = MoistureDampingCoefficient(M_f_i[1], M_x_i[1])
  eta_m_i[2] = MoistureDampingCoefficient(M_f_i[2], M_x_i[2])
  
  return(eta_m_i)
}

#Live Fuel Moisture of Extinction:
#  The live fuel moisture of extinction is calculated from the dead fuel moisture of extinction and
#in relation to the ratio of live to dead fuels.
#MORE!!!!!
#
#Input variables / parameters:
#(Mf)ij ((M sub f) sub ij) = Fuel moisture content for for each fuel type (fraction, water
#  weight/dry fuel weight)
#(Mx)1 ((M sub x) sub 1 = Dead fuel moisture of extinction (fraction, water weight/dry fuel weight).
#w_o_ij ...
#(w0)ij ((w sub o) sub ij) = Array of oven dry fuel load for each fuel class (lb/ft^2).
#σij (SAV_ij) = An array of characteristic surface-area-to-volume ratios for the fuel classes.
#liveDead = An array indicating if each index in each of the other input variables represents a
#  dead (1) or live (2) or fuel categories.
#
#Output units: fraction, water weight/dry fuel weight (M_x_2 / M_x_living)
#LiveFuelMx <- function(M_f_ij, M_x_1, w_o_ij, SAV_ij, liveDead)
LiveFuelMoistureOfExtinction <- function(M_f_ij, M_x_1, w_o_ij, SAV_ij, liveDead)
{
  #Calculate dead:live loading ratio, notated W:
  #Albini 1976 pg. 16:
  #W = Σj(w0)1jexp(-138/σ1j) / Σj(w0)2jexp(-500/σ2j)
  #
  #These sums could be done in a vector aware way but that we limit ourselves to C compatible code.
  #The loops could also be combined with a live/dead conditional.
  liveSum = 0
  for (k in which(liveDead == 1))
  {
    liveSum = liveSum + w_o_ij[k] * exp(-138 / SAV_ij[k])
  }
  
  deadSum = 0
  for (k in which(liveDead == 2))
  {
    deadSum = deadSum + w_o_ij[k] * exp(-500 / SAV_ij[k])
  }
  
  W = liveSum / deadSum
  
  #Calculate fine dead fuel moisture as:
  #Albini 1976 pg. 16:
  #Mf,dead = Σj(Mf)1j(w0)1jexp(–138/σ1j) / Σj(w0)1jexp(–138/σ1j)
  top = 0#Better names
  bottom = 0
  for (k in which(liveDead == 1))
  {
    common = w_o_ij[k] * exp(-138 / SAV_ij[k])
    top = top + M_f_ij[k] * common
    bottom = bottom + common
  }
  
  M_f_dead = top / bottom
  
  #Calcualte the live fuel moisture of extinction ((Mx)2):
  #Rothermel 1972 equation 88 with Albini 1976 pg. 16 modifications:
  #(Mx)2 = 2.9W[1 – Mf,dead⁄(Mx)1] - 0.226, (min = (Mx)1)
  M_x_2 = 2.9 * W * (1- M_f_dead / M_x_1) - 0.226
  
  if (M_x_2 < M_x_1)
  {
    M_x_2 = M_x_1
  }
  
  return(M_x_2)
}

#Mineral Damping Coefficient:
#
#Albini 1976 pg. 14 adds a maximum to Rothermel 1972 equation 30/62:
#𝜂s = 0.174Se^-0.19 (max = 1.0)
#
#Input variables / parameters:
#Se (S sub e) = effective mineral content (fuel particle property: (mass minerals – mass silica) /
#  total dry mass, unitless fraction).  For all standard fuel models this is 1% (0.01).
#
#Output units: Dimensionless coefficient
MineralDampingCoefficient <- function(S_e)
{
  eta_s = 0.174 * S_e^-0.19
  
  if (eta_s < 1.0)
  {
    return(eta_s)
  }
  else
  {
    return(1.0)
  }
}

#Mineral Damping Coefficient: (heterogeneous fuels):
#  For heterogenous fuels the mineral damping coefficient is calculated for each fuel category
#(live/dead).
#
#Input variables / parameters:
#(Se)i ((S sub e) sub ij) = effective mineral content for each fuel category.
#f_ij (f sub ij) = Weighting factors for each fuel type (dimensionless).
#
#Output units: Dimensionless coefficient (array of 2)
MineralDampingCoefficient_Het <- function(S_e_ij, f_ij, liveDead)
{
  numFuelTypes = length(S_e_ij)#Types = sum of size classes in both categories.
  
  #Calculate the weighted effective mineral content for each fuel category:
  #(Se)i = Σj fij (Se)ij
  S_e_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    S_e_i[liveDead[k]] = S_e_i[liveDead[k]] + (f_ij[k] * S_e_ij[k])
  }
  
  #Caculate the mineral damping coefficient for each fuel catagory:
  #(ηs)i = 0.174(Se)i^–0.19 (max = 1)
  eta_s_i = c(0,0)
  eta_s_i[1] = MineralDampingCoefficient(S_e_i[1])
  eta_s_i[2] = MineralDampingCoefficient(S_e_i[2])
  
  return(eta_s_i)
}


#Reaction Velocity
# The reaction velocity is the actual fuel consumption rate at the fire front.
#
#Γ = Γ'ηM ηs
#
#This is not currently needed in an of itself.  It is imbedded in ReactionIntensityAlbini().
# ReactionVelocity <- function()
# {
#   gammaPrime = OptimumReactionVelocity(meanPackingRatio, sav)
# }


#Reaction Intensity:
# The reaction intensity (IR / I sub R) is the total energy released by the fire front in 
#Btu/ft^2-min ... Btu/ft^2/min #in all forms (radiation, conduction, and convection).
#  This it not the same as fireline intensity.
#
#IR (I sub R) / I_R
#Units Btu/ft^2/min or kW/m^2/min kW/min/m^2

#Reaction Intensity Rothermel homogeneous fuel:
#
#Rothermel equation ##?:
#IR = Γ'wnhηMηs
#or: I_R = Γ' x w_n x h x η_M x η_s
#
#Input variables / parameters:
#Γ' (Gamma prime) = optimum reaction velocity
#wn (w sub n) = net fuel load (for a single fuel component...)
#h = heat content of the fuel class (Btu/lb)
#ηM (eta sub M) = moisture damping coefficient
#ηs (eta sub s) = mineral damping coefficient
#
#Output units: Btu/ft^2/min
#
#Are these the best parameters?
ReactionIntensityRothermel <- function(GammaPrime, w_n, h, eta_M, eta_s)
{
  #I_R = optimumReactionVelocity * w_n * h * η_M * η_s
  #I_R = GammaPrime * w_n * h * η_M * η_s
  I_R = GammaPrime * w_n * h * eta_M * eta_s
  
  return(I_R)
}
#function(GammaPrime, w_n, h, Mf, Mx, Se)

#Reaction Intensity for Heterogeneous Fuels:
#
#Input variables / parameters:
#Γ' (Gamma prime) = optimum reaction velocity
#(wn)i ((w sub n) sub i) = net fuel load for live/dead fuel categories (UNITS!!!!!).
#hi (h sub i) = Heat content for live/dead fuel categories.
#(ηM)i ((eta sub M) sub i) = Moisture damping coefficient for live/dead fuel categories
#(UNITS!!!!!).
#(ηs)i ((eta sub s) sub i) = mineral damping coefficient for live/dead fuel categories (UNITS!!!!!).

#Some representation of the fuel array...
#
#Output units: Btu/ft^2/min
#ReactionIntensityAlbini <- function()
ReactionIntensity_Het <- function(GammaPrime, w_n_i, h_i, eta_M_i, eta_s_i)
{
  #Rothermel equation 58 modified by Albini 1976 pg. 17:
  #IR = Γ' Σi (wn)ihi(ηM)i(ηs)i
  IR = GammaPrime * sum(w_n_i * h_i * eta_M_i * eta_s_i)
  
  return(IR)
}

#Slope factor:
#Dimensionless multiplier that accounts for the effect of slope on spread behavior.  Same for
#homegenous and heterogeneous fuels.
#
#Rothermel 1972 equation 51,80:
#ϕs (phi sub s) = 5.275β^-0.3(tan ϕ)^2
#
#Input variables / parameters:
#β = packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
#tan ϕ = slope steepness, maximum (fraction: vertical rise / horizontal distance, unitless)
#
#Output units: Dimensionless, no metric conversion required.
SlopeFactor <- function(packingRatio, slopeSteepness)
{
  phi_s = 5.275 * packingRatio^-0.3 * slopeSteepness^2
  return(phi_s)
}

#Wind factor:
#Dimensionless multiplier that accounts for the effect of wind speed on spread behavior
#(propigating flux ratio specifically).  Same for homegenous and heterogeneous fuels.
#
#Input variables / parameters:
#σ / sav = characteristic surface-area-to-volume ratio (ft^2/ft^3) 
#β = packing ratio, the fraction of the fuel bed volume occupied by fuel (dimensionless).
#βop = optimum packing ratio (Note: optimum packing ratio is a function of SAV)
#U = wind speed at midflame height (ft/min)
#
#Output units: Dimensionless
#
#Note: It is not possible to calculate if a wind limit is indicated internal to this function
#and not all authors agree that a wind limit should be used.  U should be capped, if deemed
#appropriate prior to passing it in to this function.
WindFactor <- function(sav, packingRatio, optPackingRatio, U)
{
  #C = unnamed term
  #Rothermel 1972 equation 48,82:
  #C = 7.47exp(-0.133σ^0.55)
  C = 7.47 * exp(-0.133 * sav^0.55)
  
  #B = unnamed term
  #Rothermel 1972 equation 49,83:
  #B = 0.02526σ^0.54
  B = 0.02526 * sav^0.54
  
  #E = unnamed term
  #Rothermel 1972 equation 50,84:
  #E = 0.715exp(-3.59×10^-4 σ)
  E = 0.715 * exp(-3.59 * 10^-4 * sav)
  
  #Rothermel 1972 equation 47,79:
  #ϕw = CU^B(β/βop)^-E
  phi_w = C * U^B * (packingRatio / optPackingRatio)^-E
  
  return(phi_w)
}

#Wind limit:
#  The "wind limit" or "maximum reliable wind" is used to limit the effect of wind on the fire
#spread rate as wind speed gets high.  It caps the wind speed at a value that is a function of the
#reaction intesity.
#   There is not agreement on whether the wind limit should be used.  Albini chose to not use it,
#but his code reports if the limit was reached (Albini 1976, pg 26).  More recent work finds the
#original calculation to be flawed and presents an alternate formulation from (Andrews et. al 2013).
#However, They conclude that in general neither should be used.  They state a better alternative is
#to cap the spread rate at the “effective wind speed”.
#  We implement the original formulation as an option to be able to reproduce results that do use
#the wind limit.
#
#Input variables / parameters:
#U = wind speed at midflame height (ft/min)
#IR (I sub R, I_R) = reaction intensity (Btu/ft^2/min)
#
#Output units: adjusted wind speed (U) at midflame height (ft/min)
WindLimit <- function(U, I_R)
{
  #Rothermel 1972 Equation 87:
  if (U/I_R > 0.9)
  {
    U = 0.9 * I_R
    
    #Or Andrews et al. 2013 equation 21:
    #U = 96.8 * I_R^(1/3)
  }
  #Otherwise return U unchanged.
  
  return(U)
}

#Heat Of Preignition:
#
#Rothermel 1972 equations 12,78:
#Qig = 250 + 1,116Mf
#
#Input variables / parameters:
#Mf (M sub f) = Fuel moisture content (water weight/dry fuel weight)
#
#Output units: btu/lb
#For heterogeneous fuels this is calculated for each fuel type.  In R this is array compatible but
#may need to be reworked in C++.
HeatOfPreignition <- function(Mf)
{
  Qig = 250 + 1116 * Mf
  return(Qig)
  
  #Metric: Qig = 581 + 2,594Mf
}

#Live / Dead Heat Content:
#  Calculate the weighted (low?) heat content of the live / dead fuel categories.
#Only used by reaction intensity calculation?????
#
#Input variables / parameters:
#hij (h sub ij) = Heat content for live/dead fuel categories (btu/lb).
#fij (f sub ij) = Weighting factors for each fuel type (dimensionless).
#liveDead = An array indicating if each index in each of the other input variables represents a
#  dead (1) or live (2) or fuel categories.
#
#Output units: btu/lb (technically what ever units are input, the same will come out!)
LiveDeadHeatContent <- function(h_ij, f_ij, liveDead)
{
  numFuelTypes = length(f_ij)
  
  #Rothermel 1972 equation 61:
  #hi = Σj fij hij
  h_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    h_i[liveDead[k]] = h_i[liveDead[k]] + f_ij[k] * h_ij[k]
  }
  
  return(h_i)
}

#Move bulk densities up with packing ratios?

#Bulk Density:
#  The bulk density is the mass/wt. of oven dry surface fuel per volume of fuel bed (fuel mass per
#area divided by the fuel bed depth).
#
#Rothermel 1972 equation 40:
#ρb = w_o/δ
#
#Input variables / parameters:
#w_o (w sub o) = Oven dry fuel load (lb/ft^2).  This includes combustible and mineral fractions.
#δ (delta) = fuel bed depth (ft)
#
#Output units: lb/ft^3
BulkDensity <- function(w_o, fuelBedDepth)
{
  rho_b = w_o / fuelBedDepth
  return(rho_b)
}

#EffectiveBulkDensity: Is this needed?

#Mean Bulk Density:
#  The heterogeneous fuel version of the spread equation requires a mean bulk density for the 
#fuels.
#
#Rothermel 1972 equation 74:
#ρb = 1/δ Σi Σj (wo)ij
#For i = 1 to m fuel categories (live vs. dead) and j = 1 to n fuel size classes.
#The original notation includes from and to sum subscripts and bars over rho and w.
#
#Input variables / parameters:
#(w sub o [bar]) = An array of fuel loadings for each fuel category (i: live vs. dead) and fuel size
#class (j) (lb/ft^2).
#  This function could also be written to take something like an array of fuel model objects.  It
#is not clear yet what would be most useful.
#δ (delta) = fuel bed depth (ft)
#
#Output units: lb/ft^3
MeanBulkDensity <- function(fuelLoadings, fuelBedDepth)#Change to w_o_ij
{
  #Sum the individual fuel loadings across elements:
  #No weights are needed since the w_o is expressed as mass per area.
  #The fuel loadings chould be a a 2D array / matrix but the positions have no significance in the
  #calculation  Only a sum of all elments needs to be computed.
  #If we know the size of dimensions we could apply error checking.  The 53 standard fire behavior
  #fuel models have 3 dead and 2 live classes.  That is not a fixed requirement of the Rothermel
  #in theory, but is in practice [I think].
  
  if (length(fuelLoadings) < 2)
  {
    warning("More than one fuel class expected.")
  }
  
  #Explicit matrix: R is row, column: (probably not ideal due to standard fuel models)
  # rows = nrow(fuelLoadings)
  # cols = ncol(fuelLoadings)
  # 
  # totalLoading = 0
  # for (i in 1:rows)
  # {
  #   for (j in cols)
  #   {
  #     totalLoading = fuelLoadings[i,j]
  #   }
  # }
  
  totalLoading = sum(fuelLoadings)
  #rho_b_bar = totalLoading / totalLoading#Wrong!
  rho_b_bar = totalLoading / fuelBedDepth
  
  return(rho_b_bar)
}

#Effective Heating Number:
#  This represents the proportion of a fuel type that is heated to ignition temperature in advance
#of the fire front.  It is a function of SAV.
#  I think of this as the fine fuels will be brought fully to combustion while for a large stick
#only surface would be dried and heated to burn. ?????
#
#Rothermel 1972 equation 14?
#ε = exp(-138/σ)
#
#Input variables / parameters:
#σ / sav = characteristic surface-area-to-volume ratio (ft^2/ft^3) 
#
#Units: Dimensionless.
EffectiveHeatingNumber <- function(sav)
{
  epsilon = exp(-138/sav)
  return(epsilon)
  #Note: Constant = −4.528 for metric units.
}
#For heterogeneous fuels...

#Spread Rate Calculations:

#Albini 1976 modified Rothermel spread model [for homogeneous fuels]:
#  Calculate the steady state spread rate for surface fuels and environmental conditions passed in.
#
#Input variables / parameters:
#  There are 11 input variables in total (see Andrews 2018 table 11), 4 fuel particle
#characteristics, 4 fuel array characteristics, and three environmental.
#
#Fuel particle properties: 
#h = heat content of the fuel class (Btu/lb).
#  A1l the 53 standard fuel models use 8,000 Btu/lb.
#St (S sub t) = Total mineral content (fuel particle property: mineral mass / total dry mass,
#  unitless fraction).  For all standard fuel models this is 5.55% (0.0555).
#Se (S sub e) = effective mineral content (fuel particle property: (mass minerals – mass silica) /
#  total dry mass, unitless fraction).  For all standard fuel models this is 1% (0.01).
#ρp (rho sub p) = fuel particle density (density, originally lb/ft^3)
#  A1l the 53 standard fuel models use 32 lb/ft^3.
#
#Fuel array:
#σ / sav = characteristic surface-area-to-volume ratio (ft^2/ft^3)
#w_o (w sub o) = Oven dry fuel load (lb/ft^2).  This includes combustible and mineral fractions.
#     Sometimes expressed as ton/acre?????
#δ (delta) = fuel bed depth (ft)
#Mx (M sub x) = Moisture of extinction (fraction, water weight/dry fuel weight)
#
#Environmental:
#Mf (M sub f) = Fuel moisture content (water weight/dry fuel weight)
#U = wind speed a midflame height (ft/min)
#tan ϕ = slope steepness, maximum (fraction: vertical rise / horizontal distance, unitless)
#  [For many applications this will need to be converted from degrees.]
#
#useWindLimit = Use the wind limit calculation or not.
#
#Returns: R = rate of spread in ft/min.
Albini1976_Spread <- function(heatContent = 8000, St = 0.0555, Se = 0.01,
                              fuelParticleDensity = 32,#or rho_p
                              sav, w_o, fuelBedDepth, Mx,#fuelLoad
                              Mf, U, slopeSteepness, useWindLimit = TRUE)#Update the name!
{
  #Rate of spread = heat source / heat sink
  
  #R = I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔) / ...
  
  #The heat source (numerator) term represents the heat flux from the fire front to the fuel in
  #front of it:
  #I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔)
  
  #I_Rxi = no wind, no slope propigating flux
  
  #The bulk density is needed to calculate the packing ratio and therefore is used in the numerator
  #and denominator.
  rho_b = BulkDensity(w_o, fuelBedDepth)
  
  packingRatio = PackingRatio(rho_b, fuelParticleDensity)
  optPackingRatio = OptimumPackingRatio(sav)
  
  #Heat source (numerator) = IR𝜉(1 + 𝜙w + 𝜙s)
  #IR = ReactionIntensityAlbini(...)
  
  #Reaction intensity:
  GammaPrime = OptimumReactionVelocity(packingRatio, sav)
  w_n = NetFuelLoad_Albini(w_o, St)
  eta_M = MoistureDampingCoefficient(Mf, Mx)#M_f, M_x?????
  eta_s = MineralDampingCoefficient(Se)#S_e?????
  I_R = ReactionIntensityRothermel(GammaPrime, w_n, heatContent, eta_M, eta_s)
  
  #For debugging:
  # print(paste("GammaPrime =", GammaPrime))
  # print(paste("w_n =", w_n))
  # print(paste("h =", heatContent))
  # print(paste("eta_M =", eta_M))
  # print(paste("eta_s =", eta_s))
  # print(paste("I_R =", I_R))
  
  #Other terms:
  xi = PropagatingFluxRatio(packingRatio, sav)
  phi_s = SlopeFactor(packingRatio, slopeSteepness)
  
  #Apply wind limit check:
  if (useWindLimit)
  {
    if (U/I_R > 0.9)
    {
      U = 0.9 * I_R
      
      #Or Andrews et al. 2013... ?????
      #U = 96.8 * I_R^(1/3)
    }
  }
  
  phi_w = WindFactor(sav, packingRatio, optPackingRatio, U)
  
  
  #The heat sink (denominator) represents the energy required to ignite the fuel in Btu/ft^3:
  #ρbεQig
  
  #rho_b = BulkDensity(w_o, fuelBedDepth)
  epsilon = EffectiveHeatingNumber(sav)
  Qig = HeatOfPreignition(Mf)
  
  #For debugging:
  # print(paste("rho_b =", rho_b))
  # print(paste("epsilon =", epsilon))
  # print(paste("Qig =", Qig))
  # print(paste("Heat sink = ", rho_b * epsilon * Qig))
  
  #There is a sum form for heterogeneous fuels...
  #...(Qig)ij
  #fi and fij are weights...
  
  #Full spread calculation for homogeneous fuels:
  #R = I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔) / ρbεQig
  R = (I_R * xi * (1 + phi_s + phi_w)) / (rho_b * epsilon * Qig)
  #R = (I_R * xi * (1 + phi_s + 17.8)) / (rho_b * epsilon * Qig)#Test
}

#Albini 1976 modified Rothermel spread model for heterogenous fuels:
#w_o -> w_os?
#These parameters should be converted to some object representing a fire behavior fuel mode.
#RAFBFM is a mouthful!

#Input variables / parameters:
#  Some of the input variables differ from the homogenous fuels form in that they are vectors of
#...
#Changed: h, S_t, S_e    SAV
#
#Fuel particle properties: 
#hij (h sub ij) = Heat content of the fuel types (Btu/lb).
#  A1l the 53 standard fuel models use 8,000 Btu/lb.
#(St)ij ((S sub t) sub ij) = Total mineral content (fuel particle property: mineral mass / total dry
#  mass, unitless fraction).  For all standard fuel models this is 5.55% (0.0555).
#(Se)ij ((S sub e) sub ij) = effective mineral content (fuel particle property: (mass minerals –
#  mass silica) / total dry mass, unitless fraction).  For all standard fuel models this is 1% (0.01).
#(ρp)ij ((rho sub p) sub ij) = fuel particle density (density, originally lb/ft^3)
#  A1l the 53 standard fuel models use 32 lb/ft^3.
#
#Fuel array:
#σij (sigma sub ij) / SAV = Characteristic surface-area-to-volume ratios for each fuel type
#(ft^2/ft^3).
#woij ((w sub o) sub ij) = Oven dry fuel load for each fuel type (lb/ft^2).
#     Sometimes expressed as ton/acre?????
#δ (delta) = fuel bed depth (ft)
#(Mx)1 ((M sub x) sub 1) = Dead fuel moisture of extinction (fraction, water weight/dry fuel
#  weight).
#
#Environmental:
#Mfij ((M sub f) sub i) = fuel moisture content for each fuel type (water weight/dry fuel weight)
#U = wind speed a midflame height (ft/min)
#tan ϕ = slope steepness, maximum (fraction: vertical rise / horizontal distance, unitless)
#  [For many applications this will need to be converted from degrees.]
#
#useWindLimit = Use the wind limit calculation or not.  Recent suggestion are that it not be used.
#
#Returns: R = rate of spread in ft/min.
#Name?????
SpreadRateAlbini1976_Het <- function(h_ij = 8000, S_t_ij = 0.0555, S_e_ij = 0.01, rho_p_ij = 32,
                                     SAV_ij,
                                     w_o_ij,
                                     fuelBedDepth,
                                     M_x_1,
                                     M_f_ij,
                                     U, slopeSteepness, useWindLimit = FALSE)
{
  numFuelTypes = length(SAV_ij)
  #Add checking for parameters (input lengths, etc.)...
  liveDead = c(1,1,1,2,2)[1:numFuelTypes]#Global or pass in?
  
  #For the heat content, total mineral content, effective mineral content, and fuel particle density
  #allow the value for all types to to set with a single value:
  h_ij = InitSpreadParam(h_ij, "h_ij", numFuelTypes)
  S_t_ij = InitSpreadParam(S_t_ij, "S_t_ij", numFuelTypes)
  S_e_ij = InitSpreadParam(S_e_ij, "S_e_ij", numFuelTypes)
  rho_p_ij = InitSpreadParam(rho_p_ij, "rho_p_ij", numFuelTypes)
  
  #Terms used in numerator and denominator:_______
  
  #Calculate the weights:
  weights = CalcWeightings(SAV_ij, w_o_ij, rho_p_ij, liveDead)
  #Debugging:
  #print(paste("Weights", weights))
  
  
  #The heat source (numerator) term represents the heat flux from the fire front to the fuel in
  #front of it:
  #I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔)
  #Heat source (numerator) = IR𝜉(1 + 𝜙w + 𝜙s)
  
  #(The bulk density is not used to calculate the packing ratio in the heterogeneous form:)
  meanPackingRatio = MeanPackingRatio(w_o_ij, fuelBedDepth, rho_p_ij)#AKA beta_bar
  
  #For heterogeneous fuels we need to calculate the fuel bed level SAV:
  fuelBedSAV = FuelBedSAV(SAV_ij, weights$f_ij, weights$f_i, liveDead)
  
  optPackingRatio = OptimumPackingRatio(fuelBedSAV)
  
  #Reaction intensity:
  GammaPrime = OptimumReactionVelocity(meanPackingRatio, fuelBedSAV)
  w_n_i = NetFuelLoad_Albini_Het(w_o_ij, S_t_ij, weights$g_ij, liveDead)
  
  #Heat content by live/dead fuel catagory:
  h_i = LiveDeadHeatContent(h_ij, weights$f_ij, liveDead)
  
  #The live fuel moisture of extinction must be calculated:
  M_x_i = c(NA,NA)
  M_x_i[1] = M_x_1
  M_x_i[2] = LiveFuelMoistureOfExtinction(M_f_ij, M_x_1, w_o_ij, SAV_ij, liveDead)
  
  #Damping coefficents:
  eta_M_i = MoistureDampingCoefficient_Het(M_f_ij, M_x_i, weights$f_ij, liveDead)
  eta_s_i = MineralDampingCoefficient_Het(S_e_ij, weights$f_ij, liveDead)
  
  I_R = ReactionIntensity_Het(GammaPrime, w_n_i, h_i, eta_M_i, eta_s_i)
  
  #For debugging:
  # print(paste("GammaPrime =", GammaPrime))
  # print(paste("w_n_i =", w_n_i))
  # print(paste("h_i =", h_i))
  # print(paste("eta_M_i =", eta_M_i))
  # print(paste("eta_s_i =", eta_s_i))
  # print(paste("I_R =", I_R))
  
  #Other numerator terms:
  xi = PropagatingFluxRatio(meanPackingRatio, fuelBedSAV)
  phi_s = SlopeFactor(meanPackingRatio, slopeSteepness)
  
  #Apply wind limit check:
  if (useWindLimit)
  {
    U = WindLimit(U, I_R)
  }
  
  phi_w = WindFactor(fuelBedSAV, meanPackingRatio, optPackingRatio, U)
  
  
  #The heat sink (denominator) represents the energy required to ignite the fuel in Btu/ft^3:
  #ρbεQig
  #The heat sink term is calculated with weights with out calculating epsilon explicitly...
  #Rothermel equation 77:
  #ρbεQig = ρb Σi fi Σj fij[exp(-138/σij)](Qig)ij
  
  rho_b_bar = MeanBulkDensity(w_o_ij, fuelBedDepth)
  Q_ig_ij = HeatOfPreignition(M_f_ij)
  
  #For debugging:
  # print(paste("rho_b_bar =", rho_b_bar))
  # print(paste("Q_ig_ij =", Q_ig_ij))
  
  #We'll do it in two steps:
  #Weight and size class:
  heatSink_i = c(0,0)
  for (k in 1:numFuelTypes)
  {
    heatSink_i[liveDead[k]] = heatSink_i[liveDead[k]] +
      weights$f_ij[k] * exp(-138 / SAV_ij[k]) * Q_ig_ij[k]
  }
  
  #Weight and sum by live/dead category:
  heatSink = rho_b_bar * sum(weights$f_i * heatSink_i)
  
  #For debugging:
  # print(paste("Heat sink =", heatSink))
  
  #Full spread calculation for homogeneous fuels:
  #R = I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔) / ρbεQig
  #  R = (I_R * xi * (1 + phi_s + phi_w)) / (rho_b * epsilon * Qig)
  
  #Full spread calculation for heterogeneous fuels (same as homogeneous in this form):
  #Rate of spread = heat source / heat sink
  #R = I_R𝝃(1 + 𝝓_𝒘 + 𝝓_𝒔) / ρbεQig
  #Without italics:
  #R = I_Rξ(1 + φw + φs) / ρbεQig
  R = (I_R * xi * (1 + phi_s + phi_w)) / heatSink
  
  return(R)
}

#This is a utility function to reduce code repetition in Albini1976_Spread_Het().  It checks the
#length of the parameter value passed in, converts single values to and array, and reports invalid
#lengths.
#Could add type checking.
InitSpreadParam <- function(paramVal, paramName, numFuelTypes)
{
  if (length(paramVal) == 1)
  {
    paramVal = array(data = paramVal, dim = numFuelTypes)
  }
  else if (length(paramVal) != numFuelTypes)
  {
    stop(paste(paramName, "is an unxpected length."))
  }
}
