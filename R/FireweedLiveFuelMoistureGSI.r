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
