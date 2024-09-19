/***************************************************************************************************
FireweedMetUtils.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
2024
Reference: Proj. 11 Exp. 20

Description:---------------------------------------------------------------------------------------
	This is part of the Fireweed wildfire code library.
	This file defines utility functions for calculatina and converting meteorological properties.

References:----------------------------------------------------------------------------------------

F. W. Murray.
On the Computation of Saturation Vapor Pressure.
Journal of Applied Meteorology and Climatology 6(1): 203-204, 1967. DOI: https://doi.org/10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2
	This paper compares the calculations of Goff and Gratch (1946) and Tetens (1930).  I have not
been able to find the original text of Tetens 1930, which is in German.

Arden L. Buck.
New Equations for Computing Vapor Pressure and Enhancement Factor.
Journal of Applied Meteorology and Climatology 20(12): 527-1532, 1981. https://doi-org.library.proxy.mbl.edu/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2
	This is the original source of the Buck equation.

Arden L. Buck.
Model CR-1A hygrometer with autofill operating manual.
Buck Research Instruments LLC: Aurora, CO, USA. 2012 (1996).
Obtained from https://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf.
	This manual contains an appendix in which Buck provides a set of equations converting different
humidity related quantities, including new parameters for the Buck equations and a enhancement
factor equation, which can replace the table values from the original paper.
	A manual is an unusual reference but it cited by several sources.  The manual was revised in
2012 but this appendix has was last revised in 1996 and sources cite it as Buck 1996.

***************************************************************************************************/

#include "FireweedMetUtils.h"

//Saturation Vapor Pressure:------------------------------------------------------------------------
