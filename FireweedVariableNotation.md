<!-- ---------------------------------------------------------------------------
FireweedVariableNotation.md
Joshua M. Rady
Woodwell Climate Research Center
2024

Note: This markdown document uses the GitHub markdown dialect and table extensions.
---------------------------------------------------------------------------- -->

# Fireweed Notation Cross-reference

The notation of many variables in the Rothermel and Albini fire spread model and related publications can not be represented easily in plain text code.  This file provides a cross reference between the fully formated notation of the original papers, two intermediate notations, and the notation used in the code.  See the documentation in the code itself for further details of the notation scheme.  Variable descriptions and units are provided here for convenience and are also included in the internal documentation within the code.


## Fuel Particle Properties: 

| Original | Plain Text | Longhand | Code Notation | Description | Units |
| -------- | ---------- | -------- | ------------- | ----------- | ----- |
| h | h | h | h | Heat content of the fuel type. (Byram uses H.) | Btu/lb, kJ/kg |
| h<sub>ij</sub> | hij | h sub ij | h\_ij | Heat content of the fuel types. | Btu/lb, kJ/kg |
| h<sub>i</sub> | hi | h sub i | h\_i | Heat content for live/dead fuel. | Btu/lb, kJ/kg |
| S<sub>T</sub> | ST | S sub T | S\_T | Total mineral content. | unitless fraction: mineral mass / total dry mass |
| (S<sub>T</sub>)<sub>ij</sub> | (ST)ij | (S sub T) sub ij | S\_T\_ij | An array of total mineral content for each fuel type. | unitless fraction |
| S<sub>e</sub> | Se | (S sub e) | S\_e | Effective mineral content. | unitless fraction: (mineral mass – mass silica) / total dry mass |
| (S<sub>e</sub>)<sub>ij</sub> | (Se)ij | (S sub e) sub ij | S\_e\_ij | Effective mineral content for each fuel type. | unitless fraction: (mineral mass – mass silica) / total dry mass |
| ρ<sub>p</sub> | ρp | rho sub p | rho\_p | Fuel particle density. | lb/ft^3, kg/m^3 |
| (ρ<sub>p</sub>)<sub>ij</sub>| (ρp)ij | (rho sub p) sub ij | rho\_p\_ij | Fuel particle density for each fuel type. | lb/ft^3, kg/m^3 |

## Fuel Array:
| Original | Plain Text | Longhand | Code Notation | Description | Units |
| -------- | ---------- | -------- | ------------- | ----------- | ----- |
| σ | σ | sigma | SAV | Characteristic surface-area-to-volume ratio. | ft^2/ft^3, cm^2/cm^3 |
| σ<sub>ij</sub> | σij | sigma sub ij | SAV\_ij | Characteristic surface-area-to-volume ratios for each fuel type. | ft^2/ft^3, cm^2/cm^3 |
| w<sub>o</sub> | wo | w sub o | w\_o | Oven dry fuel load.  This includes combustible and mineral fractions. | lb/ft^2, kg/m^2 |
| (w<sub>o</sub>)<sub>ij</sub> | (wo)ij | (w sub o) sub ij | w\_o\_ij | An array of oven dry fuel load for each fuel type. | lb/ft^2, kg/m^2 |
| w<sub>n</sub> | wn | w sub n | w\_n | Net fuel load for a single fuel component. | lb/ft^2, kg/m^2 |
| (w<sub>n</sub>)<sub>i</sub> | (wn)i | (w sub n) sub i | w\_n\_i | Net fuel load for live/dead fuel categories. | lb/ft^2, kg/m^2 |
| δ | δ | delta | fuelBedDepth | Fuel bed depth. | ft, m |
| M<sub>x</sub> | Mx | M sub x | M\_x | Moisture of extinction. | fraction: water weight/dry fuel weight |
| (M<sub>x</sub>)<sub>i</sub> | (Mx)i | (M sub x) sub i | M\_x\_i | Moisture of extinction each fuel category. | fraction: water weight/dry fuel weight |
| (M<sub>x</sub>)<sub>1</sub> | (Mx)1 | (M sub x) sub 1 | M\_x\_1 | Dead fuel moisture of extinction. | fraction: water weight/dry fuel weight |
| ρ<sub>b</sub>| ρb | rho sub b | rho\_b | Fuel array bulk density. | lb/ft^3, kg/m^3 |

## Environmental:
| Original | Plain Text | Longhand | Code Notation | Description | Units |
| -------- | ---------- | -------- | ------------- | ----------- | ----- |
| M<sub>f</sub> | Mf | M sub f | M\_f | Fuel moisture content. | fraction: water weight/dry fuel weight |
| (M<sub>f</sub>)<sub>ij</sub>| (Mf)ij | (M sub f) sub ij | M\_f\_ij | Fuel moisture content for for each fuel type. | fraction: water weight/dry fuel weight |
| U | U | U | U | Wind speed at midflame height. | ft/min, m/min |
| tan ϕ | tan ϕ | tan phi | slopeSteepness | Slope steepness, maximum. | unitless fraction: vertical rise / horizontal distance |

## Weights:
|  Original | Plain Text | Longhand | Code Notation | Description | Units |
| -------- | ---------- | -------- | ------------- | ----------- | ----- |
| f<sub>ij</sub> | fij | f sub ij | f\_ij | Weighting factors for each fuel type. | dimensionless |
| f<sub>i</sub> | fi | f sub i | f\_i | Weighting factors for each fuel live/dead category. | dimensionless |
| g<sub>ij</sub> | gij | g sub ij | g\_ij | Net fuel load weights (Albini 1976). | dimensionless |

## Other:
|  Original  | Plain Text | Longhand | Code Notation | Description | Units |
| -------- | ---------- | -------- | ------------- | ----------- | ----- |
| Γ' | Γ' | Gamma prime | GammaPrime | Optimum reaction velocity. | min^-1 |
| β | β | beta | packingRatio | Packing ratio, the fraction of the fuel bed volume occupied by fuel. | dimensionless |
| $`\bar{\beta}`$ | beta[bar] | beta bar | meanPackingRatio | Mean packing ratio.  | dimensionless ratio |
| β<sub>op</sub> | βop | beta sub op | optPackingRatio | Optimum packing ratio. | dimensionless |
| ϕ<sub>w</sub> | ϕw | phi sub w | phi\_w | The wind factor. | dimensionless |
| ϕ<sub>s</sub> | ϕs | phi sub s | phi\_s | The slope factor. | dimensionless |
| I<sub>R</sub> | IR | I sub R | I\_R | Reaction intensity. | Btu/ft^2/min, kj/m^2/min |
| t<sub>r</sub> | tr | t sub r | t\_r | Residence time. | min |
| I<sub>B</sub> | IB | I sub B | I\_B | Byram's fireline intensity. | Btu/ft/s, kW/m |
| η<sub>M</sub> | ηM | eta sub M | eta\_M | Moisture damping coefficient. | unitless |
| (η<sub>M</sub>)<sub>i</sub>| (ηM)i | (eta sub M) sub i | eta\_M\_i | Moisture damping coefficient for live/dead fuel categories. | unitless |
| η<sub>s</sub>| ηs | eta sub s | eta\_s | Mineral damping coefficient. | unitless |
| (η<sub>s</sub>)<sub>i</sub>| (ηs)i | (eta sub s) sub i | eta\_s\_i | Mineral damping coefficient for live/dead fuel categories. | unitless |
