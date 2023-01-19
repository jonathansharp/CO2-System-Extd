# CO2-System-Extd

<a href="https://zenodo.org/badge/latestdoi/198885961"><img src="https://zenodo.org/badge/198885961.svg" alt="DOI"></a> [![View CO2SYSv3 for MATLAB on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/78378-co2sysv3-for-matlab)

## ABOUT

This repository includes software compatible with MATLAB and GNU Octave for calculating marine CO2 system variables (CO2SYS.m), computing partial derivatives of calculated CO2 system variables with respect to inputs (derivnum.m), and propagating uncertainties for CO2 system calculations (errors.m). This software performs similarly to previously released versions of CO2SYS.m (v1: https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/; v2: https://github.com/jamesorr/CO2SYS-MATLAB), and includes the following extended capabilities, additions, and bug fixes (among other minor changes):
 
1) Can accept input parameters of [CO3], [HCO3], and [CO2], and propagate their uncertainties
2) Includes NH3 and HS as alkalinity contributors, and propagates their uncertainties
3) Uses separate inputs to specify choices for characterizations of K1K2, KSO4, KF, and TB
4) Does not evaluate input parameters equal to -999 or NaN
5) Exits pH iteration loops that do not converge and indicates where a problem occurred
6) Provides exactly identical pH results for a given input line, no matter the other lines of input parameters (this is not necessarily the case for prior versions of CO2SYS.m)
7) Uses an updated definition of the ideal gas constant (https://physics.nist.gov/cgi-bin/cuu/Value?r)
8) Fixes bugs in CO2SYS.m Revelle factor calculation and derivnum.m output conditions
9) Includes K1 and K2 constants defined by Sulpis et al. (2020), K2 constant defined by Schockman and Byrne (2021), KF constant defined by Perez and Fraga (1987), and KSO4 constant of Waters and Millero (2013)
10) Determines initial pH in iterative solvers using the approach of Munhoven (2013), detailed further in Humphreys et al. (submitted), rather than simply using an initial estimate of 8.0 each time.
11) Obtains free scale pH properly within iterative pH solvers no matter the input scale, rather than making the simplification that input pH is always on the total scale.
12) Includes substrate-inhibitor ratio (Bach, 2015) as an output argument from CO2SYS.
13) Calculates uncertainties at output conditions that are associated with equilibrium constants with respect to equilibrium constants at output conditions, rather than input conditions as previously. This essentially assume pK uncertainty is constant regardless of temperature and pressure.
14) Calculates derivatives and errors for the Revelle factor.
15) errors.m includes optional calcium concentration uncertainty input as discussed in Dillon et al. (2020)

Also included in this repository is a routine to compare CO2SYSv3 to CO2SYSv2 (compare_versions.m), a routine to calculate total concentrations of conservative elements (Na, Mg, Cl, etc.) from CO2SYS output (TOTALS.m), and an example function to run CO2SYSv3 and plot some of the output (example_CO2SYS.m).

## HISTORY

CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, later adapted for MS Excel and MATLAB by Pierrot (2006). The code was vectorized, refined, and optimized for computational speed by van Heuven (2011). Options for error propagation were added by Orr et al. (2018). This software builds upon those previous versions.

## INSTALLATION AND USE

Download the files in this repository and place them in a directory that is in the MATLAB search path. Or add the directory where they are located to the search path (https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html).

To perform CO2 system calculations, use CO2SYS.m as directed in the function's help text. All sub-routines that will be called upon are contained within the CO2SYS.m function.

To propagate uncertainties in CO2 system calculations, use errors.m as directed in the function's help text. The errors.m routine will call upon CO2SYS.m and derivnum.m. As such, this version of errors.m is compatible only with CO2SYSv3.

To compute partial derivatives of calculated CO2 system variables with respect to input parameters, use derivnum.m as directed in the function's help text. The derivnum.m routine will call upon CO2SYS.m. As such, this version of derivnum.m is compatible only with CO2SYSv3.

## CITATION

The full citation for CO2SYSv3 (Sharp et al., 2020) is given below. Cite this version if using CO2SYSv3 for CO2 system calculations or propagating errors in CO2 system calculations using the extended errors.m or derivnum.m routines provided here.

If using any CO2SYS program for CO2 system calculations, cite also the original CO2SYS DOS work of Lewis and Wallace (1998).

If using the CO2SYS MATLAB program for CO2 system calculations, cite also the work of van Heuven et al. (2011).

If using the derivnum.m and/or errors.m programs for CO2 system error propagations, cite also the work of Orr et al. (2018).

## REFERENCES

Bach, L. T. (2015). Reconsidering the role of carbonate ion concentration in calcification by marine organisms. Biogeosciences 12(16), 4939–4951.

Dillon, W. D. N., Dillingham, P. W., Currie, K. I., McGraw, C. M., 2020. Inclusion of uncertainty in the calcium-salinity relationship improves estimates of ocean acidification monitoring data quality. Marine Chemistry 226, 103872.

Humphreys, M.P., Lewis, E.R., Sharp, J.D., Pierrot, D. (2022). PyCO2SYS: marine carbonate system calculations in Python. Geoscientific Model Development 15, 15-43.

Lewis, E., Wallace, D. W. R., 1998. Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Munhoven, G., Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388, 2013

Orr, J.C., Epitalon, J.-M., Dickson, A. G., Gattuso, J.-P., 2018. Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

Sharp, J.D., Pierrot, D., Humphreys, M.P., Epitalon, J.-M., Orr, J.C., Lewis, E.R., Wallace, D.W.R. (2021, May 20). CO2SYSv3 for MATLAB (Version v3.2.0). Zenodo. http://doi.org/10.5281/zenodo.3950562

Sulpis, O., Lauvset, S. K., and Hagens, M., 2020. Current estimates of K1* and K2* appear inconsistent with measured CO2 system parameters in cold oceanic regions. Ocean Science Discussions, 1-27.

van Heuven, S., Pierrot, D., Rae, J.W.B., Lewis, E., Wallace, D.W.R., 2011. MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.
