Practical Semiquantification Strategy for Estimating Suspect Concentrations: R-script
-----------------------------------------------------------------
Description
-----------
Perform semiquantification average calibration using the area count of the targets
divided by the average area count of surrogate (deutered intern standard) in 
function of target concentration in units of nmoles/L. Use a weighted linear regression model.

The scripts in this repository organized into 2 steps
1. cali_pract_semiquantif: performed a calibration of the semiquantitive model.
2. pred_pract_semiquantif: applied the calibrated model to suspect area count to estimate the concentration.

Terminology
----------
Target: refer to the compound of interest.
Surrogate: refer to the isotopic labeled compounds used as intern standard to normalize the response factor.

Input Files Needed
------------------
> [!IMPORTANT]
> **CHECK input examples. Should be a csv file with no special character (e.g., - % ? ! /). If the data came from excel be sure to change N/A by NA.**
- Calibration curve (example: 060121_ccal.csv): first column content the exact concentration value and follow column contain area count (one per chemical species) for each concentration. 
- List of the targets and surrogates (example: target_surrogatet_list.csv): identified the name of each chemical species, formula necessary for mol calculation and type of species (_target_ or _surrogate_)
- Suspects area count (example: susp_area_count.csv): area count one column per suspect and one raw per samples
> [!IMPORTANT]
> **CHECK THAT all chemical names are correctly and similarly written in the _calibration curve_ and _list of the targets and surrogates_ files. _Suspects area count_ MUST contain the same list of surrogate than in _list of the targets and surrogates_ used to calibrate the model. All species not identified in this list are considered as suspect.**


Getting Started
----------------
1. Clone the R script/input example and save it on an empty folder (your work directory).
   - If downloaded file goes to "Downloads", then copy/move files into new empty folder.
2. 

3. The script create automatically an output folder containing all info of the calibrated model.


Citing
-------
Cao, D.; Schwichtenberg, T.; Duan, C.; Xue, L.; Muensterman, D.; Field, J. Practical semiquantification strategy for estimating suspect per- and polyfluoroalkyl substance (PFAS) concentrations. J. Am. Soc. Mass Spectrom. 2023, 34 (5), 939–947. DOI: 10.1021/jasms.3c00019.
