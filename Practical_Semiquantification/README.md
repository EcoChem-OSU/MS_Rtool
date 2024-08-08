Practical Semiquantification Strategy for Estimating Suspect Concentrations: R-script
-----------------------------------------------------------------
Description
-----------
Perform semiquantification average calibration using the area count of the targets
divided by the average area count of surrogate (deutered intern standard) in 
function of target concentration in units of mol/L.  
The model used a weighted linear regression model. 

The process is divided into 2 R scripts that should run in numerical order as identified in the file name.  
1.cali_pract_semiquantif.R: performed a calibration of the semiquantitive model.  
2.pred_pract_semiquantif.R: applied the calibrated model to suspect area count to estimate the concentration.  
The file FUN_Pract_semiquantif.R: contain home made function necessary for the ```cali_pract_semiquantif``` but should NOT be modified.

Terminology
----------
Target: refer to the compound of interest.
Surrogate: refer to the isotopic labeled compounds used as intern standard to normalize the response factor.

Parameter - input of the workflow 
-------------------------  
The script could handle input calibration curve (ccal input) data in mol/L or g/L . In case of g/L unit the formula of each compound provide in the target_surrogatet_list input is used to convert the g/L into mol/L. In.mass option defined if mol/L or g/L unit is used.  

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
2. Open the script 1.cali_pract_semiquantif.R and modified the PARAMETER needed.
3. Do not modified past the PARAMETER section in the R script.
4. Run the script (Select all and run in once).
5. The script create automatically an output folder containing:
	- Model_Perfomance: text file with the performance and detail parameter of the model
	- cali_lm_weight: calibrate model (R-object)
	- cali_curve_copy:
	- cali_dataset
	- predicted_interval
	- target_surrogatet_list_copy
	- Semiquant_cali_plot
	- Semiquant_residual_plot
6. CHECK OUTPUT 
7. Open the script 2.pred_pract_semiquantif.R and modified the PARAMETER needed.
8. Do not modified past the PARAMETER section in the R script.
9. Run the script (Select all and run in once).
10. A pop-up window should appear to select you data.
11. Check the result in: suspect_pred_conc.csv
> [!CAUTION]
> **The results in the suspect_pred_conc.csv are express in mol values similar as the calibrate model Check unit in the Model_Perfomance file.**

Citing
-------
Cao, D.; Schwichtenberg, T.; Duan, C.; Xue, L.; Muensterman, D.; Field, J. Practical semiquantification strategy for estimating suspect per- and polyfluoroalkyl substance (PFAS) concentrations. J. Am. Soc. Mass Spectrom. 2023, 34 (5), 939–947. [DOI: 10.1021/jasms.3c00019](https://doi.org/10.1021/jasms.3c00019).
