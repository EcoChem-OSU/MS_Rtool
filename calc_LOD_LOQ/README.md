Limit of Detection (LOD) and Limit of Quantification (LOQ) Method 
-----------------------------------------------------------------
Description
-----------
Perform LOD and LOQ calculation as reported in Vial et al., (1999) using the residual standard deviation of a weighted regression between area count in function of the concentration.

> [!NOTE]
> ** CHECK CAREFELY THE EXPERIMENTAL DESIGN BEFORE STARTING.**

Experimental Design
-----------
- Area count (or intensity) should be measure for a range of expected concentration.
- Vial et al. recommend to have at least 7 measurable data points to make the calculation.
- Ideally two points are lower than the expected LOD and 5 greater than LOD.
- Data should be spread within 200 fold (e.g., 1 to 200 or 200 to 40000) 
  to follow paper experimental setting. Better the data are closer to expected LOD better the evaluation is.
- Do not include concentration zero and area zero into the calculation.
- Data point do not necessarly need to be equaly separate (weight normalization).

Getting started
----------------
1. Clone the R script/input example and save it on an empty folder (your work directory).
   - If downloaded file goes to "Downloads", then copy/move files into new empty folder ("Documents" works)
2. Check Input file:
   - Structure of the input: first column (could be selected using nx in the script) with the exact concentration value
                             sequence of column containing the area count (ny in the script) or intensity measured (one per
                             chemical species) for each concentration (NA could be set for no value).  
   - Should be a csv file with no special character (e.g., - % ? ! /).
   - If the data came from excel be sure to change N/A by NA.
   - Double check given example of input if not sure ... be sure to have similar structure (header, column) file.
3. Do not modified past the PARAMETER section in the R script.
   - Ususally you need to change you work directory (workdir) only
4. Run the script (Select all and run in once). 
5. Find out an Output file with the calculation in the workdir you set in the script. The output will contain LOD, LOQ, R2 and       RMSE of the model for each species you introduce in the input. Output LOD and LOQ unit are same as the input.

Citing
-------
Vial, J.; Jardy, A. Experimental Comparison of the Different Approaches 
To Estimate LOD and LOQ of an HPLC Method. 
Anal. Chem. 1999, 71 (14), 2672-2677. [DOI:10.1021/ac981179n](https://pubs.acs.org/doi/10.1021/ac981179n).
