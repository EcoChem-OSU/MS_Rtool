# Kendrick Mass Defect Analysis for High-Resolution Mass Spectrometry data
-----------------------------------------------------------------------------------
Description
-----------
This code performs a Kendrick Mass Defect (KMD) analysis of user input mass spectral data by creating a feature list via the XCMS peak picking algorithm.
The workflow includes generating a suspect list, generating a feature list, and subsequently creating a KMD plot. 

**Code Outline**

Suspect List Generation:
- Utilizes the creat_KMD_susp_list.v1 and func_calc_residual R codes.
- Combines suspect lists from a users working directory
- Calculates the kendrick mass and kendrick mass defect for each list member for a given neutral loss (ex: CF2, C2F4)

KMD Analysis then occurs in 3 steps (3 individual scripts):  
- 1-PatRoon_XCMS_feature_opti: Optimizes the parameters used for XCMS peak picking and featurelist generation
- 2-PatRoon_featuregroup: Performs peak picking using XCMS and generates a feature list from the data (to fascilitate KMD analysis in step 3)  
- 3-KMD_analysis: Screens feature list for overlap with the suspect list based on exact mass, scans for homologous series of the same mass defect value (KMD analysis)
> [!CAUTION]  
> **DO NOT CONSIDER the KMD Suspect list matches as a comprehensive suspect screening method. Results need to be verified**


Dependency and installation
----------------
R(>=4.3.1)    

RTools 

patRoon(>=2.3.0) Check [patron handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)

**Optional**: Text editor such as Rstudio or Visual Studio Code 

MetaboCoreUtils. To install this package, start R and enter:  

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetaboCoreUtils")
```

**Getting started**
----------------
**Suspect List Generation**
1. Download the creat_KMD_susp_list.v1 and func_calc_residual R codes, store in a folder on your local machine.
2. In the folder, put all suspect lists that are to be combined (**See Example_suspect_list_PRORISKPFAS_short.csv** in Example Input directory)
3. Open the creat_KMD_susp_list.v1 code, select the repeating unit to screen for (ex: CF2, C2F4), and run the code to combine the suspect lists.
4. Output will be in the format "KMD_C2F4_neg_SuspectList_2024-07-19" where 'C2F4' is the repeating unit, 'neg' (or 'pos') reflects the polarity of the suspect list members, and 2024-07-19 is the date the list was generated.

**KMD Analysis:**
1. Download all codes to your local machine, create a directory (folder), and put all codes in the directory (**See Example Input directory**)
2. Create two folders in your directory called "input" and "ouptut"
3. In the input folder, put: (1) a sample list (default name is sample_list_KMD), (2) your suspect list generated previously, (3) all datafiles converted from .raw/.wiff2./etc file format to .mzXML (can convert using MSConvert or another software)
4. Run the scripts in order, taking care to update the parameters in script 2 (2-PatRoon_featuregroup) with the optimized peak picking parameters identified in script 1 (1-PatRoon_XCMS_feature_opti), and to have the featurelist generated from script 2 in the output folder before running script 3 (3-KMD_analysis)


References
------

The current workflow combine/translate divers code from the following sources

Peak picking with XCMS:  
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. Journal of Open Source Software 2022, 7 (71), 4029. [DOI: 10.21105/joss.04029](https://doi.org/10.21105/joss.04029); [github patRoon](https://github.com/rickhelmus/patRoon).

KMD analysis:  
Zweigle, J.; Bugsel, B.; Fabregat-Palau, J.; Zwiener, C. PFDeltaScreen - an open-source tool for automated PFAS feature prioritization in non-target HRMS data. Anal. Bioanal. Chem. 2023, 16 (2), 349–362. [DOI: 10.1007/s00216-023-05070-2](https://doi.org/10.1007/s00216-023-05070-2); [github PFAScreen](https://github.com/JonZwe/PFAScreen/blob/main/KMD_analysis.py)

Theoretical calculation of the residual based on the formula:  [gitHub NIST](https://github.com/usnistgov/NISTPFAS/tree/main/suspectlist/fn/calculate_residual.R)

Citing
------
The full work flow is described in:  
coming soon

Other examples using the workflow include:  


License
-------
This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
