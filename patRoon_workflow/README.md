patRoon workflow:
--------------------------------------------------------------------------
Description
-----------
Suspect analysis workflow using XCMS for peak picking and performing annotation using divers tools (formula, MS2 library matching, Metfrag).  

Getting started
----------------
[Tutorial on patRoon](https://rickhelmus.github.io/patRoon/articles/tutorial.html) and more advance information on  [patron Handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)

1. XCMS need to be optimized using ```PatRoon_XCMS_feature_optiv.R```
2. Suspect screening is performed using ```PatRoon_suspectXCMS.R```  
   Two versions are available:  
   ```DIA``` for data independent acquisition like Swath and   
    ```DDA``` for data dependent acquisition.   
    ```PatRoon_featuregroup.R``` code do only peack picking and filtering (quicker).  

Installation
-----------
- Check the patron [Handbook section installation](https://rickhelmus.github.io/patRoon/handbook_bd/installation.html)
- Latest version could be install from github
``` r 
install.packages("remotes") 
remotes::install_github("rickhelmus/patRoonInst")
``` 
Then, Run.. 
``` r 
patRoon::verifyDependencies()
``` 
To be sure you have all necessary package install

Input
------
- sample list --> check documentation in [patRoon Handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)
- ISTD list as an option: list of internal standarad treated as a suspect list --> check documentation in [patRoon Handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)  
>[!WARNING]
>Use ISTD list option with caution always better to check QA/QC BEFORE running the script.
- XCMS parameter --> defined by the optimization script 
- Filtering setting: ```absMinIntensity```, ```relMinReplicateAbundance``` and ```blankThreshold``` check documentation in [patRoon Handbook filtering](https://rickhelmus.github.io/patRoon/handbook_bd/filtering.html)   
- Mass Defect filtering as an option,very practical for PFAS
- Adduct and formula search parameter 
- File for annotation  
    Suspect list (e.g., [Merged NORMAN Suspect List: SusDat](https://zenodo.org/records/10510477) )  
    MetFrag data based file (e.g., [Pubchem lite](https://doi.org/10.5281/zenodo.4183801) )  
    MS2 library (e.g., [MoNA](https://mona.fiehnlab.ucdavis.edu/downloads) )  

Output
------
Output files included divers dataframes tables and text information

```AA_INFO_RUN_README.txt```: contain all parameter and setting used for the run  
```raw_unaligned_ungrouped.csv```: peak picked features without grouping and RT or m/z alignment  
```raw_aligned_grouped.csv```: peak picked features grouped and aligned  
```featureGroups.csv```: peak picked features grouped, aligned and filtered using Intensity and sample/blank ratio    
```featureGroups_averaged.csv```: peak picked features averaged through the replicate  
```SuspectScreening_all.csv```: results of the suspect screening (full results)  
```SuspectScreening_sample.csv```: results of the suspect screening for sample only including relevant information only.
```ISTD_check.csv```: In case you select ```check.istd <- "YES"``` result a table with istd values treated as suspect. 

Citing
-------
- patRoon   
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. Journal of Open Source Software 2022, 7 (71), 4029. [DOI: 10.21105/joss.04029](https://doi.org/10.21105/joss.04029).  
- Optimization script came from  
Aurich, D.; Diderich, P.; Helmus, R.; Schymanski, E. L. Non-target screening of surface water samples to identify exposome-related pollutants: a case study from Luxembourg. Environ. Sci. Eur. 2023, 35 (1), 94. [DOI: 10.1186/s12302-023-00805-5](https://doi.org/10.1186/s12302-023-00805-5).  
- Tool to retrived Kow from pubchem   
Szöcs, E.; Stirling, T.; Scott, E.; Scharmüller, A.; Schäfer, R. B. webchem : An R Package to Retrieve Chemical Information from the Web. Journal of statistical software 2020, 93. [DOI: 10.18637/jss.v093.i13](https://www.jstatsoft.org/article/view/v093i13).  

Used in
--------
Coming soon!!!
