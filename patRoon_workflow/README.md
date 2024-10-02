patron workflow:
--------------------------------------------------------------------------
Description
-----------
Suspect analysis workflow using XCMS for peak picking and performing annotation using divers tools (formula, MS2 library matching, Metfrag, ..)
Suspect false positive annotation are check using a retention time relationship with the logKow. 


Getting started
----------------
[Tutorial on patRoon](https://rickhelmus.github.io/patRoon/articles/tutorial.html) and more advance information on  [patron Handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)

1. XCMS need to be optimized using 'PatRoon_XCMS_feature_optiv.R'
2. Suspect screening is performed using 'PatRoon_suspectXCMS.R'

Installation
-----------
- Check the patron Handbook section installation
- Run 
``` r 
patRoon::verifyDependencies()
``` 
To be sure you have all necessary package install


Input
------


Output
------
Output files included divers dataframes tables

raw_unaligned_ungrouped.csv: peak picked features without grouping and RT or m/z alignment
raw_aligned_grouped.csv: peak picked features grouped and aligned
featureGroups.csv: peak picked features grouped, aligned and filtered using Intensity and sample/blank ratio  
featureGroups_averaged.csv: peak picked features averaged through the replicate
SuspectScreening_all.csv: results of the suspect screening (full results)
SuspectScreening_sample.csv: results of the suspect screening for sample only including relevant information only.
SuspectScreening_sample_RTKow_check.csv: re-assignement of _sample table to level of confidence <=3c to 4 if not consistent with logKow = function(RT)







Citing
-------
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. Journal of Open Source Software 2022, 7 (71), 4029. [DOI: 10.21105/joss.04029](https://doi.org/10.21105/joss.04029)  
Szöcs, E.; Stirling, T.; Scott, E.; Scharmüller, A.; Schäfer, R. B. webchem : An R Package to Retrieve Chemical Information from the Web. Journal of statistical software 2020, 93. [DOI: 10.18637/jss.v093.i13](https://www.jstatsoft.org/article/view/v093i13)


Used in
--------
Coming soon!!!
