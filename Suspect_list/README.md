MalbacR: Script to test batch correction methods for features normalization
--------------------------------------------------------------------------
Description
-----------
Provide metric to evaluate if normalization of the feature is needed (Ding, 2022).  
Test several normalization methods for HRMS data provide (Leach, 2023).

Script is an adapted version based on [malbacR repository](https://github.com/pmartR/malbacR) with some metric to help you over the selection of the method.

Getting started
----------------
Tutorial on [malbacR repository](https://github.com/pmartR/malbacR) and on [pmartR repository](https://pmartr.github.io/pmartR/)

Installation
-----------
``` r
install.packages("devtools")
```
Then need to install Rtools compatible with you R version.
```
devtools::install_github("pmartR/malbacR")

devtools::install_github("pmartR/pmartR@*release")
``` 


Input/Output
------------
Input:
- FeatureGroups.txt file coming from a peak picking method 
- Data need to have some QC (pooled samples)  
  
Output:  
- %Relative standard deviation (RSD) of feature intensity for all QC (QC_boplot.png)
- %RSD of feature intensity between replicate for all samples classified by time 
- PCA for each techniques
- Normalized data for each techniques
- cumulative frequency of the RSD for each techniques

Citing
-------
Leach, D. T.; et al. Anal. Chem. 2023, 95 (33), 12195-12199. [https://doi.org/10.1021/acs.analchem.3c01289](https://doi.org/10.1021/acs.analchem.3c01289).  
Ding, X.; et al. Anal. Chem. 2022, 94 (21), 7500-7509. [https://doi.org/10.1021/acs.analchem.1c05502](https://doi.org/10.1021/acs.analchem.1c05502)
