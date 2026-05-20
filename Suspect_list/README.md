Suspect List: working version of divers tool to create/merge suspect list
--------------------------------------------------------------------------
Description
-----------
Compilation of divers script used for suspect list.

> [!NOTE]
> ** CHECK CAREFELY THE DEPENDENCY.**
> Info present in section "DEPENDENCY AND FIRST INSTALL" in each script.

>[!INFO]
> 


Script catalogue
----------------
1.Fill up_list_using_smile or name_v2.2.R
2.merge_suspect_list_v4.1.R: Uniffied, fill-up, remove duplicate and create specific pos and neg mod list from suspect list found in depository (i.e., Norman, Nist, EPA, ...)
3.TPs_predV1.1.R: generate transformation suspect list based on the 
3. 
Getfrom_mol_v4.2.R: generate suspect list out of mol file drawing
PFAS_acronyms

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

REFERENCES
----------
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. J. Open Source Softw. 2022, 7 (71), 4029. DOI: 10.21105/joss.04029.

Hafner, J., Lorsbach, T., Schmidt, S., Brydon, L., Dost, K., Zhang, K., Fenner, K., & Wicker, J (2024). Advancements in Biotransformation Pathway Prediction: Enhancements, Datasets, and Novel Functionalities in enviPath. Journal of Cheminformatics, 16, 93, https://doi.org/10.1186/s13321-024-00881-6

