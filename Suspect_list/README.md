Suspect List: working version of divers tool to create/merge suspect list
--------------------------------------------------------------------------
Description
-----------
Compilation of divers script used for suspect list.

> [!NOTE]
> ** CHECK CAREFELY THE DEPENDENCY.**
<<<<<<< HEAD
> Info present in section "DEPENDENCY AND FIRST INSTALL" in each script.
=======
> Info present in section "DEPENDENCY AND FIRST INSTALL" in each script.   
> Some dependency are not present in the CRAN and need to be manually added from mention source.
>>>>>>> 57deaf3139a094906e267266af53978a72b68631

Script catalogue
----------------
- 1.Fill up_list_using_smile or name_v2.2.R 
- 2.merge_suspect_list_v4.1.R: Uniffied, fill-up, remove duplicate and create specific pos and neg mod list from suspect list found in depository (i.e., Norman, Nist, EPA, ...) 
- 3.TPs_predV1.1.R: generate transformation suspect list based on the patRoon 2.0 (Helmus, 2022). Several tools.
- 3.TP_pred_envipath -- coming soon 
- Getfrom_mol_v4.2.R: generate suspect list out of mol file drawing
- PFAS_acronyms -- coming soon 

Input
------
Input:
<<<<<<< HEAD
- Need a table with at minimal the SMILES, INCIKEY or CAS number   
=======
- Need a table with at minimal the SMILES, INCIKEY or CAS number. See input exemple.   
>>>>>>> 57deaf3139a094906e267266af53978a72b68631
>[!WARNING]
> If you work with CAS number as an input you might have to double check the results. working with SMILES or INCHIKEY are more cleaner.
  
REFERENCES
----------
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. J. Open Source Softw. 2022, 7 (71), 4029. DOI: 10.21105/joss.04029.

Hafner, J., Lorsbach, T., Schmidt, S., Brydon, L., Dost, K., Zhang, K., Fenner, K., & Wicker, J (2024). Advancements in Biotransformation Pathway Prediction: Enhancements, Datasets, and Novel Functionalities in enviPath. Journal of Cheminformatics, 16, 93, https://doi.org/10.1186/s13321-024-00881-6

