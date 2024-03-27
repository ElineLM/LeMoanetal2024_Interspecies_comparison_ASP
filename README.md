# LeMoanetal2024_Interspecies_comparison_ASP
Codes to reproduce study and figures of article "Inter-species comparison oflife traits relating to amnesic shellfish toxin kinetic infive pectinid species" submitted in Ecological Modelling.


I. Codes for the DEB parameter estimation realised for each individual species
-----------------------------------------------------------------------------
Folder: 1_Individual_species


II. Codes for the multi-species DEB parameter estimation realised
-----------------------------------------------------------------------------
Folder: 2_Multi_species


Divided into 5 folders, both 1st and 2nd folders (1_predict_5species and 2_mydata_5species) need to be added in the path of MATLAB to call them in each specific estimation. Warning: no other files with the same name should be in the different folders, or it will be erase the ones in these folders.

1. Predict
All predict files for the 5 pectinid species.

2. Mydata
All mydata files for the 5 pectinid species.

3. Same_parameters_Placopecten
Scripts for simulations with parameter set of Placopecten magellanicus and physical co-variation rules based on P. magellanicus as reference species. Simulations are done with same parameter values for all species except z which impacts p_Am and E_Hp.

4. Same_parameters_Placopecten_Estimation
From the previous parameter values for Placopecten magellanicus, estimation for all species together, with the same parameter values for all species and still physical co-variation rules.

5. Different_parameters_Estimation
From results of the 4th part, and the new set of parameters, multi-species estimation with increasing weight on v only, pM only and v and pM.
