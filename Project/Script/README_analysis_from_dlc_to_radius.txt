ReadMe for: analysis_from_dlc_to_radius.R

This script is part of the project "A toolkit for the dynamic study of spherical biological objects"
L.S.Burchardt, M. Kehy, M. Gamba, A. Ravignani, W. Pouw

Author of the script: Dr. Lara S. Burchardt

This script automatically estimates circles for all DLC tracking datafiles with the LANDAU algorithm 
as batch processing and adds relevant meta data for this specific dataset. For the pure function without 
any project specific meta data use the script: function_from_dlc_to_radius.R


Start point of the analysis (and therefore necessary prerequesites) are: 
1) a folder with the DLC trackings in csv format as produced by the corresponding DLC code

The folder with the data is chosen interactively to loop through all chosen files. 
The DLC threshold is set, indicating which tracking likelihood we deem sufficient to be included (default 0.6)

Different subfunctions are defined, preparing the data in the necessary way or adding a normalization option to the results
The main functionality is based on the circle estimation with the Landau algorithm (CircleFitByLandau, package: conicfit)
 
Succesfull circle estimation needs a minimum of three datapoints from the DLC trackings, otherwise a NA is given

After circle estimation the meta data on the respective sex and Individual ID for this specific dataset is added 
and the results interactively saved. 

The savename has to be given with the file extension .rds in the name but without any quotation marks as indicated by the 
corresponding prompt