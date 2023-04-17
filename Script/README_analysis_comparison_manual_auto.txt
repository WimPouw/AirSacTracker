ReadMe for: analysis_comparison_manual_auto.R

This script is part of the project "A toolkit for the dynamic study of spherical biological objects"
L.S.Burchardt, M. Kehy, M. Gamba, A. Ravignani, W. Pouw

Author of the script: Dr. Lara S. Burchardt

In this script automatically tracked radii from different approches are compared to manually tracked radii, 
to decide if tracking success is sufficient and which approach works best

Start point of the analysis (and therefore necessary prerequesites) are: 
1) radii values from Hough Transform or DLC trackings and subsequent 
circle estimation in csv files
2) datasubset of manually tracked radii in csv files

Those files are produced by other scripts in the project: analysis_from_dlc_to_radius.R

End point of the script are correlation values between manual and automaticlly tracked radii + visualizations of those

-------
Details:

- datafiles of manually tracked radii and dlc trackings are loaded
- loading of hough transform values is implemented in loop, to enable automatic analysis of batches of results
produced for parameter optimization of HOugh Transfom parameters (see corresponding Python script for parameter optimization)

correlations between tracked radii and manual radii are calcualted with specific conditions (i.e. radii below 100 px
are filtered from manual trackings as this indicated "no tracking")
Furthermore automatically tracked radii are filtered to include the same range of maximal radii as manually tracked radii. 

Plots of the correlation are produced
