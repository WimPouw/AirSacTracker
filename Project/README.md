# A toolkit for the dynamic study of elastic kinematics

This is the READme for the github project folder of the manuscript "A toolkit for the dynamic study of elastic kinematics".
It contains necessary code and data to re-run analyses reported on in the manuscript. If you want to run the tracking procedure described in the manuscript with your own data,
please see the toolkit folder (AirSacTracker/toolkit).

## Table of Contents
- Overview
- Folder Structure
- Proof of Concept Analyses
- ToolkitR
- Hough Example Frog
- DeepLabCut
- Contributing
- License
- Acknowledgments

## Overview

In the project, we automatically tracked the radii of Siamang airsac in two different contexts, to study airsac inflation and its influence on acoustic parameters
of produced vocalizations. We established the tracking approaches, which are reported in the toolkit folder of this repository. In this folder, we provide
all other relevant code and information to re-run all other analysis, reported on in the manuscript, such as the comparison between manually tracked circles 
and automatically tracked circles with the newly developed approaches (Hough method and and DLC+ method). 

## Folder Structure

The different folders contain information and codes for different parts of the project. It contains the following folders:

- DeepLabCut --> instructions and data of the tracked airsacs for this project
- Comparison_manual_auto --> data and code for the comparison between manually tracked airsacs and automatically tracked airsacs
- Proof_of_concept_analyses --> data and code for the first proof of concept analysis, studying boom calls
- Proof_of_concept_analyses2 --> data and code for the second proof of concept analysis, studying the radius of boom calls and the acoustic of subsequent barks
- ToolkitR --> contains equivalent R code for the key step of DLC+ (which is implemented in python instead in the DLC+ pipeline)
- Hough_example_frog --> contains a hough tracking example for a different type of airsac, namely of a frog
- Figures --> figures that were used in the manuscript

## Proof of Concept Analyses

Data snippets for the kinematic-acoustic analysis were sampled opportunistically from the full dataset and are included in the corresponding folders (see details below).
The two analyses demonstrate how the developed tracking approaches work and can be utilized for analyses.

### Proof 1

This folder contains the main R-Code for this analysis: "analysis_acoustics_radius_boom.R" as well as the relevant used data in snippets > audio > multi_checked_ok.
This is the audio data from the multisource microphones, that were manually checked not to include any loud noises or more than one siamangs vocalizations. 
Hough trackings are to be found in csv files in snippets > tracked_hough. 
And DLC trackings are stored in csv files in snippets > tracked_dlc, here you can also find the video snippets with the tracked keypoints.

### Proof 2

This folder follows the same structure as for Proof 1. It contains the main R-Code: "analysis_acoustics_bark_radius_boom_project.R" as well as the relevant used data 
in snippets > short_sequences > audio > bark_multi_checked. Again these are the brak snippets from the multisource microphones, manually checked for noises. 
The long sequences are the full great call sequences, the short sequences are the individual snippets of booms and barks in these long sequences. 

## ToolkitR

This folder contains equivalent R code for estimating circles from deeplabcut data. We dediced in the end to do everything in python for the toolkit.

## Hough Example Frog

This is a small example of how well the Hough Transformation approach of tracking circles works in a different species, here a frog. 

## DeepLabCut

This folder contains information on the tracked key points, as well as the data for the evaluation set and the training set (see manuscript for details). 

## Usage

All codes were run in R 4.2.3. No extra installation of packages is necessary, as all necesary packages are installed and loaded automatically, when running the codes. 
All file paths expect the AirSacTracker main folder to be the  working directory. 

## License
TBA

## Acknowledgments

We would like to thank the student assistant Mouaz Morad for his manual annotation of the siamang air sacs. 
We would like to thank Daniel Anthes for his comments on the toolkit code.
