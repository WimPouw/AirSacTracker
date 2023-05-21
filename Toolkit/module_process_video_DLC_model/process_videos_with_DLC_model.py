# -*- coding: utf-8 -*-
"""
Created on Wed May 10 18:16:54 2023

@author: wimpo
"""
import deeplabcut
import os
from os.path import isfile, join

# Load the pre-trained model
config_path = './trained_model_and_metainfo/config.yaml'

output_dir = './results/'

# set videofolder
videofolder = './videos/'

#version 2, using videos #################### loading in the videos
vids = [f for f in os.listdir(videofolder) if isfile(join(videofolder, f))]
vidlist = []

for i in vids: #add the image folder name to get the full path
    video_path = videofolder+i

    # Analyze the video using the pre-trained model
    deeplabcut.analyze_videos(config_path, [video_path],save_as_csv=True, videotype='.mp4', destfolder=output_dir)
    
    # Create labeled videos
    deeplabcut.create_labeled_video(config_path, [video_path], destfolder=output_dir)
    
    # Optional: Create labeled frames
    deeplabcut.extract_frames(config_path, [video_path], destfolder=output_dir)