# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 17:33:27 2023

@author: Wim Pouw
"""
import os
from os import listdir
from os.path import isfile, join
import cv2
import pandas as pd
import math

tsfol ='./timeseries/'                       #this is where your timeseries are with the same name as the complementary video
vidfol = './original_videos/'                #this is where the original videos are, that need a circle added
outfol = './output_videos_with_circles/'     #this is where you can collect your output
toprocess = os.listdir(tsfol )               #list all the time series files

for tt in toprocess:
    ts = pd.read_csv(tsfol + tt) #get the time series
    idname = tt[0:len(tt)-4] #remove the .csv
    vidloc =  vidfol+idname+'.mp4' #add mp4 (we assume we only process mp4s!)
    cap = cv2.VideoCapture(vidloc) #open the video
    frameWidth = cap.get(cv2.CAP_PROP_FRAME_WIDTH)   #get the framewidth, and use it for the new video
    frameHeight = cap.get(cv2.CAP_PROP_FRAME_HEIGHT) #get the framewidth, and use it for the new video
    fps = cap.get(cv2.CAP_PROP_FPS)   #fps = frames per second
    #what should we write to?
    out = cv2.VideoWriter(outfol+idname+'_circle.mp4',cv2.VideoWriter_fourcc(*'MP4V'), fps, 
                      (int(frameWidth), int(frameHeight)))
    print('working on video: ' +idname)
    while(cap.isOpened()):
        ret, frame = cap.read()
        if ret == False:
            break
        frame_number = int(cap.get(cv2.CAP_PROP_POS_FRAMES))
        index_var = (ts['frame'] == frame_number) #get the index of the timeseries for the current frame number
        dat = ts.loc[index_var,:] #get the slice of data for this frame
        if math.isnan(dat['radius']) == False: #only draw a circle when there a no NaN's
            cv2.circle(frame, (int(dat['x']), int(dat['y'])), int(dat['radius']), (200, 0, 0), 2) #draw circle
        out.write(frame) #save it into a new frame

# cleaning up
out.release() #release the output video
cap.release() #release the original video
print('We are all done, go look into your output folder: ' + str(os.path.abspath(outfol)))