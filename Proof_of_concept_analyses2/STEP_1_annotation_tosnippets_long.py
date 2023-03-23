import os #basic foldering functions
from moviepy.video.io.ffmpeg_tools import ffmpeg_extract_subclip #this is the video clipping function
import pandas #used here for data reading (csv)
import numpy as np #for saving data in array format
from pydub import AudioSegment

#folder info processed data
processed_data = "H:/SiamangJaderpark_Processed/OpportunisticSampling_Song"
toprocess = './annotations/combined_full_scenes.csv'
outputfol = './snippets/'
outputfol2 = outputfol

#extract videos
#what video file .mov .mp4 .avi
vidformat = ".mxf"

it = 0
#also extract files once

andata = np.array(pandas.read_csv(toprocess)) #save into an numpy array (but first read.csv using pandas)     
for ch in range(0, andata.shape[0]): #now loop through each row of your annotions
    annot = andata[ch,0]
    ID = annot.replace('_zoom_synced_boom', '')
    session = ID.split('Session',1)[1][1]
    date = ID.split('_',3)[2] +'_'+ ID.split('_',3)[1] 
    videof = processed_data+'/'+date+'_processed/Session' + session + '/zoom_cam/' #get the video that has a matching code
    video = os.listdir(videof)
    begintime = andata[ch,1] #extract begintime
    endtime = andata[ch,2] #extract end time
    annotlabel = str(andata[ch,3])
    #now take the video and select the relevant chunk and write to the outputfolder
    ffmpeg_extract_subclip(videof+video[0], begintime/1000, (endtime/1000), targetname=outputfol+"_"+ID+"_" +str(annotlabel)+vidformat)
    #exract boom
    audioboomf = processed_data+'/'+date+'_processed/Session' + session + '/boom_mic/' #get the video that has a matching code
    audioboom = os.listdir(audioboomf)
    myaudio = AudioSegment.from_file(audioboomf+audioboom[0])
    snipped = myaudio[begintime:endtime]
    snipped.export(outputfol2+"_boommic_"+ID+"_" +str(annotlabel)+".wav", format='wav')
    #exract gopros combined
    audiof = processed_data+'/'+date+'_processed/Session' + session + '/gopros_mics/' #get the video that has a matching code
    audiogo= os.listdir(audiof)
    myaudio = AudioSegment.from_file(audiof+audiogo[0])
    snipped = myaudio[begintime:endtime]
    snipped.export(outputfol2+"_multisource_"+ID +"_" +str(annotlabel)+".wav", format='wav')