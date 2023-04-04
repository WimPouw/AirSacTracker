# Acoustic Analysis and Comparison to AirSac Radius
# How does airsac inflation relate to different acoustic parameters of
# accompanying boom call
# author: Lara S. Burchardt
# version: 1.0, 27.01.2023

# Part I: Analysis

# To do -------------------------

# na.approximation?
# normalization nose bridge points
# check more acoustic features
# plots: nice plots with color code for different scenes
# clean code, add comments 
# add to github
# add saving of results
# split into different files? with saving in between?
# match enumerations, edit headings etc. 

# 00: load packages --------------------------

library(tidyverse)
library(phonTools)
library(rPraat)
library(bioacoustics)
library(tuneR)    #
library(soundgen) # to analyze spectral parameters of audio files
library(conicfit) # circle estimation
library(psych)    # correlations
library(zoo)      # na approximation
library(corrplot)
library(Hmisc)
library(seewave)  # fundamental frequency estimation

# 00b: functions -----


high_pass <- function(x, samplingrate, order, highpasscutoff)
{bf <- butter(order,highpasscutoff/samplingrate, type="high") #normalized frequency
x <<- as.numeric(signal::filtfilt(bf, x))} 

# 01: load data I------------------------------

# set path to wav file source location
path <- choose.dir()
pattern <- "wav"
list_of_files_audio <- list.files(path = path, pattern = pattern)

# test data, only one wavfile
#wavfilelocation <- "D:/PostDoc/Donders/AirSacTracker_2/AirSacTracker/Proof_of_concept_analyses/snippets/audio/_Opp_June_15_Session_1_zoom_syncedboom_37_2_2_Pelangi.wav"

# 02: setting parameters ------------



# 03: acoustic parameters --------------------

# video fs = 25 (or 50), 1 frame = 0.04 sec (or 0.02sec)
# audio fs = 48000
# 1920 audio samples are 0.04 sec (960 are 0.02)

# splitting wav file: https://stackoverflow.com/questions/39047687/fast-way-for-splitting-large-wav-file-using-r

# we want formants per video frame as well as other acoustic parameters
# fetaures we want: formant spacing, fundamental frequency 


# 02a: acoustic analysis ---------------

formants_all <- data.frame(matrix(ncol = 3, nrow = 0))
spec_param_all <- data.frame(matrix(ncol = 10, nrow = 0))
fund_freq_mean_all <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("formant_spacing", "frame", "audiofile")
y <- c("ampl_mean","specCentroid_mean","dom_mean","entropy_mean",
       "f1_freq_mean","f2_freq_mean", "duration_noSilence", "ampl_mean_noSilence",
       "frame", "audiofile")
z <- c("fundamental_mean", "frame", "audiofile")
colnames(formants_all) <- x
colnames(spec_param_all) <- y
colnames(fund_freq_mean_all) <- z


for (a in 1: length(list_of_files_audio)){
  
# include if condition for fs video, all August videos have 50, all June videos have 25  
  if(grepl('June', list_of_files_audio[a], ignore.case = TRUE) == TRUE){
    
    fs_video = 25
    frame_factor = 2
      
  }  else if (grepl('August', list_of_files_audio[a], ignore.case = TRUE) == TRUE){
    
    fs_video = 50
    frame_factor = 4
  } 

#      
formants <- data.frame(matrix(ncol = 3, nrow = 0))
spec_param <- data.frame(matrix(ncol = 8, nrow = 0))
fund_freq_mean <- data.frame(matrix(ncol = 3, nrow = 0))
  
wavfilelocation <- paste(path, list_of_files_audio[a], sep = "/")

wav <- read_wav(wavfilelocation, time_exp = 1, from = NULL, to = NULL)

fs_audio <- wav@samp.rate
samples <- length(wav@left)
duration_audio <- length(wav@left)/wav@samp.rate



# audio smaples per video frame, to cut wav file into snipptes matching video frames
audio_samples_per_frame <- fs_audio * (1/fs_video)
duration_frame <- 1/fs_video

audio_split_seq <- seq(from = 0, to = duration_audio, by = duration_frame)

number_of_subwavs <- as.integer(duration_audio/duration_frame)

# in loop, supply sequence of start and end times based on duration_frame and audio duration
# "splitting audio", to extract acoustic features per snippet matching a video frame
for (wavsplit in 1:number_of_subwavs){

subwave <- read_wav(wavfilelocation, time_exp = 1, from = audio_split_seq[wavsplit], to = (audio_split_seq[wavsplit]+(duration_frame*frame_factor))) # was *2 for video fs 25, changed to *4 for video fs 50, to have same duration?

# only use for bark analysis, comment out otherwise 
subwave <- high_pass(subwave, samplingrate = fs_video, order = 2, highpasscutoff = 300)
  

## formant spacing ----------------

spectrum <- phonTools::findformants(subwave@left, fs =subwave@samp.rate,showrejected = FALSE,verify = FALSE)
formant_spacing <- mean(diff(spectrum$formant[1:5])) 

formant_spacing_res <- c(formant_spacing, wavsplit, list_of_files_audio[a])

formants <- rbind(formants, formant_spacing_res)

colnames(formants) <- x

## fundamental frequency

fund_freq <- fund(subwave, wl = 512, ovlp = 50, fmax = 1000) #if no fmax is defined, no fundamental is found, we know the boom is around ~200 Hz, so 1000 is still a very conservative limit
fund_freq <- as.data.frame(fund_freq)

fund_freq_mean_loop <- c(mean(fund_freq$y, na.rm = TRUE), wavsplit, list_of_files_audio[a]) #as we are using short audio snippets to match to video frames, we take the mean per video frame

fund_freq_mean <- rbind(fund_freq_mean, fund_freq_mean_loop)
colnames(fund_freq_mean) <- z
## other spectral parameters------------

#spec_param_list <- analyze(subwave, roughness = list(windowLength = 15, step = 3, amRes = 100))
spec_param_list <- analyze(subwave,
                           loudness = NULL,      # no loudness analysis
                           novelty = NULL,       # no novelty analysis
                           roughness = NULL,     # no roughness analysis
                            )

spec_param_vector <- c(spec_param_list$summary$ampl_mean, spec_param_list$summary$specCentroid_mean, 
                       spec_param_list$summary$dom_mean, spec_param_list$summary$entropy_mean,
                       spec_param_list$summary$f1_freq_mean, spec_param_list$summary$f2_freq_mean,
                       spec_param_list$summary$duration_noSilence, spec_param_list$summary$ampl_noSilence_mean,
                       wavsplit, list_of_files_audio[a] )

spec_param <- rbind(spec_param, spec_param_vector)

colnames(spec_param) <- y

} 

formants_all <- rbind(formants_all, formants)
spec_param_all <- rbind(spec_param_all, spec_param)
fund_freq_mean_all <- rbind(fund_freq_mean_all, fund_freq_mean)

}

#acoustic_param <- spec_param_all
acoustic_param <- left_join(spec_param_all, formants_all, by = c('frame', 'audiofile'))
acoustic_param <- left_join(acoustic_param, fund_freq_mean_all, by = c('frame', 'audiofile'))

acoustic_param_backup <- acoustic_param
# 03: combine radius and acoustic parameters-----------------

radius_results <- readRDS('DLC_estimated_radii_meta_data_norm_scaled_2.rds')

# prepare name for joining for radius results
filename <- strsplit(radius_results$videofile, split  = "DLC")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_results$audiofile <- filename[,1]

# prepare name & micorphone ID column for joining for acoustic results

for (b in 1:nrow(acoustic_param)){
  
  if(grepl('boommic', acoustic_param$audiofile[b], ignore.case = TRUE) == TRUE){
    
    acoustic_param$mic[b] = 'boommic'
    
  }  else if (grepl('multisource', acoustic_param$audiofile[b], ignore.case = TRUE) == TRUE){
    
    acoustic_param$mic[b] = 'multisource'
    
  }
}


filename_audio <- strsplit(acoustic_param$audiofile, split = "multisource")
filename_audio <- as.data.frame(matrix(unlist(filename_audio),ncol=2,byrow=T))
filename_audio <- strsplit(filename_audio$V2, split = "boommic")
filename_audio <- as.data.frame(matrix(unlist(filename_audio),ncol=2,byrow=T))

acoustic_param$audiofile <- filename_audio[,2]

# only run the following 2 lines once per analysis!
#acoustic_param <- acoustic_param %>% 
#  mutate(audiofile = substring(audiofile, 1, nchar(audiofile)-4))

comparison_fs <- plyr::join(radius_results, acoustic_param, by= c("audiofile", "frame"), type="left", match="first")

# correlation matrix
comparison_numeric_full_fs <- comparison_fs %>% 
#comparison_numeric <- comparison %>% 
  #select(-audiofile, -frame) %>% 
  mutate(ampl_mean = as.numeric(ampl_mean),
         specCentroid_mean = as.numeric(specCentroid_mean),
         dom_mean = as.numeric(dom_mean),
         entropy_mean = as.numeric(entropy_mean),
         f1_freq_mean = as.numeric(f1_freq_mean),
         f2_freq_mean = as.numeric(f2_freq_mean),
         formant_spacing = as.numeric(formant_spacing),
         fundamental_mean = as.numeric(fundamental_mean)*1000,      #fundamental frequency is given in kHz, transform to Hertz
         duration_noSilence = as.numeric(duration_noSilence),
         ampl_mean_noSilence = as.numeric(ampl_mean_noSilence),
         radius = as.numeric(radius))

#saving comparison for use in other scripts, etc.
# full dataframe, including frame and videofile/audiofile name
saveRDS(comparison_numeric_full_fs,'radius_acoustic_param_comparison_proof_Concept_1_fs.rds')



