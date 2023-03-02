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

# 01: load data I------------------------------

# set path to wav file source location
path <- choose.dir()
pattern <- "wav"
list_of_files_audio <- list.files(path = path, pattern = pattern)

# test data, only one wavfile
#wavfilelocation <- "D:/PostDoc/Donders/AirSacTracker_2/AirSacTracker/Proof_of_concept_analyses/snippets/audio/_Opp_June_15_Session_1_zoom_syncedboom_37_2_2_Pelangi.wav"

# 02: setting parameters ------------

fs_video <- 25 #sometimes 50, check how to make dynamic, I think it is noted in the tracking csv

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

subwave <- read_wav(wavfilelocation, time_exp = 1, from = audio_split_seq[wavsplit], to = (audio_split_seq[wavsplit]+(duration_frame*2)))

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

acoustic_param <- left_join(spec_param_all, formants_all, by = c('frame', 'audiofile'))
acoustic_param <- left_join(acoustic_param, fund_freq_mean_all, by = c('frame', 'audiofile'))


# 03: radius extraction ------------------

# 01: load data II ----------------
path <- choose.dir()
pattern <- "csv"
list_of_files_radius <- list.files(path = path, pattern = pattern)

# maximum radius of video snippets to correlate to formant spacing
radius_all <- data.frame(matrix(ncol = 3, nrow = 0))
z <- c("radius", "frame", "videofile")
colnames(radius_all) <- z


for (b in 1: length(list_of_files_radius)){
  
  radius <- data.frame()
  
  auto_data <- read_delim(paste(path, list_of_files_radius[b], sep = "\\"), delim = "," )
  
  # redefine column names to distinct names, delete first two rows with identifiers afterwards
  colnames <- c()
  colnames[1] <- "frames"
  for (x in 2: ncol(auto_data)){
    colnames[x] <- paste(auto_data[1,x], auto_data[2,x], sep = '_')
  }
  
  colnames(auto_data) <- colnames
  df_all <- auto_data[-1:-2,]
  
  threshold <- 0.6 #threshold for likelihood for DLC tracking of points
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("x", "y", "likelihood")
  colnames(df) <- x
  
  for (c in 1:nrow(df_all)){
    
    #subset DLC results to necessary columns  
    
    df_sub <- df_all %>% 
      select(Start_outline_outer_left_x,Start_outline_outer_left_y, Start_outline_outer_left_likelihood,
             Start_outline_outer_right_x, Start_outline_outer_right_y, Start_outline_outer_right_likelihood,
             LowestPoint_outline_x, LowestPoint_outline_y, LowestPoint_outline_likelihood,
             MidLowleft_outline_x, MidLowleft_outline_y, MidLowleft_outline_likelihood,
             MidLowright_outline_x, MidLowright_outline_y, MidLowright_outline_likelihood)
    
    df_sub <- as.data.frame(apply(df_sub,2,as.numeric))
    
    # save DLC points in necessary format to estimate circles
    
    for (d in 1:5){
      end_col <- d*3
      start_col <- end_col - 2
      
      df[d, 1:3] <-  df_sub[c, start_col:end_col]
    } 
    
    # check for likelihood
    # do we need to make sure, to somehow save information, which points are used? 
    
    df_filter <- df %>% 
      dplyr::filter(likelihood >= threshold)
  
  if(nrow(df_filter)>=3){
  circles_LAN <- CircleFitByLandau(df_filter[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
  } else {
    circles_LAN <- c(NA, NA, NA)
  }  
  circles_res <- c(circles_LAN[3],c, list_of_files_radius[b])
  
  radius <- rbind(radius, circles_res)
  colnames(radius) <- z
}

  radius_all <- rbind(radius_all, radius)
  
}

# interpolate radii

#radius_all_inter <- radius_all %>% 
#  mutate (radius_inter <- na.approx(radius_all$radius, maxgap = 2))


# 04: combine radius and acoustic parameters-----------------


filename <- strsplit(radius_all$videofile, split  = "DLC")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_all <- radius_all %>% 
  mutate(audiofile = filename[,1],
         radius = as.numeric(radius)) %>% 
  select(radius, frame, audiofile)


# only run the following 2 lines once per analysis!
#acoustic_param <- acoustic_param %>% 
#  mutate(audiofile = substring(audiofile, 1, nchar(audiofile)-4))

comparison <- left_join(radius_all, acoustic_param, by = c('frame', 'audiofile'))

# correlation matrix
comparison_numeric_full <- comparison %>% 
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
         ampl_mean_noSilence = as.numeric(ampl_mean_noSilence))

#saving comparison for use in other scripts, etc.
# full dataframe, including frame and videofile/audiofile name
saveRDS(comparison_numeric_full,'radius_acoustic_param_comparison_2.rds')



