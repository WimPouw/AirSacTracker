# functions for acoustic analysis for different subsequent analyses
# 

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: folder with audio snippets
# end point: dataset of acoustic features for the chosen audio snippets for subsequent
# comparison to airsac radii

############

# helpful links:
# http://cran.nexr.com/web/packages/soundgen/vignettes/acoustic_analysis.html

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse", "signal", "phonTools", "rPraat", "bioacoustics", "tuneR",
             "soundgen", "psych", "zoo", "corrplot", "Hmisc", "seewave", "plyr")

# 01a: load data ----

path <- choose.dir()
pattern <- "wav"

# 02: functions -----

acoustic_analysis <- function(path, pattern){

list_of_files_audio <- list.files(path = path, pattern = pattern)
  
  
## sub function a ----

# fundamental frequency and other parameters are calculated with different functions
# and therefore need to be saved separately at first 
spec_param_all <- data.frame()
fund_freq_mean_all <- data.frame(matrix(ncol = 3, nrow = 0))

z <- c("fundamental_mean", "frame", "audiofile")

#colnames(spec_param_all) <- y
colnames(fund_freq_mean_all) <- z

for (a in 1: length(list_of_files_audio)){
  
  #spec_param <- data.frame()
  fund_freq_mean <- data.frame(matrix(ncol = 3, nrow = 0))
  
  wavfilelocation <- paste(path, list_of_files_audio[a], sep = "/")
  
  wav <- read_wav(wavfilelocation, time_exp = 1, from = NULL, to = NULL)
  
  fs_audio <- wav@samp.rate
  samples <- length(wav@left)
  duration_audio <- length(wav@left)/wav@samp.rate
  
  
  # audio smaples per video frame, to cut wav file into snipptes matching video frames
  audio_samples_per_frame <- fs_audio * (1/50)  #50 is fps that we use as reference here
  # for other situations, we use the dynamic fs_video to differentiate between the 25 and 50 fps videos
    
    ## fundamental frequency
 #!!   # windowlength should match videoframe duration, check notes on how we wanted to achieve that
    
    fund_freq <- fund(wav, wl = audio_samples_per_frame, ovlp = 50, fmax = 1000, plot = FALSE) #if no fmax is defined, no fundamental is found, we know the boom is around ~200 Hz, so 1000 is still a very conservative limit
    fund_freq <- as.data.frame(fund_freq)
    
    #fund_freq_mean_loop <- c(mean(fund_freq$y, na.rm = TRUE), wavsplit, list_of_files_audio[a]) #as we are using short audio snippets to match to video frames, we take the mean per video frame
    
    fund_freq_mean <- rbind(fund_freq_mean, fund_freq_mean_loop)
    colnames(fund_freq_mean) <- z
    ## other spectral parameters------------
    
    #spec_param_list <- analyze(subwave, roughness = list(windowLength = 15, step = 3, amRes = 100))
    spec_param_list <- analyze(wav,
                               loudness = NULL,      # no loudness analysis
                               novelty = NULL,       # no novelty analysis
                               roughness = NULL,     # no roughness analysis
                               windowLength = 20,  # windowLength in msec, matched to videoframe length of 50fps video
                               step = 10          # step size of sliding window, defining overlap, for wl 0.02 and step 0.01 we have an 50% overlap
    )
    
    spec_param_list <- spec_param_list$detailed
    
    spec_param_list$audiofile <- list_of_files_audio[a]
    
    spec_param_all <- rbind(spec_param_all, spec_param_list)
    fund_freq_mean_all <- rbind(fund_freq_mean_all, fund_freq_mean)
    
  } 
# to reduce size, we delete all columns that only have NA
# after loop, because different files have different NA columns, so we first need to leave all columns in
not_any_na <- function(x) any(!is.na(x))

spec_param_all <- spec_param_all %>%
  select(where(not_any_na))

# we define a "frame" column, to be matched to fundamental frequency and video frames, this is a proxy! This is not exactly a frame,
# but a reasonable way to match correct snippets with each other
spec_param_all <- spec_param_all %>% 
  mutate(frame = time/10)

acoustic_param <- plyr::join(spec_param_all, fund_freq_mean_all, by= c("audiofile", "frame"), type="left", match="first")

acoustic_param_backup <- acoustic_param

return(datasheet_with_acoustic_features_per_audio_and_frame)
}
