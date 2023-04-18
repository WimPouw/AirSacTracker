# # Acoustic Analysis and Comparison to AirSac Radius
# How does airsac inflation of boom relate to different acoustic parameters of
# accompanying bark in great call sequences of female siamangs
# author: Lara S. Burchardt
# version: 1.0, 24.03.2023

# Part I: Analysis

# Input:

# 1) acoustic parameters of barks from script: circular_tracker_mock_analysis
# 2) radii of previous boom call from script: circular_estimations_DLC_for_acoustics

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

acoustics_bark_proof_2 <- readRDS("acoustic_features_bark_great_call_proof_2.rds")

radius_boom_proof_2 <- readRDS("proof_2_boom_DLC_estimated_radii_normalized_new.rds")

# 02: setting parameters ------------



# 03: acoustic parameters --------------------

# we need a column with the date info and a column with the id info, then match according to those two

# were in the name is the boom/bark? can we just delete that and then match instead of the other way around?

# yes! example name: _multisource_Opp_August_14_Session_1_zoom_syncedboom_1_1_bark
# _Opp_August_14_Session_1_zoom_syncedboom_1_1_boomDLC_resnet101_Deep_AirSacTrackingV1Jan1shuffle1_500000_labeled

# we split at "bark" for audio and at "boomDLC" for radii this is combined in a completely new column called "match"
# we leave the audiofilename in there and the videofilename, so that we can track errors if necessary

# video match name
filename <- strsplit(radius_boom_proof_2$videofile, split  = "boomDLC")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_boom_proof_2$match <- filename[,1]

# audio match name

filename <- strsplit(acoustics_bark_proof_2$audiofile, split  = "bark")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))
filename <- strsplit(filename$V1, split  = "multisource")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

acoustics_bark_proof_2$match <- filename[,2]

# radius data preparation
## we only want to compare the last radius of each boom with the averaged amplitude or ma amplitude of the following bark

radius_boom_last <- radius_boom_proof_2 %>% 
  group_split(videofile)

for(z in 1:length(radius_boom_last)){

  data <- radius_boom_last[[z]]
  data <- data %>% 
    mutate(radius = as.numeric(radius),
           frame = as.numeric(frame))
  data <- data[!is.na(data$radius),]
  
  if(nrow(data) > 1){
  data_max[z,] <- data[data$frame == max(data$frame),]
  } else {
    data_max[z,] <- NA
  }
}

# amplitude data preparation
# we want to compare to maximum and or average of the acoustic parameters of the bark

acoustic_bark_summary <- acoustics_bark_proof_2 %>%
  mutate(ampl_mean = as.numeric(ampl_mean),
         specCentroid_mean = as.numeric(specCentroid_mean),
         dom_mean = as.numeric(dom_mean),
         entropy_mean = as.numeric(entropy_mean),
         f1_freq_mean = as.numeric(f1_freq_mean),
         f2_freq_mean = as.numeric(f2_freq_mean),
         formant_spacing = as.numeric(formant_spacing),
         #fundamental_mean = as.numeric(fundamental_mean)*1000,      #fundamental frequency is given in kHz, transform to Hertz
         duration_noSilence = as.numeric(duration_noSilence),
         ampl_mean_noSilence = as.numeric(ampl_mean_noSilence)) %>% 
  group_by(audiofile) %>% 
  summarise_at(c("ampl_mean", "specCentroid_mean", "dom_mean"), funs(mean, max), na.rm = TRUE)

filename <- strsplit(acoustic_bark_summary$audiofile, split  = "bark")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))
filename <- strsplit(filename$V1, split  = "multisource")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

acoustic_bark_summary$match <- filename[,2]


# joining together
# check again, there is a lot of NAs, why?

proof_2_radius_acoustics <- left_join(data_max, acoustic_bark_summary)

proof_2_radius_acoustics %>% 
  ggplot(aes(x= radius, y = ampl_mean_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Amplitude of subsequent bark')+
  xlab('Maximum Airsac Radius of previous Boom [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)
        #axis.title.x = element_blank())

proof_2_radius_acoustics %>% 
  ggplot(aes(x= radius, y = ampl_mean_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Maximum Amplitude of subsequent bark')+
  xlab('Maximum Airsac Radius of previous Boom [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)


proof_2_radius_acoustics %>% 
  ggplot(aes(x= radius, y = dom_mean_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Dominant Frequency of subsequent bark')+
  xlab('Maximum Airsac Radius of previous Boom [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)


# proof_2_radius_acoustics_filter <- proof_2_radius_acoustics %>% 
#   dplyr::filter(radius < 400)
# 
# proof_2_radius_acoustics_filter_numeric <- proof_2_radius_acoustics_filter %>% 
#   ungroup() %>% 
#   select(-audiofile, -frame, -videofile, -match) %>% 
#   mutate(ampl_mean = as.numeric(ampl_mean),
#          specCentroid_mean = as.numeric(specCentroid_mean),
#          dom_mean = as.numeric(dom_mean),
#          entropy_mean = as.numeric(entropy_mean),
#          f1_freq_mean = as.numeric(f1_freq_mean),
#          f2_freq_mean = as.numeric(f2_freq_mean),
#          formant_spacing = as.numeric(formant_spacing),
#          #fundamental_mean = as.numeric(fundamental_mean)*1000,      #fundamental frequency is given in kHz, transform to Hertz
#          duration_noSilence = as.numeric(duration_noSilence),
#          ampl_mean_noSilence = as.numeric(ampl_mean_noSilence),
#          radius = as.numeric(radius))

# 03: correlation matrix ----

cor_all_proof2 <- cor(proof_2_radius_acoustics)
cor_all_proof2_2 <-rcorr(as.matrix(proof_2_radius_acoustics))


#correlation plot, show all
corrplot(cor_all_proof2_2$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_proof2_2$r, type="lower", order="hclust", 
         p.mat = cor_all_proof2_2$P, sig.level = 0.05, insig = "blank", diag = FALSE)
