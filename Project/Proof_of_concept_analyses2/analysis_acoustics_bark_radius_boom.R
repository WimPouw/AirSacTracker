# analysis: how do acoustic features of bark relate to the radius of a previous boom call
# in female great call sequences?

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: dataset with acoustic features and dataset with radius information
# end point: matched dataset of acoustics and radii + analysis and visualization

############

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse","conicfit", "scales", "spiro", "signal", "foreach", "Hmisc", "cowplot", "corrplot")
library("foreach")

# 01: data ----
## 01a: load data ----

acoustics_bark <- readRDS("acoustics_bark_multi_checked_proof_2.rds")

radius_boom_bark_sequences <- readRDS("circle_estimation_boom_DLC_proof2.rds")

## 01b: data preparation ----

### video match name
filename <- strsplit(radius_boom_bark_sequences$videofile, split  = "boomDLC")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_boom_bark_sequences$match <- filename[,1]

### radius data preparation: we only want to compare the last radius of each boom with the averaged amplitude or max amplitude of the following bark
### last available radius of each boom video is saved in data_max (one row per video)


radius_boom_bark_sequences_last <- radius_boom_bark_sequences %>% 
  group_split(videofile)

for(z in 1:length(radius_boom_bark_sequences_last)){
  
  data <- radius_boom_bark_sequences_last[[z]]
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

### amplitude data preparation: we want to compare to maximum and or average of the acoustic parameters of the bark

#acoustic_bark_summary <- acoustics_bark %>%
#  group_by(audiofile) %>% 
#  summarise_at(c("ampl", "specCentroid", "dom", "entropy", "entropySh", "f1_freq", "f2_freq", "specSlope", "pitch", "harmEnergy",
#                 "peakFreq", "fmPurity", "HNR"), list(mean = mean, max = max, min = min), na.rm = TRUE)

acoustic_bark_summary <- acoustics_bark %>%
  group_by(audiofile) %>% 
  summarise_at(c("ampl", "specCentroid", "dom", "entropy", "entropySh", "f1_freq", "f2_freq", "specSlope", "pitch", "harmEnergy",
                 "peakFreq", "fmPurity", "HNR"), list(~mean(.x,na.rm = TRUE), ~max(.x,na.rm = TRUE), ~min(.x,na.rm = TRUE)))

### audio match name

filename <- strsplit(acoustic_bark_summary$audiofile, split  = "bark")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))
filename <- strsplit(filename$V1, split  = "multisource")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

acoustic_bark_summary$match <- filename[,2]

### data joining: match acoustics_bark with last radius of preceeding boom

combined_radius_acoustics_boom_bark <- left_join(data_max, acoustic_bark_summary)
#combined_radius_acoustics <- left_join(combined_radius_acoustics, duration_bark)

# 02: data visualization ----

plot1<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>% 
  ggplot(aes(x= radius, y = ampl_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Amplitude subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)
#axis.title.x = element_blank())

plot2<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = ampl_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max Amplitude subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)


plot3<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = dom_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Dom Freq subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot4<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = dom_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max Dom Freq subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot5<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = specCentroid_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean spec Cent subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot6<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = specCentroid_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max spec Cent subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot7<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = entropy_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Entropy (Wiener) subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot8<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = entropy_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max Entropy (Wiener) subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot9<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = entropySh_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Entropy (Shannon) subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot10<- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = entropySh_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max Entropy (Shannon) subsequent bark')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot11 <- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = f1_freq_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean F1 subsequent bark [Hz]')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',)

plot12 <- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius>80) %>%
  ggplot(aes(x= radius, y = f2_freq_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean F2 subsequent bark [Hz]')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

plot13 <- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius > 80) %>%
  dplyr::filter(pitch_mean > 300) %>% #to filter out part of snippets that might involve boom sound
  ggplot(aes(x= radius, y = pitch_mean))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean F0 subsequent bark [Hz]')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

plot14 <- combined_radius_acoustics_boom_bark %>%
  dplyr::filter(radius > 80) %>%
  dplyr::filter(pitch_max > 300) %>% #to filter out part of snippets that might involve boom sound
  ggplot(aes(x= radius, y = pitch_max))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Max F0 subsequent bark [Hz]')+
  #xlab('Max Radius previous Boom [px]')+
  xlab('Max Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

cowplot::plot_grid(plot1,plot5,plot7,plot13, plot2, plot6, plot8,plot14, ncol = 4, nrow = 2)

# 03: correlation matrix ----

acoustic_correlation <- combined_radius_acoustics_boom_bark %>% 
  select(radius, ampl_mean, pitch_mean, entropy_mean, specCentroid_mean, f1_freq_mean, f2_freq_mean,
         harmEnergy_mean, peakFreq_mean, HNR_mean)#,
         #ampl_max, pitch_max, entropy_max, specCentroid_max, f1_freq_max, f2_freq_max,
         #harmEnergy_max, peakFreq_max, HNR_max,
         #ampl_min, pitch_min, entropy_min, specCentroid_min, f1_freq_min, f2_freq_min,
         #harmEnergy_min, peakFreq_min, HNR_min)


  #dplyr::filter(radius>80) %>% 
  #select(-videofile, -audiofile, -match)

cor_all_proof2 <- cor(acoustic_correlation)
cor_all_proof2_2 <-rcorr(as.matrix(acoustic_correlation))


#correlation plot, show all
corrplot(cor_all_proof2_2$r[1,2:19, drop = FALSE], type="upper", 
         tl.col = "black")

#correlation plot, only show significant correlations
corrplot(cor_all_proof2_2$r[1,2:19, drop = FALSE], type="upper", 
         p.mat = cor_all_proof2_2$P[1,2:19, drop = FALSE], sig.level = 0.05,
         tl.col = "black",
         insig = "blank")

