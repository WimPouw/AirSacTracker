# Acoustic Analysis and Comparison to AirSac Radius
# How does airsac inflation relate to different acoustic parameters of
# accompanying boom call
# author: Lara S. Burchardt
# version: 1.0, 27.01.2023
# Part II: Plotting

# 00: load packages ----

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

# 01: load data ----

comparison_data <- readRDS('radius_acoustic_param_comparison.rds')

# 02: data preparations ----

## radii > 400 are unlikely to be correct and are therefore deleted

comparison_data_filtered <- comparison_data %>% 
  filter(radius < 400)

## for correlation matrix we need a matrix only containing numeric variables

comparison_data_filtered_numeric <- comparison_data_filtered %>% 
  select(-audiofile, -frame) 

## for comparison we would like to add the sex of the siamang as another column
## the only female is Pelangi

index_female <- grep('Pelangi', comparison_data_filtered$audiofile)
comparison_data_filtered$index <- 1:nrow(comparison_data_filtered)

for(a in 1: nrow(comparison_data_filtered)){
if(comparison_data_filtered$index[a] %in% index_female == TRUE){

  comparison_data_filtered$sex[a] = 'f'
  
} else {
  
  comparison_data_filtered$sex[a] = 'm'
  
}
}

# 03: correlation matrix ----

cor_all <- cor(comparison_data_filtered_numeric)
cor_all_2<-rcorr(as.matrix(comparison_data_filtered_numeric))

#correlation plot, show all
corrplot(cor_all_2$r, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_2$r, type="lower", order="original", 
         p.mat = cor_all_2$P, sig.level = 0.05, insig = "blank", diag = FALSE)

# 04: individual plots ----

## plotting acoustic features over radius, color-coding for scene

## amplitude mean

comparison_data_filtered %>%
  ggplot(aes(x = ampl_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  #geom_smooth(method = 'lm')+
  xlab('Mean Amplitude')+
  ylab('Airsac Radius [px]')
  #theme(legend.position = 'none')


