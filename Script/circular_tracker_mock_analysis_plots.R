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
library(cowplot)

# 01: load data ----

comparison_data <- readRDS('radius_acoustic_param_comparison_2.rds')

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

# names: fajar, baju, roger, pelangi
for (b in 1:nrow(comparison_data_filtered)){
  
 if(grepl('Pelangi', comparison_data_filtered$audiofile[b], , ignore.case = TRUE) == TRUE){
   
   comparison_data_filtered$ID[b] = 'Pelangi'
   
 }  else if (grepl('Fajar',comparison_data_filtered$audiofile[b], ignore.case = TRUE) == TRUE){
   
   comparison_data_filtered$ID[b] = 'Fajar'
   
 } else if (grepl('Baju',comparison_data_filtered$audiofile[b], ignore.case = TRUE) == TRUE){
   
   comparison_data_filtered$ID[b] = 'Baju'
   
 } else if (grepl('Roger',comparison_data_filtered$audiofile[b], ignore.case = TRUE) == TRUE){
   
   comparison_data_filtered$ID[b] = 'Roger'
   
 } else {
   
   comparison_data_filtered$ID[b] = 'NA'
 }
  
}


# 03: correlation matrix ----

cor_all <- cor(comparison_data_filtered_numeric)
cor_all_2<-rcorr(as.matrix(comparison_data_filtered_numeric))


#correlation plot, show all
corrplot(cor_all_2$r, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_2$r, type="upper", order="original", 
         p.mat = cor_all_2$P, sig.level = 0.05, insig = "blank", diag = FALSE)

# 04: individual plots ----

## plotting acoustic features over radius, color-coding for scene

## amplitude mean - sex

scatter_ampl_mean <- comparison_data_filtered %>%
  ggplot(aes(x = ampl_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Mean Amplitude')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## amplitude mean - ID

scatter_ampl_mean_ID <- comparison_data_filtered %>%
  ggplot(aes(x = ampl_mean, y = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Mean Amplitude')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## entropy mean -sex

scatter_entropy_mean <- comparison_data_filtered %>%
  ggplot(aes(x = entropy_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Mean Entropy')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## entropy mean - ID
scatter_entropy_mean_ID <- comparison_data_filtered %>%
  ggplot(aes(x = entropy_mean, y = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Mean Entropy')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## dominant frequency mean - sex

scatter_dom_mean <- comparison_data_filtered %>%
  ggplot(aes(x = dom_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Dominant Frequency mean [Hz]')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())

## dominant frequency mean - ID

scatter_dom_mean_ID <- comparison_data_filtered %>%
  ggplot(aes(x = dom_mean, y = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Dominant Frequency mean [Hz]')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())

## spectral Centroid mean - sex

scatter_spec_Centroid <- comparison_data_filtered %>%
  ggplot(aes(x = specCentroid_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Spectral Centroid Mean')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())

## spectral Centroid mean - ID

scatter_spec_Centroid_ID <- comparison_data_filtered %>%
  ggplot(aes(x = specCentroid_mean, y = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Spectral Centroid Mean')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())

## fundamental frequcny mean - sex

scatter_f0 <- comparison_data_filtered %>%
  ggplot(aes(x = fundamental_mean, y = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Fundamental [Hz]')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))#,
        #axis.title.y = element_blank())

## fundamental frequcny mean - ID

scatter_f0_ID <- comparison_data_filtered %>%
  ggplot(aes(x = fundamental_mean, y = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab('Fundamental [Hz]')+
  ylab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))#,
#axis.title.y = element_blank())




## cowplot

cowplot::plot_grid(scatter_ampl_mean, scatter_dom_mean, scatter_entropy_mean, scatter_spec_Centroid,
                   rel_widths = c(0.45,0.55), ncol = 2)
                   #labels = c("A", "B", "C", "D"), ncol = 2)

cowplot::plot_grid(scatter_ampl_mean_ID, scatter_dom_mean_ID, scatter_entropy_mean_ID, scatter_spec_Centroid_ID,
                   rel_widths = c(0.45,0.55), ncol = 2)
#labels = c("A", "B", "C", "D"), ncol = 2)

ggsave("acoustic_parameter_radius_ID.jpg", 
       dpi = 300, 
       width = 25,
       height = 16)
