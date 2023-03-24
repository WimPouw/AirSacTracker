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
# comparison datasheet including both microphone versions (boom & multisource)

comparison <- readRDS('radius_acoustic_param_comparison_proof_Concept_1_corrected.rds')
# sound source: 4 go pros microphone
#comparison_data_new <- readRDS('radius_acoustic_param_comparison_proof_Concept_1.rds')

#sound source: boom? 
comparison_data_old <- readRDS('radius_acoustic_param_comparison_2.rds')
#comparison_data <- readRDS('radius_acoustic_param_comparison_scaled.rds')

# 02: data preparations ----

## radii > 400 are unlikely to be correct and are therefore deleted
comparison_data_fs_filtered <- comparison_numeric_full_fs %>% 
  dplyr::filter(radius < 400)



comparison_filtered <- comparison %>% 
  dplyr::filter(radius < 400)

# split into two microphone groupes
# grouped_mic_data <- comparison_boom_multisource_filtered %>%
#   group_by(mic, .add = TRUE) %>% 
#   group_split(mic)
# 
# 
# boom_mic <- grouped_mic_data[[1]]
# multisource <- grouped_mic_data[[2]]


## for correlation matrix we need a matrix only containing numeric variables

comparison_data_old_filtered_numeric <- comparison_data_Old_filtered %>% 
  ungroup() %>% 
  select(-audiofile, -frame, -sex) 


data_fs_filtered_numeric <- comparison_data_fs_filtered %>% 
  ungroup() %>% 
  select(-audiofile, -frame, -videofile, -sex, - ID, - index) 

# 03: correlation matrix ----

cor_all_old <- cor(comparison_data_old_filtered_numeric)
cor_all_Old_2 <-rcorr(as.matrix(comparison_data_old_filtered_numeric))

cor_all <- cor(data_fs_filtered_numeric)
cor_all_2 <-rcorr(as.matrix(data_fs_filtered_numeric))

#correlation plot, show all
corrplot(cor_all_2$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_2$r, type="lower", order="hclust", 
         p.mat = cor_all_2$P, sig.level = 0.05, insig = "blank", diag = FALSE)

# 04: individual plots ----

## plotting acoustic features over radius, color-coding for scene

## amplitude mean - sex

#scatter_ampl_mean <- comparison_data_filtered %>%
scatter_ampl_mean <- comparison_data_fs_filtered %>% 
  #dplyr::filter(sex == "m") %>% 
  ggplot(aes(y = ampl_mean, x = norm_radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Amplitude')+
  xlab('Airsac Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

#scatter_ampl_mean

## amplitude mean - ID

scatter_ampl_mean_ID <- comparison_data_fs_filtered %>%
  #dplyr::filter(sex == "m") %>% 
  #ggplot(aes(y = ampl_mean, x = radius, fill = ID, color = ID))+
  ggplot(aes(y = ampl_mean, x = norm_radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Amplitude')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

## entropy mean -sex

scatter_entropy_mean <- comparison_data_fs_filtered %>%
  ggplot(aes(y = entropy_mean, x = norm_radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Entropy')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## entropy mean - ID
scatter_entropy_mean_ID <- comparison_data_fs_filtered %>%
  #ggplot(aes(y = entropy_mean, x = radius, fill = ID, color = ID))+
  ggplot(aes(y = entropy_mean, x = norm_radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Mean Entropy')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## dominant frequency mean - sex

scatter_dom_mean <- comparison_data_fs_filtered %>%
  ggplot(aes(y = dom_mean, x = norm_radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Dominant Frequency mean [Hz]')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.x = element_blank())

## dominant frequency mean - ID

scatter_dom_mean_ID <- comparison_data_fs_filtered %>%
  ggplot(aes(y = dom_mean, x = norm_radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Dominant Frequency mean [Hz]')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20),
        axis.title.x = element_blank())

## spectral Centroid mean - sex

scatter_spec_Centroid <- comparison_data_fs_filtered %>%
  ggplot(aes(y = specCentroid_mean, x = norm_radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Centroid Mean')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

## spectral Centroid mean - ID

scatter_spec_Centroid_ID <- comparison_data_fs_filtered %>%
  ggplot(aes(y = specCentroid_mean, x = norm_radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Centroid Mean')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

## fundamental frequcny mean - sex

scatter_f0 <- comparison_data_fs_filtered %>%
  ggplot(aes(y = fundamental_mean, x = norm_radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Fundamental [Hz]')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))#,
        #axis.title.y = element_blank())

## fundamental frequcny mean - ID

scatter_f0_ID <- comparison_data_fs_filtered %>%
  ggplot(aes(y = fundamental_mean, x = norm_radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Fundamental [Hz]')+
  xlab('Airsac Radius [px]')+
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

#cowplot::plot_grid(scatter_ampl_mean_ID, scatter_entropy_mean_ID, nrow = 2)

ggsave("normalized_radius_01_acousticparam_per_ID.jpg", 
       dpi = 300, 
       width = 15,
       height = 10)
# 05: timeseries plot of inflation ----

# timeseries plots Fajar
plot_fajar <- comparison_data_filtered %>%
  mutate(frame = as.numeric(frame)) %>%
  filter(ID == 'Fajar') %>% 
  ggplot(aes(x = frame, y = radius, group = audiofile, color = audiofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Fajar')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Pelangi
plot_pelangi <- comparison_data_filtered %>%
  mutate(frame = as.numeric(frame)) %>%
  filter(ID == 'Pelangi') %>% 
  ggplot(aes(x = frame, y = radius, group = audiofile, color = audiofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Pelangi')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Baju
plot_baju <- comparison_data_filtered %>%
  mutate(frame = as.numeric(frame)) %>%
  filter(ID == 'Baju') %>% 
  ggplot(aes(x = frame, y = radius, group = audiofile, color = audiofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Baju')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Roger
plot_roger <- comparison_data_filtered %>%
  mutate(frame = as.numeric(frame)) %>%
  filter(ID == 'Roger') %>% 
  ggplot(aes(x = frame, y = radius, group = audiofile, color = audiofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Roger')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

cowplot::plot_grid(plot_fajar, plot_pelangi, plot_baju, plot_roger,
                   ncol = 2)
