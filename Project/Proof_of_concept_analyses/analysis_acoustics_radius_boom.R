# analysis: how do acoustic features of boom call relate to its radius?

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

install_load("tidyverse","conicfit","corrplot", "scales", "spiro", "signal", "foreach", "Hmisc")
library("foreach")

# 01: data ----

## 01a: load data ----
acoustics_boom <- readRDS("acoustic_analysis_boom_multi_checked.rds")

radius_boom <- readRDS("radius_estimation_boom_proof1.rds")

## 01b: data preparation ----

### get video sampling rate from video name (in June recordings were made with 25 fps, in August with 50)

for (a in 1: nrow(acoustics_boom)){
  
  if(grepl('June', acoustics_boom$audiofile[a], ignore.case = TRUE) == TRUE){
    
    acoustics_boom$fs_video[a] <- 25
    
  }  else if (grepl('August', acoustics_boom$audiofile[a], ignore.case = TRUE) == TRUE){
    
    acoustics_boom$fs_video[a] <- 50
    
  } 
}

## add original video sampling rate to radius_boom Data as well  
for (a in 1: nrow(radius_boom)){
  
  if(grepl('June', radius_boom$videofile[a], ignore.case = TRUE) == TRUE){
    
    radius_boom$fs_video[a] <- 25
    
  }  else if (grepl('August', radius_boom$videofile[a], ignore.case = TRUE) == TRUE){
    
    radius_boom$fs_video[a] <- 50
    
  } 
}

radius_boom <- radius_boom %>% 
  mutate(frame = as.numeric(frame))

subset_radius_fps50 <- radius_boom %>% 
  dplyr::filter(fs_video == 50)

subset_radius_fps25 <- radius_boom %>% 
  dplyr::filter(fs_video == 25) 


### prepare match names to match audio to radius snippets

## 01b_2: data prep version 2 ----
## downsampling approach: instead of averaging acoustic data for every 2 rows to match 25fps, we now run acoustic analysis with 25 fps and
# downsample 50fps video data
# so we need to downsample the radius data, we do that by taking every second row of the dataset (per video) and then need a new index as new frame
# variable

downsampled_50fps_radius <- subset_radius_fps50 %>%
  group_by(videofile) %>%
  dplyr::slice(seq(1, n(), by = 2))

# add new frame index
downsampled_50fps_radius <- downsampled_50fps_radius %>% 
  group_by(videofile) %>%
  mutate( original_frame = frame,
          frame = row_number(videofile))

# combine transformed 50fps with 25fps again

radius_boom_downsampled <- rbind(downsampled_50fps_radius, subset_radius_fps25)

# video match name
filename <- strsplit(radius_boom_downsampled$videofile, split  = "DLC_resnet")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_boom_downsampled$match <- filename[,1]

# audio match name

filename <- strsplit(acoustics_boom$audiofile, split  = "multisource")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))
acoustics_boom$match <- filename[,2]
acoustics_boom$match <- str_sub(acoustics_boom$match, 1, -5)

### matching radius and acoustics dataframe by match name and frame

comparison_radius_boom <- left_join(radius_boom_downsampled,acoustics_boom,  by= c("match", "frame"))

# filter radii that are definitly wrong, i.e. very large
# we set the maximum to the approx. maximum we found when tracking manually

comparison_radius_boom <- comparison_radius_boom %>%
  mutate(radius = as.numeric(radius)) %>% 
  dplyr::filter(radius < 270)


# 02: analysis ----

# 02a: adding metadata ----

# adding column with information on sex
index_female <- grep('Pelangi', comparison_radius_boom$videofile)
comparison_radius_boom$index <- 1:nrow(comparison_radius_boom)

for(a in 1: nrow(comparison_radius_boom)){
  if(comparison_radius_boom$index[a] %in% index_female == TRUE){
    
    comparison_radius_boom$sex[a] = 'f'
    
  } else {
    
    comparison_radius_boom$sex[a] = 'm'
    
  }
}

# names: fajar, baju, roger, pelangi
for (b in 1:nrow(comparison_radius_boom)){
  
  if(grepl('Pelangi', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Pelangi'
    comparison_radius_boom$ageclass[b] = "Adult"
    
  }  else if (grepl('Fajar', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Fajar'
    #comparison_radius_boom$ageclass[b] = "Juvenile"
    comparison_radius_boom$ageclass[b] = "non-Adult"
    
  } else if (grepl('Baju', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Baju'
    #comparison_radius_boom$ageclass[b] = "Subadult"
    comparison_radius_boom$ageclass[b] = "non-Adult"
    
  } else if (grepl('Roger', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Roger'
    comparison_radius_boom$ageclass[b] = "Adult"
    
  } else {
    
    comparison_radius_boom$ID[b] = 'NA'
    comparison_radius_boom$ageclass[b] = "NA"
  }
  
}


# 02b: correlation matrix ----

comparison_radius_boom_numeric_adult <- comparison_radius_boom %>% 
  ungroup() %>%
  dplyr::filter(ageclass == "Adult") %>% 
  select(radius, ampl, pitch, entropy, specCentroid)
  #select(-audiofile, -match,
  #       -videofile, -voiced,
  #       -fs_video, -frame) 


cor_all <- cor(comparison_radius_boom_numeric)
cor_all_2 <- rcorr(as.matrix(comparison_radius_boom_numeric_adult))


#correlation plot, show all
#corrplot(cor_all_2$r[1,1:55, drop=FALSE], type = "upper", order = "original", 
#         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_2$r[1,1:5, drop=FALSE], 
         p.mat = cor_all_2$P[1,1:5, drop=FALSE],
         sig.level = 0.05, insig = "blank", diag = FALSE,
         tl.col = "black", #tl.srt = 90,
         addCoef.col = 'black',
         cl.pos = 'r')#, col = COL2('BrBG'))


# 03: visualization -----

# 03a: adding meta data for plotting: sex + ID -----



# 03b: scatter plots acoustics vs. radius ----
 
# scatter amplitude - sex

scatter_ampl_sex <- comparison_radius_boom %>% 
  #ggplot(aes(y = ampl_mean, x = radius, fill = sex, color = sex))+
  ggplot(aes(y = ampl, x = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Amplitude')+
  xlab('Airsac Radius [px]')+
  #coord_cartesian(xlim = c(0,10))+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

# amplitude mean - ID

scatter_ampl_ID <- comparison_radius_boom %>%
  ggplot(aes(y = ampl, x = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Amplitude')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

# amplitude - age

scatter_ampl_age <- comparison_radius_boom %>%
  ggplot(aes(y = ampl, x = radius, fill = ageclass, color = ageclass))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Amplitude')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none',
        axis.title.x = element_blank())

## pitch (f0) -sex

scatter_pitch_sex <- comparison_radius_boom %>%
  ggplot(aes(y = pitch, x = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Pitch')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_blank())

## pitch (f0) -ID

scatter_pitch_ID <- comparison_radius_boom %>%
  ggplot(aes(y = pitch, x = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Pitch')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(axis.title.x = element_blank())


## pitch (f0) - age

scatter_pitch_age <- comparison_radius_boom %>%
  ggplot(aes(y = pitch, x = radius, fill = ageclass, color = ageclass))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Pitch')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))#+
  theme(axis.title.x = element_blank())

## entropy  -sex

scatter_entropy_sex <- comparison_radius_boom %>%
  ggplot(aes(y = entropy, x = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Entropy')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## entropy - ID
scatter_entropy_ID <- comparison_radius_boom %>%
  #ggplot(aes(y = entropy_mean, x = radius, fill = ID, color = ID))+
  ggplot(aes(y = entropy, x = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Entropy')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## entropy - age
scatter_entropy_age <- comparison_radius_boom %>%
  ggplot(aes(y = entropy, x = radius, fill = ageclass, color = ageclass))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Entropy')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))+
  theme(legend.position = 'none')

## spectral Centroid - sex

scatter_spec_Centroid_sex <- comparison_radius_boom %>%
  ggplot(aes(y = specCentroid, x = radius, fill = sex, color = sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Centroid')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

## spectral Centroid - ID

scatter_spec_Centroid_ID <- comparison_radius_boom %>%
  ggplot(aes(y = specCentroid, x = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Centroid')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

## spectral Centroid  - age

scatter_spec_Centroid_age <- comparison_radius_boom %>%
  ggplot(aes(y = specCentroid, x = radius, fill = ageclass, color = ageclass))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Centroid')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

## spectral slope
scatter_spec_Slope_ID <- comparison_radius_boom %>%
  ggplot(aes(y = specSlope, x = radius, fill = ID, color = ID))+
  geom_point()+
  geom_smooth(method = 'lm')+
  ylab('Spectral Slope')+
  xlab('Airsac Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))


## plot grids ----

# plot grid sex
cowplot::plot_grid(scatter_ampl_sex, scatter_pitch_sex,
                   scatter_entropy_sex, scatter_spec_Centroid_sex,
                   rel_widths = c(0.45,0.55), ncol = 2)
#labels = c("A", "B", "C", "D"), ncol = 2)

#plot grid ID
cowplot::plot_grid(scatter_ampl_ID, scatter_pitch_ID,
                   scatter_entropy_ID, scatter_spec_Centroid_ID,
                   rel_widths = c(0.45,0.55), ncol = 2)

# plot grid ageclass

cowplot::plot_grid(scatter_ampl_age, scatter_pitch_age,
                   scatter_entropy_age, scatter_spec_Centroid_age,
                   rel_widths = c(0.45,0.55), ncol = 2)

# 03c: radius inflation over time ----

# timeseries plots Fajar
plot_fajar <- comparison_radius_boom %>%
  dplyr::filter(ID == 'Fajar') %>% 
  ggplot(aes(x = frame, y = radius, group = videofile, color = videofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Fajar')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Pelangi
plot_pelangi <- comparison_radius_boom %>%
  dplyr::filter(ID == 'Pelangi') %>% 
  ggplot(aes(x = frame, y = radius, group = videofile, color = videofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Pelangi')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Baju
plot_baju <- comparison_radius_boom %>%
  dplyr::filter(ID == 'Baju') %>% 
  ggplot(aes(x = frame, y = radius, group = videofile, color = videofile))+
  geom_line(size = 1.4)+
  xlab('Frame')+
  ylab('Airsac Radius [px]')+
  ggtitle('Baju')+
  theme_minimal()#+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

#timeseries plots Roger
plot_roger <- comparison_radius_boom %>%
  dplyr::filter(ID == 'Roger') %>% 
  ggplot(aes(x = frame, y = radius, group = videofile, color = videofile))+
  geom_line(size = 1.4)+
  xlab('Time [ms]')+
  ylab('Airsac Radius [px]')+
  ggtitle('Roger')+
  theme_minimal()+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())

cowplot::plot_grid(plot_fajar, plot_pelangi, plot_baju, plot_roger,
                   ncol = 2)

# deflation-inflation example Baju

baju_example <- "_Opp_June_13_Session_1_zoom_syncedboom_7_3_BajuDLC_resnet101_Deep_AirSacTrackingV1Jan1shuffle1_500000.csv"

detail_plot_baju <- comparison_radius_boom %>%
  dplyr::filter(videofile == baju_example) %>% 
  ggplot()+
  geom_line(aes(x = frame, y = radius), color = "black",size = 1.4)+
  geom_line(aes(x = frame, y = ampl), color = "red", size = 1.4)+
  geom_line(aes(x = frame, y = pitch), color = "green", size = 1.4)+
  geom_line(aes(x = frame, y = entropy), color = "yellow", size = 1.4)+
  xlab('Frame')+
  ylab('Parameters')+
  ggtitle('Baju')+
  theme_minimal()#+
  theme(legend.position = 'none',
        text = element_text(size = 20))#,
#axis.title.y = element_blank())