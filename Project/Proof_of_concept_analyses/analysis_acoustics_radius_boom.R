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

install_load("tidyverse","conicfit","ggcorrplot", "scales", "spiro",
             "signal", "foreach", "Hmisc", "ggplotify", "cowplot")

# 01: data ----

## 01a: load data ----
acoustics_boom <- readRDS("acoustic_analysis_boom_multi_checked.rds")

radius_boom <- readRDS("radius_estimation_boom_proof1.rds")

# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# we only use colors from a colorblind-safe palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

# colors for sex  (male/female) :  "#44AA99", "#CC6677"
# colors for age classes (juvenile/ sub-adult) : "#882255", "#117733"
# color adults pooled: black
# color non-adults pooled: black

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
    comparison_radius_boom$ageclass_detailed[b] = "Adult"
    comparison_radius_boom$ageclass[b] = "Adult"
    
  }  else if (grepl('Fajar', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Fajar'
    comparison_radius_boom$ageclass_detailed[b] = "Juvenile"
    comparison_radius_boom$ageclass[b] = "non-Adult"
    
  } else if (grepl('Baju', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Baju'
    comparison_radius_boom$ageclass_detailed[b] = "Subadult"
    comparison_radius_boom$ageclass[b] = "non-Adult"
    
  } else if (grepl('Roger', comparison_radius_boom$videofile[b], ignore.case = TRUE) == TRUE){
    
    comparison_radius_boom$ID[b] = 'Roger'
    comparison_radius_boom$ageclass_detailed[b] = "Adult"
    comparison_radius_boom$ageclass[b] = "Adult"
    
  } else {
    
    comparison_radius_boom$ID[b] = 'NA'
    comparison_radius_boom$ageclass[b] = "NA"
  }
  
}

# 02b: stats ----

# number of call sequences

adults <- comparison_radius_boom %>% 
  dplyr::filter(ageclass == "Adult")

n_calls <- unique(adults$videofile)
n_calls <- length(n_calls)

# number of call sequences: 25

# number of datapoints

n_datapoints <- nrow(adults)

# number of datapoints: 176 

# 02c: correlation matrix ----

## 02b-1: adults ----

comparison_radius_boom_numeric_adult <- comparison_radius_boom %>% 
  ungroup() %>%
  dplyr::filter(ageclass == "Adult") %>% 
  select(radius, ampl, pitch, entropy, specCentroid, f1_freq, f2_freq,
         harmEnergy, peakFreq, fmPurity, HNR)
  
#corrplots with ggcorrplot: http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2
corr <- as.data.frame(round(cor(comparison_radius_boom_numeric_adult, use = "pairwise.complete.obs", method = "pearson"),1))

corr_sub <- corr %>% 
  dplyr::filter(radius != 1) %>% # to filter out radius v radius, as we can't do that in the ggcorrplot function
  select("radius") 

row.names(corr_sub) <- c("Amplitude", "F0", "Entropy", "Spectral Centroid", 
                         "F1", "F2", "harmonic Energy", "peak Frequency", "fm Purity",
                         "HNR")
colnames(corr_sub) <- c("Radius")

p.mat <- cor_pmat(corr)

corrplot_adults <- ggcorrplot(corr_sub, method = "circle", lab = TRUE, p.mat = p.mat[1,2:11], insig = "blank") 

ggsave("corrplot_adults.jpg", dpi = 300, width = 10, height = 3) 


# old version with corrplot

#cor_all_2 <- rcorr(as.matrix(comparison_radius_boom_numeric_adult))

# corrplot(cor_all_2$r[1,1:11, drop=FALSE], 
#         p.mat = cor_all_2$P[1,1:11, drop=FALSE],
#         sig.level = 0.01, insig = "blank", diag = FALSE,
#         tl.col = "black", #tl.srt = 90,
#         addCoef.col = 'black',
#         cl.pos = 'r') #, col = COL2('BrBG'))
 


 
 ## 02b-2: non-adults ----
 
 comparison_radius_boom_numeric_nonadult <- comparison_radius_boom %>% 
   ungroup() %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   select(radius, ampl, pitch, entropy, specCentroid, f1_freq, f2_freq,
          harmEnergy, peakFreq, fmPurity, HNR)


#corrplots with ggcorrplot: http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2
corr_na <- as.data.frame(round(cor(comparison_radius_boom_numeric_nonadult, use = "pairwise.complete.obs", method = "pearson"),1))

corr_sub_na <- corr_na %>% 
  dplyr::filter(radius != 1) %>% # to filter out radius v radius, as we can't do that in the ggcorrplot function
  select("radius") 

row.names(corr_sub_na) <- c("Amplitude", "F0", "Entropy", "Spectral Centroid", 
                         "F1", "F2", "harmonic Energy", "peak Frequency", "fm Purity",
                         "HNR")
colnames(corr_sub_na) <- c("Radius")

p.mat_na <- cor_pmat(corr_na)

corrplot_nonadults <- ggcorrplot(corr_sub_na, method = "circle", lab = TRUE, p.mat = p.mat_na[1,2:11], insig = "blank") 

ggsave("corrplot_nonadults.svg", dpi = 300, width = 10, height = 3) 
 
# old version with corrplot

 cor_all_2_nonadult <- rcorr(as.matrix(comparison_radius_boom_numeric_nonadult))
 
 #correlation plot, only show significant correlations for non-adults
 
 corrplot(cor_all_2_nonadult$r[1,1:11, drop=FALSE], 
          p.mat = cor_all_2_nonadult$P[1,1:11, drop=FALSE],
          sig.level = 0.05, insig = "blank", diag = FALSE,
          tl.col = "black", #tl.srt = 90,
          addCoef.col = 'black',
          cl.pos = 'r' , col = COL2("BuRd"))
 
# correlation_plot <-recordPlot()
 
# 03: visualization -----

## 03a: adults ----
## acoustics vs. radius for adults only (Pelangi & Roger)
## plotted acoustic parameters: amplitude, pitch (f0), entropy and Spectral Centroid 
 
 
### adults - pooled - amplitude
 
 scatter_adults_ampl_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = ampl, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Amplitude')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=120, y=0.4, label= " R² = 0.45")
 
### adults - pooled - pitch (f0) 
 
 scatter_adults_pitch_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = pitch, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Pitch [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(200,450),
                      breaks= c(200,250,300,350,400))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=120, y=400, label= " R² = 0.82")

### adults - pooled - entropy
 
 scatter_adults_entropy_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = entropy, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Entropy')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=120, y=0.4, label= " R² = -0.31")

### adults - pooled - spectral Centroid
 
 scatter_adults_specCentroid_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = specCentroid, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Spectral Centroid [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(1000,6000),
                      breaks= c(1000,2000,3000,4000, 5000))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x = 120, y = 5000, label = " R² = -0.55")

 ### adults divided by sex
 ### adults - by sex - amplitude
 
 scatter_adults_ampl_sex <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = ampl, x = radius, fill = sex, color = sex))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Amplitude')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Sex",
                     labels = c("Female", "Male"),
                     values = c("#CC6677", "#44AA99"))+
   scale_color_manual(name = "Sex",
                      labels = c("Female", "Male"),
                      values = c("#CC6677", "#44AA99"))#+
   #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### adults - by sex - pitch (f0)
 
 scatter_adults_pitch_sex <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = pitch, x = radius, fill = sex, color = sex))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Pitch [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(200,450),
                      breaks= c(200,250,300,350,400))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Sex",
                     labels = c("Female", "Male"),
                     values = c("#CC6677", "#44AA99"))+
   scale_color_manual(name = "Sex",
                      labels = c("Female", "Male"),
                      values = c("#CC6677", "#44AA99"))#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### adults - by sex - entropy 
 
 scatter_adults_entropy_sex <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = entropy, x = radius, fill = sex, color = sex))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Entropy')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Sex",
                     labels = c("Female", "Male"),
                     values = c("#CC6677", "#44AA99"))+
   scale_color_manual(name = "Sex",
                      labels = c("Female", "Male"),
                      values = c("#CC6677", "#44AA99"))#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### adults - by sex - spectral Centroid
 
 scatter_adults_specCentroid_sex <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "Adult") %>% 
   ggplot(aes(y = specCentroid, x = radius, fill = sex, color = sex))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Spectral Centroid [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(1000,6000),
                      breaks= c(1000,2000,3000,4000, 5000))+
   theme(text = element_text(size = 20))+
   scale_fill_manual(name = "Sex",
                     labels = c("Female", "Male"),
                     values = c("#CC6677", "#44AA99"))+
   scale_color_manual(name = "Sex",
                      labels = c("Female", "Male"),
                      values = c("#CC6677", "#44AA99"))
   #scale_colour_manual(name = "Sex",
   #                    labels = c("Male", "Female"))#,
   #                     values = c("#44AA99","#CC6677"))#+
   #theme(legend.position = 'none')#+
  #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 legend_adults_sex <- get_legend(scatter_adults_specCentroid_sex)
 
 scatter_adults_specCentroid_sex <- scatter_adults_specCentroid_sex + theme(legend.position = "none")
 
## 03b: non-adults ----

 ### non- adults - pooled - amplitude
 
 scatter_nonadults_ampl_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = ampl, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', , color = "grey15")+
   ylab('Amplitude')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=80, y=0.4, label= " R² = -0.03")
 
 ### non-adults - pooled - pitch (f0) 
 
 scatter_nonadults_pitch_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = pitch, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', , color = "grey15")+
   ylab('Pitch [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(200,450),
                      breaks= c(200,250,300,350,400))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=80, y=400, label= " R² = 0.01")
 
 ### non-adults - pooled - entropy
 
 scatter_nonadults_entropy_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = entropy, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Entropy')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x=80, y=0.4, label= " R² = 0.03")
 
 ### non-adults - pooled - spectral Centroid
 
 scatter_nonadults_specCentroid_pooled <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = specCentroid, x = radius))+
   geom_point()+
   geom_smooth(method = 'lm', color = "grey15")+
   ylab('Spectral Centroid [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(1000,6000),
                      breaks= c(1000,2000,3000,4000, 5000))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   annotate("text", x = 80, y = 5000, label = " R² = -0.01")
 
 ### non-Adutls divided by detailed age class (juvenile/sub-adult)
 ### nonadults - by ageclass - amplitude
 
 scatter_nonadults_ampl_ageclass <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = ampl, x = radius, fill = ageclass_detailed, color = ageclass_detailed))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Amplitude')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Ageclass",
                     labels = c("Subadult", "Juvenile"),
                     values = c("#117733", "#882255"))+
   scale_color_manual(name = "Ageclass",
                      labels = c("Subadult", "Juvenile"),
                      values = c("#117733", "#882255"))#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### nonadults - by ageclass - pitch (f0)
 
 scatter_nonadults_pitch_ageclass <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = pitch, x = radius, fill = ageclass_detailed, color = ageclass_detailed))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Pitch [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(200,450),
                      breaks= c(200,250,300,350,400))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Ageclass",
                     labels = c("Subadult", "Juvenile"),
                     values = c("#117733", "#882255"))+
   scale_color_manual(name = "Ageclass",
                      labels = c("Subadult", "Juvenile"),
                      values = c("#117733", "#882255"))#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### nonadults - by ageclass - entropy 
 
 scatter_nonadults_entropy_ageclass <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = entropy, x = radius, fill = ageclass_detailed, color = ageclass_detailed))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Entropy')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(0,0.5),
                      breaks= c(0,0.1,0.2,0.3,0.4))+
   theme(text = element_text(size = 20))+
   theme(legend.position = 'none')+
   scale_fill_manual(name = "Ageclass",
                     labels = c("Subadult", "Juvenile"),
                     values = c("#117733", "#882255"))+
   scale_color_manual(name = "Ageclass",
                      labels = c("Subadult", "Juvenile"),
                      values = c("#117733", "#882255"))#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 
 ### nonadults - by ageclass - spectral Centroid
 
 scatter_nonadults_specCentroid_ageclass <- comparison_radius_boom %>%
   dplyr::filter(ageclass == "non-Adult") %>% 
   ggplot(aes(y = specCentroid, x = radius, fill = ageclass_detailed, color = ageclass_detailed))+
   geom_point()+
   geom_smooth(method = 'lm')+
   ylab('Spectral Centroid [Hz]')+
   xlab('Airsac Radius [px]')+
   #coord_cartesian(xlim = c(0,10))+
   theme_minimal()+
   scale_y_continuous(expand = c(0, 0),
                      limits= c(1000,6000),
                      breaks= c(1000,2000,3000,4000, 5000))+
   theme(text = element_text(size = 20))+
   scale_fill_manual(name = "Ageclass",
                     labels = c("Subadult", "Juvenile"),
                     values = c("#117733", "#882255"))+
   scale_color_manual(name = "Ageclass",
                      labels = c("Subadult", "Juvenile"),
                      values = c("#117733", "#882255"))
 #scale_colour_manual(name = "Sex",
 #                    labels = c("Male", "Female"))#,
 #                     values = c("#44AA99","#CC6677"))#+
 #theme(legend.position = 'none')#+
 #annotate("text", x=100, y=0.4, label= " R² = 0.45, p = ")
 legend_nonadults_ageclass <- get_legend(scatter_nonadults_specCentroid_ageclass)
 
 scatter_nonadults_specCentroid_ageclass <- scatter_nonadults_specCentroid_ageclass + theme(legend.position = "none")
 
 
## 03c: plot grids ----

# plot grid adults pooled


 cowplot::plot_grid(scatter_adults_ampl_pooled, scatter_adults_pitch_pooled,
                   scatter_adults_entropy_pooled, scatter_adults_specCentroid_pooled,
                   ncol = 4, align = "h", labels = c("A", "B", "C", "D"))

#plot grid adults by sex
 
cowplot::plot_grid(scatter_adults_ampl_sex, scatter_adults_pitch_sex,
                   scatter_adults_entropy_sex, scatter_adults_specCentroid_sex,legend_adults_sex,
                   ncol = 5, labels = c("F", "G", "H", "I", ""))

# plot grid ageclass

cowplot::plot_grid(scatter_nonadults_ampl_pooled, scatter_nonadults_pitch_pooled,
                   scatter_nonadults_entropy_pooled, scatter_nonadults_specCentroid_pooled,
                   ncol = 4, align = "h", labels = c("J", "K", "L", "M"))

# plot grid ageclass by ageclass detailed

cowplot::plot_grid(scatter_nonadults_ampl_ageclass, scatter_nonadults_pitch_ageclass,
                   scatter_nonadults_entropy_ageclass, scatter_nonadults_specCentroid_ageclass,legend_nonadults_ageclass,
                   ncol = 5, align = "h", labels = c("N", "O", "P", "Q", ""))

  # 04: old plots -----
  
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
    theme(legend.position = 'none')
  
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
          axis.title.x = element_blank())+
    scale_colour_manual(name = "Age Class",
                        labels = c("Adult", "Juvenile"),
                        values = c("darkred", "grey60"))+
    scale_fill_manual(name = "Age Class",
                      labels = c("Adult", "Juvenile"),
                      values = c("darkred", "grey40"))
  
  ## pitch (f0) -sex
  
  scatter_pitch_sex <- comparison_radius_boom %>%
    ggplot(aes(y = pitch, x = radius, fill = sex, color = sex))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Pitch')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20))+
    theme(legend.position = 'none')
  
  ## pitch (f0) -ID
  
  scatter_pitch_ID <- comparison_radius_boom %>%
    ggplot(aes(y = pitch, x = radius, fill = ID, color = ID))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Pitch')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20))+
    theme(axis.title.x = element_blank(),
          legend.position = 'none')
  
  
  ## pitch (f0) - age
  
  scatter_pitch_age <- comparison_radius_boom %>%
    ggplot(aes(y = pitch, x = radius, fill = ageclass, color = ageclass))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Pitch')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20))+
    theme(axis.title.x = element_blank(),
          legend.position = 'none')+
    scale_colour_manual(name = "Age Class",
                        labels = c("Adult", "Juvenile"),
                        values = c("darkred", "grey60"))+
    scale_fill_manual(name = "Age Class",
                      labels = c("Adult", "Juvenile"),
                      values = c("darkred", "grey40"))
  
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
    theme(legend.position = 'none',
          axis.title.x = element_blank())
  
  ## entropy - age
  scatter_entropy_age <- comparison_radius_boom %>%
    ggplot(aes(y = entropy, x = radius, fill = ageclass, color = ageclass))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Entropy')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20))+
    theme(legend.position = 'none',
          axis.title.x = element_blank())+
    scale_colour_manual(name = "Age Class",
                        labels = c("Adult", "Juvenile"),
                        values = c("darkred", "grey60"))+
    scale_fill_manual(name = "Age Class",
                      labels = c("Adult", "Juvenile"),
                      values = c("darkred", "grey40"))
  
  ## spectral Centroid - sex
  
  scatter_spec_Centroid_sex <- comparison_radius_boom %>%
    ggplot(aes(y = specCentroid, x = radius, fill = sex, color = sex))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Spectral Centroid')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20),
          legend.position = 'none')
  
  ## spectral Centroid - ID
  
  scatter_spec_Centroid_ID <- comparison_radius_boom %>%
    ggplot(aes(y = specCentroid, x = radius, fill = ID, color = ID))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Spectral Centroid')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20),
          axis.title.x = element_blank(),
          legend.position = 'none')
  
  ## spectral Centroid  - age
  
  scatter_spec_Centroid_age <- comparison_radius_boom %>%
    ggplot(aes(y = specCentroid, x = radius, fill = ageclass, color = ageclass))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Spectral Centroid')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20),
          axis.title.x = element_blank(),
          legend.position = 'none')+
    scale_colour_manual(name = "Age Class",
                        labels = c("Adult", "Juvenile"),
                        values = c("darkred", "grey60"))+
    scale_fill_manual(name = "Age Class",
                      labels = c("Adult", "Juvenile"),
                      values = c("darkred", "grey40"))
  
  ## spectral slope
  scatter_spec_Slope_ID <- comparison_radius_boom %>%
    ggplot(aes(y = specSlope, x = radius, fill = ID, color = ID))+
    geom_point()+
    geom_smooth(method = 'lm')+
    ylab('Spectral Slope')+
    xlab('Airsac Radius [px]')+
    theme_minimal()+
    theme(text = element_text(size = 20))
  
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
  
  
  ##03e: Spectral Centroid over time vs radius over time
  
  plot_roger <- comparison_radius_boom %>%
    dplyr::filter(ID == 'Roger') %>% 
    ggplot()+
    geom_line(aes(x = frame, y = radius, group = videofile), size = 1.4, color = "black")+
    geom_line(aes(x = frame, y = specCentroid, group = videofile), size = 1.4, color = "green")+
    xlab('Time [ms]')+
    ylab('Airsac Radius [px]')+
    ggtitle('Roger')+
    theme_minimal()+
    scale_y_log10()+
    theme(legend.position = 'none',
          text = element_text(size = 20))#,
  