# comparison of automatically tracked radii from different approaches
# to manually tracked radii, via correlation it can be decided if tracking success
# is sufficient and which approach works best

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: radii values from Hough Transform or DLC trackings and subsequent 
# circle estimation & datasubset of manually tracked radii
# end point: correlation between manual and automaticlly tracked radii 

############

# To do:
# 1) check that column names match, are understandable and changed accordingly everywhere in code

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse", "effsize", "psych", "signal", "foreach", "kza", "plyr")

# 01a: functions ----

# function explanation goes here
butter.it <- function(x, samplingrate, order, lowpasscutoff)
{bf <- butter(order,lowpasscutoff/samplingrate, type="low") #normalized frequency
x <<- as.numeric(signal::filtfilt(bf, x))} 

# 01b: load data ----

# we need: manual trackings combined (we don't need to have the code of how we combined them,just load the data)

manual_radius <- read_delim("manually_tracked_airsac_radii.csv", delim = ',')

# we need the corresponding dlc radii

# where does this come from? result of from dlc to circle?
dlc_radius <- readRDS("dlc_for_comparison_to_manualtracks.rds")

# hough data is not just a single datafile, which is why it is loaded in a loop later

# 02: Hough Transform ----

# example of what to read in: folder D:\PostDoc\Donders\AirSacTracker_2\Tracker_Backup_8May2023\module_hough\results\example 1
# in this folder we have the tracked videos and csv files of one parameter combination
df <- data.frame()
joined_radii_all <- data.frame()

path <- choose.dir()
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

for (a in 1: length(list_of_files)){
  
  hough_radius <- read_delim(paste(path, list_of_files[a], sep = "\\"), delim = "," )
  
  hough_radius <- hough_radius %>% 
    select(-...1)
  
  #retrieve videoname from parameter optimization csv
  videoname <- list_of_files[a]
  videoname <- substring(videoname, 1, nchar(videoname)-4)
  #videoname <- str_split(videoname, pattern = "_")
  #videoname <- videoname[[1]]
  #videoname <- paste(videoname[4], videoname[5], sep = "_")
  
  hough_radius <- hough_radius %>%
    drop_na() 
  hough_radius <- hough_radius %>% 
    mutate(smoothed_hough_radius_butter = butter.it(hough_radius$r,
                                             samplingrate =  25,
                                             order = 2,
                                             lowpasscutoff = 10),
           smoothed_hough_radius_kolmogorov = kza(hough_radius$r, k = 4,m = 3)$kza )
  
  joined_radii <- plyr::join(hough_radius, manual_radius, by= c("name"), type="left", match="first")
  joined_radii <- joined_radii %>% 
    mutate(videoname = videoname)
  
  cor_r <- corr.test(joined_radii$r, joined_radii$radius_man)
  cor_r_sm_butter <- corr.test(joined_radii$smoothed_hough_radius_butter, joined_radii$radius_man)
  cor_r_sm_kol <- corr.test(joined_radii$smoothed_hough_radius_kolmogorov, joined_radii$radius_man)
  cor_x <- corr.test(joined_radii$x, joined_radii$X)
  cor_y <- corr.test(joined_radii$y, joined_radii$Y)

  # saving relevant parameters 
  
  df[a,1] <- videoname
  df[a,2] <- cor_r$r
  df[a,3] <- cor_r_sm_kol$r
  df[a,4] <- cor_x$r
  df[a,5] <- cor_y$r
  df[a,6] <- cor_r_sm_butter$r
  df[a,7] <- joined_radii$examplenr[1]
  #df[a,6] <- hough_radius$alpha[1]
  #df[a,7] <- hough_radius$beta[1]
  #df[a,8] <- hough_radius$threh_div1[1]
  #df[a,9] <- hough_radius$threh_div2[1]
  #df[a,10] <- hough_radius$dilation[1]
  #df[a,11] <- hough_radius$medianblur[1]
  #df[a,12] <- nrow(joined_radii)
  #df[a,13] <- nrow(hough_radius)

  joined_radii_all <- rbind(joined_radii_all, joined_radii)
  
}

colnames(df) <- c("videoname", "cor_radius","cor_radius_smoothed_kol", "cor_x",  "cor_y","cor_radius_smoothed_butter", "example_nr")#, 
  #                "alpha","beta","thresh_div1","thresh_div2","dilation",
  #                "medianblur", "nr_frames_used", "nr_frames_in_video")

# saved every of the 5 best examples in a seperate dataframe df_exp1, df_exp2 and so on by renaming the resulting df
# of former step by df_exp1 etc. corresponidng to data

df_full<- rbind(df_exp1, df_exp2, df_exp3, df_exp4, df_exp5)

saveRDS(df_full, file = "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.rds")
write.table(df_full, "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.csv", sep = ",")

# combine correlation grouped by parameter combination and grouped by video for statistics to report in manuscript

correlation_per_video <- df_full %>% 
  dplyr::group_by(videoname) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE)) 

# per example 
correlation_per_setting <- df_full %>% 
  dplyr::group_by(example_nr) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE), .groups = 'drop') 


# 03: dlc combining datasets + preparations ----

# 03a: preparations
dlc_radius <- dlc_radius[!is.na(dlc_radius$radius),]

# to match with manual data we need the correct name format 
# we need something like: framenr_52_framevid_June09_03, we have frame in column frame and videoname in filename as first part before DLC
# -> split filename at DLC, use first in paste with frame and the respective words
match_name <- data.frame()
for(i in 1:length(dlc_radius$videofile)){
  names <- str_split(dlc_radius$videofile[i], pattern = "DLC")
  names <- names[[1]]
  match_name[i,1] <-names[[1]]
  match_name[i,2] <- dlc_radius$frame[i]
}

dlc_radius <- dlc_radius %>% 
  dplyr::mutate(name = paste("framenr", match_name$V2, "framevid", match_name$V1, sep = "_"),
                radius = as.numeric(radius),
                videoname = match_name$V1) 

# smoothing per video: kolmogorov-zurbenko 
dlc_radius <- dlc_radius %>%
  group_by(videoname) %>% 
  mutate(kza_smoothed_radii = kza(radius, k = 4,m = 3)$kza)

joined_radii_dlc <- plyr::join(dlc_radius, manual_radius, by= c("name"), type="left", match="first")

# removing all lines, where no circle was tracked manually (indicated by very small circles)
joined_radii_dlc <- joined_radii_dlc %>% 
  dplyr::filter(radius_man > 100 & radius < 270 & kza_smoothed_radii < 270) # we apply the filter to only take into accoutn automatically tracked radii below 270 px, to make it comparable to Hough circles, that can maximally find circles of 270 px

# 03b: correlation to test if tracking success is sufficient

# dlc correl
dlc_correlations <- joined_radii_dlc %>% 
  group_by(videoname) %>% 
  summarise(cor_radius = cor(radius, radius_man),
            cor_radius_sm = cor(kza_smoothed_radii, radius_man, use= "pairwise.complete.obs"))

dlc_correlations
#   cor_radius cor_radius_sm
# 1  0.8579891     0.8543827




# 04: visualization both approaches ----

# plotting radius comprison

# visualization DLC
dlc <- joined_radii_dlc %>% 
  ggplot(aes(x= radius_man, y = radius))+
  geom_point(aes( fill = videoname, color = videoname), size = 2, alpha = 0.6)+
  geom_smooth(method = 'lm')+
  ylab('Automatically tracked Radius [px], DLC')+
  xlab('Manually labeled Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

# visualization hough
hough <- joined_radii_all %>%
  dplyr::filter(radius_man >= 100) %>% 
  ggplot(aes(x= radius_man, y = smoothed_hough_radius_kolmogorov))+
  geom_point(aes( fill = videoname, color = videoname), size = 2, alpha = 0.6)+
  geom_smooth(method = 'lm')+
  ylab('Automatically tracked Radius [px], Hough')+
  xlab('Manually labeled Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

# labels legend ändern! can we put legend on the bottom of plot 2? 
cowplot::plot_grid(dlc, hough, rel_heights = c(0.5, 0.5), ncol = 1, labels = c("A", "B"))

ggsave("comparison_dlc_hough_manual_row.jpg", dpi = 300, width= 10, height = 14)
