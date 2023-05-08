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

install_load("tidyverse", "effsize", "psych", "signal", "foreach", "kza")

# 01a: functions ----

# function explanation goes here
butter.it <- function(x, samplingrate, order, lowpasscutoff)
{bf <- butter(order,lowpasscutoff/samplingrate, type="low") #normalized frequency
x <<- as.numeric(signal::filtfilt(bf, x))} 

# 01b: load data ----

# we need: manual trackings combined (we don't need to have the code of how we combined them,just load the data)

manual_radius <- read_delim("manually_tracked_airsac_radii.csv", delim = ',')

# we need the corresponding dlc radii

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

# saved every of the 5 best examples in a seperate dataframe df_exp1, df_exp2 and so on

df_full<- rbind(df_exp1, df_exp2, df_exp3, df_exp4, df_exp5)

saveRDS(df_full, file = "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.rds")
write.table(df_full, "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.csv", sep = ",")

# combine correlation grouped by parameter combination and grouped by video for statistics to report in manuscript

correlation_per_video <- df_full %>% 
  group_by(videoname) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE)) 

# per example 
correlation_per_setting <- df_full %>% 
  group_by(example_nr) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE)) 
# reported in manuscript are: mean/sd/min/max of the max_cor and max_cor_sm, to indicate variability between videos 
# max_cor and max_cor_sm are indicative of best parameter combination per video 

best_combos_all <- df %>% 
  group_by(alpha, beta, thresh_div1, thresh_div2, dilation, medianblur) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed), median_cor_sm = median(cor_radius_smoothed),
            min_cor_sm = min(cor_radius_smoothed), max_cor_sm = max(cor_radius_smoothed),
            sd_cor_sm = sd(cor_radius_smoothed, na.rm = TRUE), sd_cor = sd(cor_radius, na.rm = TRUE),
            mean_cor = mean(cor_radius, na.rm = TRUE), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE),
            mean_cor_x = mean(cor_x, na.rm = TRUE), mean_cor_y = mean(cor_y, na.rm = TRUE)) 

which.max(best_combos_all$mean_cor_sm) #best parameter combo, reported in manuscript are mean, sd, min, max for both smoothed and raw radii

# looking at correlations in detail for the best setting combination for all videos

# search for string: c1_8_c2_10_al_2_b_35_dil_7_blur_35 (best parameter combination overall)
list_of_files_df <- as.data.frame(list_of_files)

list_of_files_subset <- list_of_files_df %>% 
  dplyr::filter(grepl('c1_8_c2_10_al_2_b_35_dil_7_blur_35', list_of_files, ignore.case = TRUE) == TRUE)

temp_df <- data.frame()

#reading each file within the range and append them to create one file

for (i in 1:nrow(list_of_files_subset)){
  #read the file
  current_File <- read_delim(paste(path, list_of_files_subset[i,1], sep = "\\"), delim = ",")
  
  videoname <- list_of_files_subset[i,1]
  videoname <- substring(videoname, 1, nchar(videoname)-4)
  videoname <- str_split(videoname, pattern = "_")
  videoname <- videoname[[1]]
  videoname <- paste(videoname[14], videoname[15], sep = "_")
  
  current_File <- current_File %>% 
  dplyr::mutate(smoothed_hough_radius = butter.it(current_File$r,
                                             samplingrate =  25,
                                             order = 2,
                                             lowpasscutoff = 10),
                videoname = videoname)
  #Append the current file
  temp_df = rbind(temp_df, current_File)    
}

joined_radii_best_hough <- plyr::join(temp_df, manual_radius, by= c("name"), type="left", match="first")

joined_radii_best_hough <- joined_radii_best_hough %>% 
  dplyr::filter(radius_man > 100)
# Visualization

joined_radii_best_hough %>% 
  #dplyr::filter(videoname == "June16_07" | videoname == "June09_03") %>% 
  ggplot(aes(x= radius_man, y = smoothed_hough_radius))+
  geom_point(aes( fill = videoname, color = videoname), size = 2)+
  geom_smooth(method = 'lm')+
  ylab('Automatically tracked Radius [px], Hough')+
  xlab('Manually labeled Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))

# 03: dlc combining datasets + preparations ----

# 03a: preparations
dlc_radius <- dlc_radius[!is.na(dlc_radius$radius),]
##
dlc_test <- dlc_radius %>%
  group_by(videoname) %>%
  group_modify(
      mutate(kol_smooth = kza(radius, k = 4,m = 3)$kza))
      
  group_map(~ kza(hough_radius$r, k = 4,m = 3)$kza)

dlc_radius <- dlc_radius %>% 
  mutate(smoothed_dlc_radius = butter.it(dlc_radius$radius,
                                         samplingrate =  25,
                                         order = 2,
                                         lowpasscutoff = 10))

# we need: framenr_52_framevid_June09_03, we have frame in column frame and videoname in filename as first part before DLC
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

joined_radii_dlc <- plyr::join(dlc_radius, manual_radius, by= c("name"), type="left", match="first")

# removing all lines, where no circle was tracked manually (indicated by very small circles)
joined_radii_dlc <- joined_radii_dlc %>% 
  dplyr::filter(radius_man > 100 & radius < 270) # we apply the filter to only take into accoutn automatically tracked radii below 250 px, to make it comparable to Hough circles, that can maximally find circles of 250 px
# double check with Wim and check if 250 was the limit in hough transform 

# 03b: correlation to test if tracking success is sufficient

dlc_correlations <- joined_radii_dlc %>% 
  group_by(videoname) %>% 
  summarise(cor_radius = cor(radius, radius_man),
            cor_radius_sm = cor(smoothed_dlc_radius, radius_man))

# plotting radius comprison

dlc <- joined_radii_dlc %>% 
  ggplot(aes(x= radius_man, y = radius))+
  geom_point(aes( fill = videoname, color = videoname), size = 2, alpha = 0.6)+
  geom_smooth(method = 'lm')+
  ylab('Automatically tracked Radius [px], DLC')+
  xlab('Manually labeled Radius [px]')+
  theme_minimal()+
  theme(text = element_text(size = 20))#,
        legend.position = "none")

cor_r <- corr.test(joined_radii_dlc$radius, joined_radii_dlc$radius_man)$r
cor_r_sm <- corr.test(joined_radii_dlc$smoothed_dlc_radius, joined_radii_dlc$radius_man)

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
