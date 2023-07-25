# comparison of automatically tracked radii from different approaches
# to manually tracked radii, via correlation so it can be decided if tracking success
# is sufficient and which approach works best

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: radii values from Hough Transform or DLC trackings and subsequent 
# circle estimation & datasubset of manually tracked radii
# end point: correlation between manual and automaticlly tracked radii + visualization

############

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse", "effsize", "psych", "signal", "foreach", "kza", "plyr")

# 01a: functions ----


# 01b: load data ----

# we need: manual trackings combined (we don't need to have the code of how we combined them,just load the data)

# path in Github structure: AirSacTracker\Project\comparison_manual_auto
manual_radius <- read_delim("manually_tracked_airsac_radii.csv", delim = ',')

# we need the corresponding dlc radii

# path in Github structure: AirSacTracker\Project\comparison_manual_auto
dlc_radius <- readRDS("dlc_for_comparison_to_manualtracks.rds")

# hough tracking data is not just a single datafile, which is why it is loaded in a loop later

# 02: Hough Transform ----

## 02a: preparations Hough ----
df <- data.frame()
joined_radii_all <- data.frame()

# path used for this project (Github folder structure): AirSacTracker\Project\comparison_manual_auto\hough_trackings_for_comparison_to_manual

path <- getwd()
path <- paste0(path, '/hough_trackings_for_comparison_to_manual/')
pattern <- "*.csv"
list_of_files <- list.files(path = path, pattern = pattern)

for (a in 1: length(list_of_files)){
  
  hough_radius <- read_delim(paste(path, list_of_files[a], sep = "\\"), delim = "," )
  
  hough_radius <- hough_radius %>% 
    select(-...1)
  
  #retrieve videoname from parameter optimization csv
  videoname <- list_of_files[a]
  videoname <- substring(videoname, 1, nchar(videoname)-4)
  
  hough_radius <- hough_radius %>%
    drop_na() 
  hough_radius <- hough_radius %>% 
    mutate(smoothed_hough_radius_kolmogorov = kza(hough_radius$r, k = 4,m = 3)$kza)
  
  joined_radii <- plyr::join(hough_radius, manual_radius, by= c("name"), type="left", match="first")
  joined_radii <- joined_radii %>% 
    mutate(videoname = videoname)
  
  cor_r <- corr.test(joined_radii$r, joined_radii$radius_man)
  cor_r_sm_kol <- corr.test(joined_radii$smoothed_hough_radius_kolmogorov, joined_radii$radius_man)
  cor_x <- corr.test(joined_radii$x, joined_radii$X)
  cor_y <- corr.test(joined_radii$y, joined_radii$Y)

  # saving relevant parameters 
  
  df[a,1] <- videoname
  df[a,2] <- cor_r$r
  df[a,3] <- cor_r_sm_kol$r
  df[a,4] <- cor_x$r
  df[a,5] <- cor_y$r
  df[a,6] <- joined_radii$examplenr[1]

  joined_radii_all <- rbind(joined_radii_all, joined_radii)
  
}

colnames(df) <- c("videoname", "cor_radius","cor_radius_smoothed_kol", "cor_x",  "cor_y", "example_nr")#, 
  #                "alpha","beta","thresh_div1","thresh_div2","dilation",
  #                "medianblur", "nr_frames_used", "nr_frames_in_video")

saveRDS(df, file = "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.rds")
write.table(df, "hough_vs_manuallytracked_radii_correlation_all_parameter_combinations_5examples.csv", sep = ",")

## 02b: correlations Hough -----
# combine correlation grouped by parameter combination and grouped by video for statistics to report in manuscript

detach(package:plyr) # needs to be detached, otherwise grouping does not work in function below

correlation_per_video <- df %>%
  dplyr::group_by(videoname) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE)) 

# per example 
correlation_per_setting <- df %>% 
  group_by(example_nr) %>%                            
  summarise(mean_cor_sm = mean(cor_radius_smoothed_kol), median_cor_sm = median(cor_radius_smoothed_kol),
            min_cor_sm = min(cor_radius_smoothed_kol), max_cor_sm = max(cor_radius_smoothed_kol),
            mean_cor = mean(cor_radius, na.rm = TRUE ), median_cor = median(cor_radius, na.rm = TRUE),
            min_cor = min(cor_radius, na.rm = TRUE), max_cor = max(cor_radius, na.rm = TRUE)) 
# numbers from example 4 are reported in the manuscript, because of best mean correlation 

# 03: DLC ----

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

# dlc correlations
dlc_correlations <- joined_radii_dlc %>% 
  group_by(videoname) %>% 
  summarise(cor_radius = cor(radius, radius_man),
            cor_radius_sm = cor(kza_smoothed_radii, radius_man, use= "pairwise.complete.obs"))

dlc_correlations
#   cor_radius cor_radius_sm
# 1  0.8579891     0.8543827

# 04: visualization both approaches ----

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

best_two_videos_hough <- df %>% 
  dplyr::filter(videoname == "June16_20" | videoname == "June16_02")

June16_20_hough_exp4 <- joined_radii_all %>% 
  dplyr::filter(examplenr == "exp4") %>% 
  dplyr::filter(videoname == "June16_20")

June16_02_hough_exp4 <- joined_radii_all %>% 
  dplyr::filter(examplenr == "exp4") %>% 
  dplyr::filter(videoname == "June16_02") 

# plotting radius comparison

# visualization DLC
dlc <- joined_radii_dlc %>%
  ggplot(aes(x= radius_man, y = radius))+
  geom_point(aes( fill = videoname, color = videoname), size = 2, alpha = 0.6)+
  geom_smooth(method = 'lm', color = "grey15")+
  scale_color_manual(values = safe_colorblind_palette)+
  ylab('Automatically Tracked Radius [px], DLC')+
  xlab('Manually Labeled Radius [px]')+
  scale_y_continuous(limits = c(75, 270),
                     breaks = c(50, 100, 150, 200, 250))+
  scale_x_continuous(limits = c(100, 270),
                     breaks = c(100, 150, 200, 250))+
  #scale_x_continuous(limits = c(120, 270),
  #                   breaks = c( 120, 160, 200, 240))+
  theme_minimal()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  annotate("text", x=200, y= 75, label= " R² = 0.86", size = 5)


# visualization hough

hough <- joined_radii_all %>%
  dplyr::filter(examplenr == "exp4") %>% 
  dplyr::filter(radius_man >= 100) %>% 
  ggplot(aes(x= radius_man, y = smoothed_hough_radius_kolmogorov))+
  geom_point(aes( fill = videoname, color = videoname), size = 2, alpha = 0.6)+
  geom_smooth(method = 'lm', color = "grey15")+
  geom_smooth(data = June16_02_hough_exp4, method = "lm", color = "#44AA99")+
  geom_smooth(data = June16_20_hough_exp4, method = "lm", color = "#882255")+
  scale_color_manual(values = safe_colorblind_palette)+
  ylab('Automatically Tracked Radius [px], Hough')+
  xlab('Manually Labeled Radius [px]')+
  scale_y_continuous(limits = c(75, 270),
                     breaks = c(50, 100, 150, 200, 250))+
  scale_x_continuous(limits = c(100, 270),
                     breaks = c(100, 150, 200, 250))+
  #scale_x_continuous(limits = c(120, 270),
  #                   breaks = c( 120, 160, 200, 240))+
  theme_minimal()+
  theme(text = element_text(size = 20),
        legend.position="none")+
  #guides(fill = FALSE)+
  annotate("text", x=200, y= 75, label= " R² = 0.23", size = 5)+
  annotate("text", x=200, y= 85, label= " R² = 0.80", size = 5, color = "#44AA99")+
  annotate("text", x=200, y= 95, label= " R² = 0.53", size = 5, color = "#882255")

cowplot::plot_grid(dlc, hough, ncol = 2, labels = c("A", "B"), label_size = 16)

ggsave("comparison_dlc_hough_manual_new.jpg", dpi = 300, width= 10, height = 8)
