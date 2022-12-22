# Intro ------------------------------------
# circular tracking, data analysis
# author: Lara B.
# to analyse and compare the manually tracked data of the Siamang air sacs
# and the automatically tracked data with the python airsac tracker
# 

# 00: packages & functions ----------------------

library(tidyverse)
library(effsize)
library(psych)
library(signal)

butter.it <- function(x, samplingrate, order, lowpasscutoff)
{bf <- butter(order,lowpasscutoff/(samplingrate/2), type="low") #normalized frequency
x <<- as.numeric(signal::filtfilt(bf, x))} 

# 01: loading data --------------------------

#auto_data <- read_delim("sens_analysis_random_comb_tracked_06092022_final.csv", delim = ",")
#auto_data <- read_delim("AirSacTracker/TestResults/c1_8_c2_13_al_2.5_b_35_dil_7_blur_5_vid_June09_01.mp4.csv", delim = ";")

man_data <- read_delim("manual_results_circular_tracking_deleted_duplicate_second_version.csv", delim = ";")

#add transformation from perimeter to radius in manual data
man_data$...6 <- man_data$Perim./(2*pi)
names(man_data)[names(man_data) == "...6"] <- "radius_man"
#names(man_data)[names(man_data) == "Label"] <- "input_id"
names(man_data)[names(man_data) == "Label"] <- "name"


# set path for comparison
path <- choose.dir()
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

# set up results dataframe

df <- data.frame()

# 02: main -------------------------------------------
# start loop

for (a in 1:length(list_of_files)) {
  
auto_data <- read_delim(paste(path, list_of_files[a], sep = "\\"), delim = "," )

#getting videoname
videoname <- list_of_files[a]
videoname <- substring(videoname, 1, nchar(videoname)-8)
videoname <- str_split(videoname, pattern = "_")
videoname <- videoname[[1]]
videoname <- paste(videoname[14], videoname[15], sep = "_")

#only necessary for wrong files, where name within file was saved incorrecttly
# switching wrong part with right part, wrong part is "June16_02", right part is videoname

# this needs to be in loop, because we also need the frame
# also substring should not be used with character index but again with seperator, because number format changess
correct_name <- data.frame()
for (b in 1:nrow(auto_data)){
  
match_name <- auto_data$name[b]
match_name <- str_split(match_name, pattern = "_")
match_name <- match_name[[1]]
match_name <- paste(match_name[1],match_name[2], match_name[3], videoname, sep = "_")
#match_name <- paste(substring(match_name, 1, 19), videoname, sep = "")

correct_name[b,1] <- match_name
}

auto_data <- auto_data %>% 
  mutate(name = correct_name$V1)

#smoothing of auto tracked radius, adding to auto_data
# smoothing: butterworth filter, function defined in 00

smoothed_auto <- as.data.frame(
                  butter.it(auto_data$r, samplingrate =  25, order = 2, lowpasscutoff = 15))

auto_data <- auto_data %>% 
  mutate(smoothed_r = smoothed_auto[1:nrow(smoothed_auto), 1])

# joining datasets
joined_data <- left_join(auto_data, man_data, by = "name")

# removing all lines, where no circle was tracked manually

joined_data <- joined_data %>% 
  filter(radius_man > 100 | !is.na(radius_man)) %>% 
  filter(r < 250)

# run correlation

cor_r <- corr.test(joined_data$r, joined_data$radius_man)
cor_r_sm <- corr.test(joined_data$smoothed_r, joined_data$radius_man)
cor_x <- corr.test(joined_data$x, joined_data$X)
cor_y <- corr.test(joined_data$y, joined_data$Y)

# difference in tracking
diff <- c()
diff_sqr <- c() 
for (c in 1: nrow(joined_data)){
  
  diff[c] <- round(joined_data$r[c] - joined_data$radius_man[c], digits = 2)
  diff_sqr[c] <- diff[c]^2
  
}

diff_min <- min(abs(diff))
diff_max <- max(abs(diff))
diff_avg <- mean(abs(diff))
diff_median <- median(abs(diff))

diff_sqr_min <- min(abs(diff_sqr))
diff_sqr_max <- max(abs(diff_sqr))
diff_sqr_avg <- mean(abs(diff_sqr))
diff_sqr_median <- median(abs(diff_sqr))

# save relevant parameters 
# what do we need to save? 
# video name, correlation, preprocessing parameters, how many frames in correlation from how many in the video

df[a,1] <- videoname
df[a,2] <- cor_r$r
df[a,3] <- cor_r_sm$r
df[a,4] <- cor_x$r
df[a,5] <- cor_y$r
df[a,6] <- nrow(joined_data)
df[a,7] <- nrow(auto_data)
df[a,8] <- joined_data$alpha[1]
df[a,9] <- joined_data$beta[1]
df[a,10] <- joined_data$threh_div1[1]
df[a,11] <- joined_data$threh_div2[1]
df[a,12] <- joined_data$dilation[1]
df[a,13] <- joined_data$phase1_medianblur[1]
df[a,14] <- diff_min 
df[a,15] <- diff_max
df[a,16] <- diff_avg 
df[a,17] <- diff_median
df[a,18] <- diff_sqr_min 
df[a,19] <- diff_sqr_max
df[a,20] <- diff_sqr_avg
df[a,21] <- diff_sqr_median 



}

colnames(df) <- c("videoname", "cor_radius","cor_radius_smoothed", "cor_x",  "cor_y", "nr_frames_used", "nr_frames_in_video",
                               "alpha","beta","thresh_div1","thresh_div2","dilation",
                               "medianblur", "diff_r_min", "diff_r_max","diff_r_avg","diff_r_median","diff_r_sqr_min",
                                "diff_r_sqr_max","diff_r_sqr_avg","diff_r_sqr_median")

write.csv(df, 'comparison_man_auto_tracker.csv', row.names = FALSE)

#### end of comparison loop 

# 03: comparison statistics --------------

# compare correlation of smoothed radius per video

#sm_r_per_video_low_15 <- aggregate(cor_radius_smoothed ~ videoname,data = df, FUN = median)

best_combos_low_5_video <- df_low_5 %>% 
  group_by(videoname) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 

best_combos_low_10_video <- df_low_10 %>% 
  group_by(videoname) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 

best_combos_low_15_video <- df_low_15 %>% 
  group_by(videoname) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 


#write.csv(sm_r_videos, 'correlations_different_filter.csv', row.names = FALSE)



# aggregate per parameter setting

best_combos_low_5 <- df_low_5 %>% 
  group_by(alpha, beat, thresh_div1, thresh_div2, dilation, medianblur) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 

best_combos_low_10 <- df_low_10 %>% 
  group_by(alpha, beat, thresh_div1, thresh_div2, dilation, medianblur) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 

best_combos_low_15 <- df_low_10 %>% 
  group_by(alpha, beat, thresh_div1, thresh_div2, dilation, medianblur) %>%                            # multiple group columns
  summarise(mean_cor = mean(cor_radius_smoothed), median_cor = median(cor_radius_smoothed),
            min_cor = min(cor_radius_smoothed), max_cor = max(cor_radius_smoothed)) 

###############################################################################
# 04: some other stuff ----------------------------
joined_data %>% 
  ggplot(aes(x=smoothed_r, y = radius_man ))+
  geom_point()+
  geom_smooth(method = lm)

#does video length influence tracking success?

df_low_5 %>% 
  ggplot(aes(x= nr_frames_in_video, y = cor_radius_smoothed, fill = videoname, color = videoname))+
  geom_jitter()

df_low_5 %>% 
  ggplot(aes(x= nr_frames_used, y = cor_radius_smoothed, fill = videoname, color = videoname))+
  geom_jitter()
# add differences in radius tracked

joined_data$diff <- round(joined_data$r - joined_data$radius, digits = 2)
joined_data$diff_sqr <- joined_data$diff^2


# group per parameter combo and calculate some statistcis

joined_data %>% 
  group_by(alpha, beta, canny_min, canny_max, erode_it_1, erode_it_2, param2) %>% 
  ggplot(aes(y = diff_sqr))+
  geom_boxplot()
  
  
joined_data_stats <- joined_data %>% 
  group_by(alpha, beta, canny_min, canny_max, erode_it_1, erode_it_2, param2) %>%
  summarise_at(vars(diff), list( ~ mean(., na.rm = TRUE), ~median(., na.rm = TRUE),
                                ~ min(., na.rm = TRUE), ~max(., na.rm = TRUE),
                                ~sd(., na.rm = TRUE) , ~cor(r, radius, method = "pearson")))



# plot data with highest correlation on r and radius

joined_data %>% 
  filter(alpha == 2, beta == 20, canny_min == 20, canny_max == 30, erode_it_1 == 10,
         erode_it_2 == 5, param2 == 22) %>% 
  ggplot()+
  geom_point(aes(x = radius, y = r))+
  xlim(c(100, 250))+
  ylim(c(100, 250))


# simple linear model 

#on scaled data? that induces errors and issues, I don't understand so maybe not

auto_data_scaled <- auto_data %>% 
  select(-"inputfile") %>% 
  scale() %>% 
  as.data.frame()

# 
model <- lm(r ~ alpha  +  param2 + erode_it_1 , auto_data)
summary(model)

# look at NA values for tracking

na_radius <- auto_data %>% 
  filter(is.na(r))
