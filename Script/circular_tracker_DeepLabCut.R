# Airsac Tracker Project
# DeepLabCut approach (DLC)
# estimate circles (r and center) from >=three points estimated with DLC
# 


# version one: analytical
# second version: Least-squares fitting of circles
# https://doi.org/10.1007/s10851-005-0482-8 : Least Square fitting of circle, Chernov & Lesort 2005
# https://lucidar.me/en/mathematics/least-squares-fitting-of-circle/ 
# https://stackoverflow.com/questions/67058935/least-squares-regression-to-fit-circle-with-constrained-centre-point
# https://people.cas.uab.edu/~mosya/cl/cl1.pdf 
# https://cran.r-project.org/web/packages/conicfit/conicfit.pdf 
# geometric fit method, is well known too apparently
# Levenberg-Marquardt for non-linear OLS


# code to test difference in calculation time
# start.time <- Sys.time()
# ...Relevent codes...
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

# 00: load packages -------------------------------------

library(tidyverse)
library(conicfit)
library(psych)
library(signal)
library(zoo)
library(scales)

# 01: load data -----------------------------------------

# set path for comparison
path <- choose.dir()
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

#df_all <- read_delim("AirSacTracker/DeepLabCut/evaluation_set/June09_01DLC_resnet101_Deep_AirSacTrackingV1Jan1shuffle1_280000_filtered.csv", delim = ',')
#videoname <- "June09_01"

# load manually tracked data, as combined in script "circular_tracker_snippets"
man_data <- read_delim("manually_tracked_airsac_radii.csv", delim = ',')


# 02: define functions  ------------

# https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points
# https://stackoverflow.com/questions/55357724/plot-circle-segment-defined-by-three-points-with-ggplot2 

# this approach only works with exactly three points

### get circles using linear algebra function
get_circle <- function(df){
  # df: three-row data frame containing columns x and y
  mat <- 
    df %>% 
    transmute(ss = x^2 + y^2, x, y, ones = 1) %>% 
    as.matrix
  
  center <- 
    c(x = det(mat[,c('ss', 'y', 'ones')]), y = -det(mat[,c('ss', 'x', 'ones')])
    )/(2*det(mat[,c('x', 'y', 'ones')]))
  
  r <- sqrt(sum((unlist(df[1, c('x', 'y')]) - center)^2))
  
  list(center = center, r = r)
}

### butterworth filter function 
butter_it <- function(x, samplingrate =  25, order = 2, lowpasscutoff = 10)
{bf <- butter(order,lowpasscutoff/(samplingrate/2), type="low") #normalized frequency by the nyquist limit
x <<- as.numeric(signal::filtfilt(bf, x))} 

# start main loop ------------------------------------------------

comparison_all <- data.frame()

for (i in 1:length(list_of_files)) {
  
auto_data <- read_delim(paste(path, list_of_files[i], sep = "\\"), delim = "," )

videoname <- substring(list_of_files[i], 1, 9)

# redefine column names to distinct names, delete first two rows with identifiers afterwards
colnames <- c()
colnames[1] <- "frames"
for (x in 2: ncol(auto_data)){
  colnames[x] <- paste(auto_data[1,x], auto_data[2,x], sep = '_')
}

colnames(auto_data) <- colnames
df_all <- auto_data[-1:-2,]


# 03: circle transformation -------------------------------
## 03a: data preparation  -----------------------------------

# define df for radius estimation and estimate 

# we need a clean df per frame with two columns x and y
# with at least three rows, the points we tracked with DLC
# so we need to loop through df_all, one row is one frame
# per frame we need the five tracked points of the airsac with x,y and likelihood
# when we have that, we filter for likelihood
# then we calculate the circle radius and center from the remaining points (at least 3)
# if too few points remain after likelihood condition, we assign NA to that frame

threshold <- 0.6 #threshold for likelihood for DLC tracking of points 

df <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("x", "y", "likelihood")
colnames(df) <- x

radius_results_lan <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("videoname", "frame","radius_lan")
colnames(radius_results_lan) <- x

for (a in 1:nrow(df_all)){

#subset DLC results to necessary columns  
    
  df_sub <- df_all %>% 
    select(Start_outline_outer_left_x,Start_outline_outer_left_y, Start_outline_outer_left_likelihood,
           Start_outline_outer_right_x, Start_outline_outer_right_y, Start_outline_outer_right_likelihood,
           LowestPoint_outline_x, LowestPoint_outline_y, LowestPoint_outline_likelihood,
           MidLowleft_outline_x, MidLowleft_outline_y, MidLowleft_outline_likelihood,
           MidLowright_outline_x, MidLowright_outline_y, MidLowright_outline_likelihood)
  
  df_sub <- as.data.frame(apply(df_sub,2,as.numeric))
  
# save DLC points in necessary format to estimate circles
  
for (b in 1:5){
  end_col <- b*3
  start_col <- end_col - 2
  
  df[b, 1:3] <-  df_sub[a, start_col:end_col]
    } 

# check for likelihood
# do we need to make sure, to somehow save information, which points are used? 
  
  df_filter <- df %>% 
    dplyr::filter(likelihood >= threshold)

  
## 03b: NA interpolation  

    
## 03c: linear algebra approach --------------------  
# circle estimation and saving of results
# using function defined in 02, using matrices and determinants

if(nrow(df_filter) >= 3){
  
 ## 03c: Least-squares fitting of circles--------------------------
 
 #using package conifcit, Landau algorithm 
 
 circles_LAN <- CircleFitByLandau(df_filter[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
 
 radius_results_lan[a,1] <- videoname
 radius_results_lan[a,2] <- a
 radius_results_lan[a,3] <- circles_LAN[3]
 radius_results_lan[a,4] <- circles_LAN[1]
 radius_results_lan[a,5] <- circles_LAN[2]
 
} else {
  
  radius_results_lan[a,1] <- videoname
  radius_results_lan[a,2] <- a
  radius_results_lan[a,3] <- NA
  radius_results_lan[a,4] <- NA
  radius_results_lan[a,5] <- NA
  
        }
}

comparison_method <- radius_results_lan

comparison_all <- rbind(comparison_all, comparison_method)

}

## 03c: interpolating NAs in DLC data introduced by likelihood criterion
# not useful, because doesn't change the data used in correlation because of filters in lines 301 and 303
#comparison_all <- comparison_all[1:1637,] #there was an NA in the last line, that na.approx couldn't work with

#comparison_all <- comparison_all %>% 
#  mutate(radius_lina_interpol = na.approx(comparison_all$radius_lina, maxgap = 4),
#         radius_lan_interpol = na.approx(comparison_all$radius_lan, maxgap = 4),
#         radius_lm_interpol = na.approx(comparison_all$radius_lm, maxgap = 4)) 

# need to add to be able to use join, to join with manual data: framenr_1_framevid_June16_20
match_name <- c()

for (i in 1:nrow(comparison_all)){
  match_name[i] <- paste("framenr", comparison_all[i,2], "framevid", comparison_all[i,1], sep = "_")
}

col_nr <- ncol(comparison_all)+1

comparison_all[,col_nr] <- match_name
colnames(comparison_all)[col_nr] <- "name"


# 03d: comparison to manually tracked data---------------------------------------

joined_data_DLC <- left_join(comparison_all, man_data, by = "name")

joined_data_DLC_rounded <- joined_data_DLC %>% 
  mutate(radius_lan = round(radius_lan, 0))

joined_data_DLC_rounded <- joined_data_DLC_rounded %>% 
  dplyr::filter(radius_man > 100) %>% 
  dplyr::filter( radius_lan < 400)

cor_r_lan <- corr.test(joined_data_DLC_rounded$radius_lan, joined_data_DLC_rounded$radius_man)
cor_x_lan <- corr.test(joined_data_DLC_rounded$center_x_lan, joined_data_DLC_rounded$X)
cor_y_lan <- corr.test(joined_data_DLC_rounded$center_y_lan, joined_data_DLC_rounded$Y)

cor_r_lan
cor_x_lan
cor_y_lan


# 04: smoothing -----------------------------------
comparison_all$time <- comparison_all$frame*(1/25)
joined_data_DLC_rounded$time <- joined_data_DLC_rounded$frame*(1/25)

#auto_data$smoothed_r <- butter.it(auto_data$r, samplingrate =  25, order = 4, lowpasscutoff = 10)
#auto_data$smoothed_r_euro <-  OneEurofiltfilt(t = auto_data$time, x = auto_data$r)


smoothed_radii <- data.frame(matrix(ncol = 3, nrow = nrow(joined_data_DLC_rounded)))
  
smoothed_radii[,1] <- butter.it(joined_data_DLC_rounded$radius_lan, samplingrate =  25, order = 4, lowpasscutoff = 10)

colnames(smoothed_radii) <- c("sm_butter_r_lan")

# correlation smoothing

# for some reason, there are very high smoothed values in the filtered data
# so we filter again (which actually shouldn't be necessary?)


joined_data_DLC_sm <- cbind(joined_data_DLC_rounded, smoothed_radii)


# butterworth smoothing correlation

cor_sm_butter_lan_r <- corr.test(joined_data_DLC_sm$radius_man, joined_data_DLC_sm$sm_butter_r_lan)

cor_sm_butter_lan_r 

# butterworth filter works much better on the data, order = 2, low cutoff= 10

# 05: normalization ----

comparison_data <- readRDS('radius_acoustic_param_comparison_2.rds')

# we run two different normalizations: rescaling radius to ly between 0 and 1
# and standardizing by nose bridge (only for suitable subset, calcualte nosebridge distance using euclidean distance
# in two dimensions and dividing every radius by this value)

## 05a: rescaling 0-1 ----

comparison_data_scaled <- comparison_data %>% 
  group_by(audiofile) %>% 
  mutate(radius_scaled = scales::rescale(radius, to = c(0,1)))

saveRDS(comparison_data_scaled, file = "radius_acoustic_param_comparison_scaled.rds")

## 05b: standadrize by nose bridge ----


