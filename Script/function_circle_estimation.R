# circle estimation as an individual function

# start point: DLC data sheet
# end point: circle estimation per frame

# try to avoid loops, use apply as often as possible

############

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

# 01: load data ------------

path <- choose.dir()
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

# 02: function

# data: data to process
# threshold: likelihood threshold (DLC)
# frame: row number to process in datasubset -> the frame to estimate circles for

circle_format <- function(df_sub){
  
  df_all <- data.frame()
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("x", "y", "likelihood", "frame")
  colnames(df) <- x
    
for (frame in 1:nrow(df_sub)){  
  for (b in 1:5){
    end_col <- b*3
    start_col <- end_col - 2
    
    df[b, 1:3] <-  df_sub[frame, start_col:end_col]
    df[b, 4] <- frame
  
  }   
  
  df_all <- rbind(df_all, df)
  
  }
  return(df_all)
}
  
# version 2, circle_format, can we use apply instead of outer loop?

# circle_format <- function(df_sub){
#   
#   df_all <- data.frame()
#   df <- data.frame(matrix(ncol = 4, nrow = 0))
#   x <- c("x", "y", "likelihood", "frame")
#   colnames(df) <- x
#   
#   df_all <- sapply(df_sub, MARGIN = 1, FUN = circle_format_1(df_sub))
#     
#     df_all <- rbind(df_all, df)
#     
#   return(df_all)
# }
# 
# circle_format_1 <- function(df_sub){
#   
#   for (b in 1:5){
#     end_col <- b*3
#     start_col <- end_col - 2
#     
#     df[b, 1:3] <-  df_sub[frame, start_col:end_col]
#     df[b, 4] <- frame
#   }
#   return(df)
# }

# main function data preparation
  
data_prep_radius_estim_DLC <- function(data, threshold, frame){
  
  colnames <- c()
  colnames[1] <- "frames"
  for (x in 2: ncol(data)){
    colnames[x] <- paste(data[1,x],data[2,x], sep = '_')
  }
  
  colnames(data) <- colnames
  df_all <- data[-1:-2,]
  
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("x", "y", "likelihood")
  colnames(df) <- x
  
  df_sub <- df_all %>% 
    select(Start_outline_outer_left_x,Start_outline_outer_left_y, Start_outline_outer_left_likelihood,
           Start_outline_outer_right_x, Start_outline_outer_right_y, Start_outline_outer_right_likelihood,
           LowestPoint_outline_x, LowestPoint_outline_y, LowestPoint_outline_likelihood,
           MidLowleft_outline_x, MidLowleft_outline_y, MidLowleft_outline_likelihood,
           MidLowright_outline_x, MidLowright_outline_y, MidLowright_outline_likelihood)
  
  df_sub <- as.data.frame(apply(df_sub,2,as.numeric))
  
  # save DLC points in necessary format to estimate circles
  # circle estimation needs x as one column and y as second column, we save likelihood as third column
  # to be able to filter for that later 
  # b in 1:5 -> 5 potentially tracked points in circle
  # every row of df_sub is a frame, always three columns per point (x,y,likelihood)
  
# rewrite this part as a function, so that we can use apply to use this function on every row of df_sub
# then saving all frames in same datafile is also easier
circle_format_data <- circle_format(df_sub)
  
  
  return(circle_format_data)  #which variables need to be returned by the function? 
}


# 03: analysis ----


comparison_all <- data.frame()

for (i in 1:length(list_of_files)) {
  
  auto_data <- read_delim(paste(path, list_of_files[i], sep = "\\"), delim = "," )
  
  videoname <- substring(list_of_files[i], 6, 14) #needs to be checked individually, depends on filenames
  
  circle_format_data <- data_prep_radius_estim_DLC(auto_data, threshold = 0.6)
  
  list_circles <- circle_format_data %>% 
  group_split(frame)
  
}

  df_filter <- df %>% 
    dplyr::filter(likelihood >= threshold)
  
  # after function:
  # possible with apply? can we have the save DLC points in a continous datasheet with id column to group? then we could use group walk
  # for the circle estimation? 
  if(nrow(df_filter) >= 3){
    
    circles_LAN <- CircleFitByLandau(df_filter[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
    
  }