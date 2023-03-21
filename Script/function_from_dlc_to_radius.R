# batch processing: automatically estimate circles for all DLC tracking datafiles
# of interest for subsequent analysis and comparison of radii with acoustic parameters

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: folder with DLC trackings in csv format
# end point: circle estimations per frames of all videos chosen, structure: radius, frame, videofile

############

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse","conicfit", "scales", "spiro", "signal")

# 01a: load data ----

path <- choose.dir()
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

# 01b: set parameters, dataframes ----

# DLC tracking likelihood threshold
threshold <- 0.6
threshold_normalization <- 0.9 

# 01c: functions ----------------
## sub function data preparation ----

data_prep_radius_estim_DLC <- function(data){
  
  # data: data to process, a dataframe
  
  ############
  
  # 00: load packages ----
  
  if (!require(install.load)) {
    install.packages("install.load")
  }
  
  library(install.load)
  
  install_load("tidyverse")
  
  # 01: functionality ----
  
  ## 01a : sub-function circle_format ----
  # sub-function to transform already pre-processed datasheet into format
  # necessary for circle estimation, which needs x and y for all points in two
  # columns, 1 point combination per row
  
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
  
  ## 01b: main ----
  
  # list of columns needed for circle estimation used later in function
  list_airsac_points <- c('Start_outline_outer_left_x','Start_outline_outer_left_y', 'Start_outline_outer_left_likelihood',
                          'Start_outline_outer_right_x', 'Start_outline_outer_right_y', 'Start_outline_outer_right_likelihood',
                          'LowestPoint_outline_x', 'LowestPoint_outline_y', 'LowestPoint_outline_likelihood',
                          'MidLowleft_outline_x', 'MidLowleft_outline_y', 'MidLowleft_outline_likelihood',
                          'MidLowright_outline_x', 'MidLowright_outline_y', 'MidLowright_outline_likelihood')
  
  colnames <- c()
  colnames[1] <- "frames"
  for (x in 2: ncol(data)){
    colnames[x] <- paste(data[1,x],data[2,x], sep = '_')
  }
  
  colnames(data) <- colnames
  df_all <- data[-1:-2,]
  
  df_sub <- df_all %>%
    select(all_of(list_airsac_points))
  
  df_sub <- as.data.frame(apply(df_sub,2,as.numeric))
  
  circle_format_data <- circle_format(df_sub)
  
  return(circle_format_data)
}

## complete function ----
from_DLC_to_circle <- function(path, list_of_files){
  
  radius_all <- data.frame(matrix(ncol = 3, nrow = 0))
  z <- c("radius", "frame", "videofile")
  colnames(radius_all) <- z
  
  normalization_value <- data.frame(matrix(ncol = 2, nrow = 0))
  y <- c("normalization_value", "videofile")
  colnames(normalization_value) <- y
  ## starting main loop ----
  for (a in 1:length(list_of_files)) {
    
    radius <- data.frame()
    
    auto_data <- read_delim(paste(path, list_of_files[a], sep = "\\"), delim = "," )
    
    # 02: pre-processing data ----
    
    data_circle_estimation <- data_prep_radius_estim_DLC(data = auto_data)
    grouped_data_circle_estimation <- data_circle_estimation %>%
      dplyr::filter(likelihood > threshold) %>% 
      group_split(frame)
    
    for(n in 1:length(grouped_data_circle_estimation)){
      
      frame_data <- as.data.frame(grouped_data_circle_estimation[[n]])
      
      if(nrow(frame_data)>=3){
        
        circles_LAN <- CircleFitByLandau(frame_data[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
        
      }else {
        circles_LAN <- c(NA, NA, NA)
      }  
      circles_res <- c(circles_LAN[3],n, list_of_files[a])
      
      radius <- rbind(radius, circles_res)
      colnames(radius) <- z
    }
    
    radius_all <- rbind(radius_all, radius)
    
  }
  return(radius_all)
}

# 02: running analysis -----

results <- from_DLC_to_circle(path = path, list_of_files = list_of_files)