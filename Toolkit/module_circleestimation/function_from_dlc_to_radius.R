# functions to get circle (i.e. radius) estimations from DLC trackings 
# automatically estimate circles for all DLC tracking datafiles
# of interest for subsequent analysis and comparison of radii with acoustic parameters
# with Landau circle estimation

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

install_load("tidyverse","conicfit", "scales", "spiro", "signal", "foreach")
library("foreach")
# 01a: load data ----

# path in Github structure depends on dataset
# For individual boom snippets: \AirSacTracker\Project\Proof_of_concept_analyses\snippets\tracked_dlc
# For boom - bark - sequences: \AirSacTracker\Project\Proof_of_concept_analyses2\snippets\short_sequences\tracked\boom

path <- choose.dir()
# path <- getwd() 
pattern <- "csv"
list_of_files <- list.files(path = path, pattern = pattern)

# 01b: set parameters, dataframes ----

# DLC tracking likelihood threshold

threshold <- 0.6

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

## sub function normalization ----

nose_eye_normalization <- function(auto_data, a = a,  min_frames = 2, threshold_normalization = 0.8){  
  norm_data <- c()
  list_normalization_points <- c('Nose_x','Nose_y', 'Nose_likelihood',
                                 'EyeBridge_x', 'EyeBridge_y', 'EyeBridge_likelihood')
  
  colnames <- c()
  colnames[1] <- "frames"
  for (x in 2: ncol(auto_data)){
    colnames[x] <- paste(auto_data[1,x],auto_data[2,x], sep = '_')
  }
  
  colnames(auto_data) <- colnames
  df <- auto_data[-1:-2,]
  
  df_sub <- df %>%
    select(all_of(list_normalization_points))
  
  df_sub <- as.data.frame(apply(df_sub,2,as.numeric))
  
  euc_dist <- function(xbridge, xnose, ybridge, ynose){
    return(sqrt(((xbridge-xnose)^2+((ybridge-ynose)^2))))
  }
  
  df_sub_normalization <- df_sub %>% 
    dplyr::filter(Nose_likelihood >= threshold_normalization && 
                    EyeBridge_likelihood >= threshold_normalization)
  
  if(nrow(df_sub_normalization) >= min_frames){  
    
    distance <- foreach::foreach(i= 1:nrow(df_sub_normalization), .combine = c) %do% 
      euc_dist(df_sub_normalization$EyeBridge_x, df_sub_normalization$Nose_x,
               df_sub_normalization$EyeBridge_y, df_sub_normalization$Nose_y)
  } else {
    distance <- NA
  }

  norm_data <- mean(distance, na.rm = TRUE)
  return(norm_data)
}

## complete function ----

from_DLC_to_circle <- function(path, list_of_files){
  
  radius_all <- data.frame(matrix(ncol = 5, nrow = 0))
  z <- c("radius","x","y", "frame", "videofile")
  colnames(radius_all) <- z
  
  normalization_value <- data.frame(matrix(ncol = 2, nrow = length(list_of_files)))
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
      
      if (length(grouped_data_circle_estimation) > 0){
      frame_data <- as.data.frame(grouped_data_circle_estimation[[n]])
      
      if(nrow(frame_data)>=3){
        
        circles_LAN <- CircleFitByLandau(frame_data[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
        
      }else {
        circles_LAN <- c(NA, NA, NA, NA, NA)
      }  
      circles_res <- c(circles_LAN[3],circles_LAN[1], circles_LAN[2], n, list_of_files[a])
      
      radius <- rbind(radius, circles_res)
      colnames(radius) <- z
      } else {circles_LAN <- c(NA, NA, NA, NA, NA)
    
      circles_res <- c(circles_LAN[3],circles_LAN[1], circles_LAN[2], n, list_of_files[a])
      
      radius <- rbind(radius, circles_res)
      colnames(radius) <- z}
    }
    print(a)
    radius_all <- rbind(radius_all, radius)
    
    normalization_value$normalization_value[a] <- nose_eye_normalization(auto_data)
    normalization_value$videofile[a] <- list_of_files[a]
    
  }
  
  results_I <- list(radius_all, normalization_value)
  results <- left_join(results_I[[1]], results_I[[2]], by = "videofile")
  
  results$norm_radius <- as.numeric(results$radius)/results$normalization_value
  
  return(results)
}

# 02: running analysis -----

results <- from_DLC_to_circle(path = path, list_of_files = list_of_files)



# 03: saving ----
#save to csv
csv_savename <- readline(prompt = "Enter a savename for the dataset for the csv file without any quote signs, the phrase 'DLCtoRadii' will be added at the end:")

write.table(results, paste0(path, '/', csv_savename, '_DLCtoRadii.csv'), row.names = FALSE, sep = ",")

savename <- readline(prompt = "Enter a savename for the dataset, including the fileending .rds but without any quote signs:")

saveRDS(results, file = savename)
