# function to estimate circle radii from DeepLabCut trackings
# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: DLC data sheet
# end point: circle estimation per frame
# used in the following scripts: 

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
