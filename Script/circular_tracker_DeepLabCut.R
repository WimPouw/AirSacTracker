# Airsac Tracker Project
# DeepLabCut approach (DLC)
# estimate circles (r and center) from >=three points estimated with DLC
# 


# version one: analytical
# second version: Least-squares fitting of circles
# https://lucidar.me/en/mathematics/least-squares-fitting-of-circle/ 
# https://stackoverflow.com/questions/67058935/least-squares-regression-to-fit-circle-with-constrained-centre-point
# https://people.cas.uab.edu/~mosya/cl/cl1.pdf 
# https://cran.r-project.org/web/packages/conicfit/conicfit.pdf 
# geometric fit method, is well known too apparently
# Levenberg-Marquardt for non-linear OLS

# 00: load packages -------------------------------------

library(tidyverse)

# 01: load data -----------------------------------------

df_all <- read_delim("AirSacTracker/DeepLabCut/evaluation_set/June09_01DLC_resnet101_Deep_AirSacTrackingV1Jan1shuffle1_280000_filtered.csv", delim = ',')

videoname <- "June09_01"
# redefine column names to distinct names, delete first two rows with identifiers afterwards
colnames <- c()
colnames[1] <- "frames"
for (x in 2: ncol(df_all)){
  colnames[x] <- paste(df_all[1,x], df_all[2,x], sep = '_')
}

colnames(df_all) <- colnames
df_all <- df_all[-1:-2,]


# 02: define functions for circle estimation ------------

# https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points
# https://stackoverflow.com/questions/55357724/plot-circle-segment-defined-by-three-points-with-ggplot2 

# question: does it work as well, when df contains more than 3 points? results better or worse?

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

# define df for radius estimation and estimate ------------------

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

radius_results <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("videoname", "frame","radius", "center_x", "center_y")
colnames(radius_results) <- x

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
    filter(likelihood >= threshold)
  
# circle estimation and saving of 
  
if(nrow(df_filter) >= 3){
  
 circles <- get_circle(df_filter[1:3,])
 
 radius_results[a,1] <- videoname
 radius_results[a,2] <- a
 radius_results[a,3] <- circles$r
 radius_results[a,4] <- circles$center[1]
 radius_results[a,5] <- circles$center[2]
 
} else {
  radius_results[a,1] <- videoname
  radius_results[a,2] <- a
  radius_results[a,3] <- NA
  radius_results[a,4] <- NA
  radius_results[a,5] <- NA
  
        }
}

