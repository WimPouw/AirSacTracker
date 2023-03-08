# batch processing: automatically estimate circles for all DLC tracking datafiles
# of interest for subsequent analysis and comparison of radii with acoustic parameters

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: folder with DLC trackings in csv format
# end point: circle estimations for all frames of all files chosen

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

# 03: circle estimation ----
  
## 03a: nose-eyebridge normalization ----
  # to be run on the subset of data, suitable for nose-eyebridge normalization! 
  # comment out for full run
  
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
  
  euc_dist <- function(x1, x2){
    return(sqrt(sum((x1 - x2)^2)))
  }
  
  df_sub_normalization_nose <- df_sub %>%
    dplyr::filter(Nose_likelihood >= threshold) %>% 
    select(Nose_x,Nose_y) %>% 
    transmute(Nose_x = as.numeric(Nose_x),
              Nose_y = as.numeric(Nose_y))
  
  df_sub_normalization_eyebridge <- df_all %>% 
    dplyr::filter(EyeBridge_likelihood > threshold) %>% 
    select(EyeBridge_x, EyeBridge_y) %>% 
    transmute(EyeBridge_x = as.numeric(EyeBridge_x),
              EyeBridge_y = as.numeric(EyeBridge_y))
  
  
  distances <- foreach(i = 1:nrow(df_sub_normalization_nose), .combine = c) %do% 
    euc_dist(df_sub_normalization_nose[i,], df_sub_normalization_eyebridge[i,])
  
  normalization_value[a,] <- c(mean(distances, na.rm = TRUE), list_of_files[a])  

# 03b: circle estimations
grouped_data_circle_estimation <- data_circle_estimation %>%
  dplyr::filter(likelihood > threshold) %>% 
  group_split(frame)

    for(n in 1:length(grouped_data_circle_estimation)){
      
      frame_data <- as.data.frame(grouped_data_circle_estimation[[n]])
      
      if(nrow(frame_data)>=3){
      
      circels_LAN <- CircleFitByLandau(frame_data[,1:2], ParIni = NA, epsilon = 1e-06, IterMAX = 500)
      
      }else {
        circles_LAN <- c(NA, NA, NA)
      }  
        circles_res <- c(circels_LAN[3],n, list_of_files[a])
        
        radius <- rbind(radius, circles_res)
        colnames(radius) <- z
      }
      

  
      radius_all <- rbind(radius_all, radius)
      
}


# 04a: post-processing I: adding meta-data, sex + ID + normalization value (if applicable) ----

# adding column with information on sex
index_female <- grep('Pelangi', radius_all$videofile)
radius_all$index <- 1:nrow(radius_all)

for(a in 1: nrow(radius_all)){
  if(radius_all$index[a] %in% index_female == TRUE){
    
    radius_all$sex[a] = 'f'
    
  } else {
    
    radius_all$sex[a] = 'm'
    
  }
}

# names: fajar, baju, roger, pelangi
for (b in 1:nrow(radius_all)){
  
  if(grepl('Pelangi', radius_all$videofile[b], ignore.case = TRUE) == TRUE){
    
    radius_all$ID[b] = 'Pelangi'
    
  }  else if (grepl('Fajar', radius_all$videofile[b], ignore.case = TRUE) == TRUE){
    
    radius_all$ID[b] = 'Fajar'
    
  } else if (grepl('Baju', radius_all$videofile[b], ignore.case = TRUE) == TRUE){
    
    radius_all$ID[b] = 'Baju'
    
  } else if (grepl('Roger', radius_all$videofile[b], ignore.case = TRUE) == TRUE){
    
    radius_all$ID[b] = 'Roger'
    
  } else {
    
    radius_all$ID[b] = 'NA'
  }
  
}

radius_all <- left_join(radius_all, normalization_value, by = 'videofile')

# 04b: post-processing II: smoothing ----

radius_all <- radius_all %>%
  mutate(radius = as.numeric(radius))

# the first and potentially last radius of each video is flawed using the function like this
radius_all <- radius_all %>%
  dplyr::mutate(radius_smoothed =  butter.it(radius_all$radius, samplingrate =  25, order = 4, lowpasscutoff = 10))


# 04b: post-processing II: add re-scaled radius ---- 

radius_all <- radius_all %>% 
  group_by(videofile) %>% 
  mutate(radius_scaled = scales::rescale(radius, to = c(0,1)))

# 05: saving data ----

saveRDS(radius_all, file = "DLC_estimated_radii_meta_data.rds")

