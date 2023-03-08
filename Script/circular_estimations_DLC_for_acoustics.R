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

## starting main loop ----
for (a in 1:length(list_of_files)) {
  
  radius <- data.frame()
  
  auto_data <- read_delim(paste(path, list_of_files[a], sep = "\\"), delim = "," )

# 02: pre-processing data ----

  data_circle_estimation <- data_prep_radius_estim_DLC(data = auto_data)

# 03: circle estimation ----

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


# 04a: post-processing I: adding meta-data, like sex and ID ----

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

saveRDS(comparison_data_scaled, file = "DLC_estimated_radii_meta_data.rds")
  
  