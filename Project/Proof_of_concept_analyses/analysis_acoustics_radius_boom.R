# analysis: how do acoustic features of boom call relate to its radius?

# Manuscript: A toolkit for the dynamic study of spherical biological objects
# author: Dr. Lara S. Burchardt
# start point: dataset with acoustic features and dataset with radius information
# end point: matched dataset of acoustics and radii + analysis and visualization

############

# 00: load packages ----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse","conicfit", "scales", "spiro", "signal", "foreach")
library("foreach")

# 01: data ----

## 01a: load data ----
acoustics_boom <- readRDS("acoustic_analysis_boom_multi_checked.rds")

radius_boom <- readRDS("radius_estimation_boom_proof1.rds")

## 01b: data preparation ----


### get video sampling rate from video name (in June recordings were made with 25 fps, in August with 50)

for (a in 1: nrow(acoustics_boom)){
  
  if(grepl('June', acoustics_boom$audiofile[a], ignore.case = TRUE) == TRUE){
    
    acoustics_boom$fs_video[a] <- 25
    
  }  else if (grepl('August', acoustics_boom$audiofile[a], ignore.case = TRUE) == TRUE){
    
    acoustics_boom$fs_video[a] <- 50
    
  } 
}

radius_boom <- radius_boom %>% 
  mutate(frame = as.numeric(frame))

### for 25 fps videos, 2 audio snippets need to be combined to match
### one video frame: every pair is combined, how? just take the first one? or average? 

subset_fps25 <- acoustics_boom %>% 
  dplyr::filter(fs_video == 25) %>% 
  group_split(audiofile)

subset_fps50 <- acoustics_boom %>% 
  dplyr::filter(fs_video == 50)

averaged_fps25 <- data.frame()

for(a in 1:length(subset_fps25)){
  
  data <- subset_fps25[[a]]
  
  if((nrow(data) %% 2) == 0) {
    
    data <- data
    
  } else {
    
    data <- head(data, -1)
  }
  
  counter = 0
  matchname <- data$match[1]
  audiofilename <- data$audiofile[1]
  
for( b in seq(from = 1, to = nrow(data), by = 2)){
  
  counter = counter + 1
  data$group[b] <- counter
  data$group[b+1] <- counter
  
}
  #! check if means are actually calculated, if so move on with checking the combination with bind_rows with fps 50 parts
  # which information are we losing? 
  data_average_per_group <- data %>% 
  group_by(group) %>% 
  summarise_at(vars(amEnvDep:fundamental_mean),list(mean))
  
  data <- data_average_per_group %>% 
    mutate(audiofile = audiofilename,
           match = matchname,
           fs_video = 25,
           frame = group)
  
  averaged_fps25 <- rbind(averaged_fps25, data)
  
}


##problem: unmatching column names now
acoustics_boom_averagedfps25 <- rbind(averaged_fps25, subset_fps50)

acoustics_boom_averagedfps25<- dplyr::bind_rows(averaged_fps25, subset_fps50)

### prepare match names to match audio to radius snippets

# video match name
filename <- strsplit(radius_boom$videofile, split  = "DLC_resnet")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))

radius_boom$match <- filename[,1]

# audio match name

filename <- strsplit(acoustics_boom$audiofile, split  = "multisource")
filename <- as.data.frame(matrix(unlist(filename),ncol=2,byrow=T))
acoustics_boom$match <- filename[,2]
acoustics_boom$match <- str_sub(acoustics_boom$match, 1, -5)

# for (b in 1: nrow(acoustics_boom)){
#   
#   if(grepl('Fajar', acoustics_boom$match[b], ignore.case = TRUE) == TRUE){
#     
#     # acoustics_boom$match string -9 (fajar.wav is 9 characters)
#     acoustics_boom$match[b] <- str_sub(acoustics_boom$match[b], 1, -10)
#     
#     
#   }  else if (grepl('Baju', acoustics_boom$match[b], ignore.case = TRUE) == TRUE){
#     
#     # acoustics_boom$match string -8 (baju.wav is 8 characters)
#     acoustics_boom$match[b] <- str_sub(acoustics_boom$match[b], 1, -9)
#     
#   } else if (grepl('Pelangi', acoustics_boom$match[b], ignore.case = TRUE) == TRUE){
#     
#     # acoustics_boom$match string -11 (pelangi.wav is 11 characters)
#     acoustics_boom$match[b] <- str_sub(acoustics_boom$match[b], 1, -12)
#     
#   } else if (grepl('Roger', acoustics_boom$match[b], ignore.case = TRUE) == TRUE) {
#     
#     # acoustics_boom$match string -9 (roger.wav is 9 characters)
#     acoustics_boom$match[b] <- str_sub(acoustics_boom$match[b], 1, -10)
#     
#   }
# }

### matching radius and acoustics dataframe by match name and frame

comparison_radius_boom <- left_join(radius_boom,acoustics_boom,  by= c("match", "frame"))

# filter radii that are definitly wrong, i.e. very large
# we set the maximum to the approx. maximum we found when tracking manually

comparison_radius_boom <- comparison_radius_boom %>%
  mutate(radius = as.numeric(radius)) %>% 
  dplyr::filter(radius < 270)


# 02: analysis ----
# comparison boom acoustics and radius 

# 02a: correlation matrix ----

comparison_radius_boom_numeric <- comparison_radius_boom %>% 
  select(-audiofile, -match, -videofile, -voiced, -fundamental_mean, -fs_video, -frame) 


cor_all <- cor(comparison_radius_boom_numeric)
cor_all_2 <-rcorr(as.matrix(comparison_radius_boom_numeric))

test <- cor_all_2$r


#correlation plot, show all
corrplot(cor_all_2$r, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

#correlation plot, only show significant correlations
corrplot(cor_all_2$r, type="lower", order="original", 
         p.mat = cor_all_2$P, sig.level = 0.05, insig = "blank", diag = FALSE)
