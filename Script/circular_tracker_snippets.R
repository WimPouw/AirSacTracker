# circular tracker
# snippets

curfol <- getwd() #set to source location by going to session -> set working directory
basefol <- dirname(curfol)

dataset_to_be_stored <- paste0(curfol, '/Manuscript/datasets_for_statistics/comparison_manual_tracked_v1.csv')

#merge mOuaz's data
modata <- list.files(paste0(curfol, '/AirSacTracker/TestFrames/manually_tracked_results/batches/'), pattern = '.csv')

mdat <- data.frame()
for(i in modata){
  mo <- read.csv(paste0(curfol, '/AirSacTracker/TestFrames/manually_tracked_results/batches/', i))
  mdat<- bind_rows(mdat, mo)
}
mdat$name <- sub("^.*:", "", mdat$Label)
mdat$name <- sub(".jpeg", "", mdat$name)

#if(!file.exists(dataset_to_be_stored)){
  man_data <- as.data.frame(mdat)
  #add transformation from perimeter to radius in manual data
  man_data$radius_man <- round(man_data$Perim./(2*pi))
  
  
  #only radii over 100 pixels are admitted (as below are non-tracked)
  man_data <- man_data %>% 
    dplyr::filter(radius_man > 100 | !is.na(radius_man))
#}

write.table(man_data, file = "manually_tracked_airsac_radii.csv", sep = ',', col.names = TRUE, row.names = FALSE)
