library(rstudioapi)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(conicfit)

parentfolder <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))  #what is the current folder?
audiof <- paste0(parentfolder, '/audio/snippets/_Opp_August_20_Session_1_zoom_syncedboom_0_1.wav')

#load in tracking
tracking  <- read.csv(paste0(parentfolder, '/snippets/tracked/_Opp_August_20_Session_1_zoom_syncedboom_0_1DLC_resnet101_Deep_AirSacTrackingV1Jan1shuffle1_500000_labeled'))
tracking <- tracking[-c(1,2),]
tracking <- as.data.frame(apply(tracking, 2, FUN = function(x)as.numeric(x)))

colnames(tracking) <- c('frames' ,
                        'upperlip_x', 'upperlip_y', 'upperlip_likelihood', 
                        'lowerlip_x', 'lowerlip_y', 'lowerlip_likelihood',
                        'nose_x', 'nose_y', 'nose_likelihoods',
                        'eyebridge_x', 'eyebridge_y', 'eyebridge_likelihood',
                        'startoutline_left_x', 'startoutline_left_y', 'startoutline_left_likelihood',
                        'startoutline_right_x', 'startoutline_right_y', 'startoutline_left_likelihood',
                        'lowestpoint_x', 'lowestpoint_y', 'lowestpoint_x_likelihood',
                        'midlowleft_x', 'midlowleft_y', 'midlowleft_likelihood',
                        'midlowright_x', 'midlowright_y', 'midlowright_likelihood')

totransformx <- tracking[, c(14,15, 20,23, 24)]
totransformy <- tracking[, c(15,18, 21,24, 27)]
confidences  <- tracking[, c(16,19, 22,25, 28)]
#calculate some articulatory kinematics
#distance lip movements


tracking$radius <- NA
for(i in 1:nrow(totransform)) #loop through frames and estimate circles
{
  XY <- cbind(as.numeric(as.vector(totransformx[i,])), as.numeric(as.vector(totransformy[i,])))
  #only keep points which have .60 confidence
  confidence <- confidences[i,]
  XY <- XY[confidence>.6,]
  
  if(sum(confidence>.6)> 2)
  {
  cic <- LMcircleFit(XY, LambdaIni = 1, epsilon = 1e-06, IterMAX = 50)
  tracking$radius[i] <- cic[,3]
  }
}

sub <- cbind.data.frame(tracking$frames, tracking$radius)
colnames(sub) <- c('frames', 'radius')
ggplot(sub, aes(x=frames)) + geom_point(aes(y=radius))+theme_cowplot()


