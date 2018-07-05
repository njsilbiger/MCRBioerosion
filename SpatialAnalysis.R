
### Spatial analysis of MBB and parrotfish bites within a coral #####
### Created by Nyssa Silbiger ##########
### Edited on 7/5/2018 ##########

######## Clear the workspace ###########
rm(list=ls())

##### load the libraries ###########
library(tidyverse)

######## load the data ###############
XYData<-read.csv('Data/XYCoordData.csv')

#### quick plot to visualize MBB and parrotfish bites on one coral through time ####
### pull ou the coral outline
subData<-XYData %>%
  filter(ID=='coral 1A' & Item == 'coral outline')
### pull out the MBB and bite coords
subData.points<-XYData %>%
  filter(ID=='coral 1A' & Item != 'coral outline')

year<-unique(subData$Year)
# plot images of individual coral through time
par(mfrow=c(3,2))
for(i in 1:length(unique(subData$Year))){
plot(subData$x[subData$Year==year[i]], subData$y[subData$Year==year[i]], type = 'l', xlim = c(min(subData$x), max(subData$x)), 
     ylim = c(min(subData$y), max(subData$y)), main = year[i], xlab = "", ylab = "", xaxt='n', yaxt= 'n')
points(subData.points$x[subData.points$Year==year[i]& subData.points$Item=='bite'], subData.points$y[subData.points$Year==year[i]& subData.points$Item=='bite'], col = 'red', pch = 6)
points(subData.points$x[subData.points$Year==year[i]& subData.points$Item=='macroborer'], subData.points$y[subData.points$Year==year[i]& subData.points$Item=='macroborer'], col = 'blue', pch = 19)

}
legend('topleft', legend = c('Coral outline','Macroborer',' Fish bite'), col = c('black','blue','red'), lty=c(1,0,0), pch = c(0,19,6), bty = 'n')


