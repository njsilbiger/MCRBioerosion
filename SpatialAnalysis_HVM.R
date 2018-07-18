
### Spatial analysis of MBB and parrotfish bites within a coral #####
### Created by Nyssa Silbiger ##########
### Edited on 7/5/2018 ##########

######## Clear the workspace ###########
rm(list=ls())

##### load the libraries ###########
library(tidyverse)

##HOLLYADDED
library(SDMTools)
library(sp)
##ENDHOLLYADDED

######## load the data ###############

##HOLLYADDED
setwd('~/Box Sync/MCRBioerosion-master')
##ENDHOLLYADDED

XYData<-read.csv('Data/XYCoordData.csv')

#### quick plot to visualize MBB and parrotfish bites on one coral through time ####
### pull ou the coral outline
levels(XYData$ID)
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



##HOLLYADDED

LoopOverCorals <- function(){}
# Generate a dataframe to store t-stats and p-values
XYData$ID.Year <- as.factor(paste(XYData$ID,XYData$Year,sep='.'))
length(levels(XYData$ID.Year))
boots <- 1000 # Number of bootstrapping runs (generations of null bite distributions, followed by t-test)
full.t <- as.data.frame(matrix(c(rep(NA,boots*length(levels(XYData$ID.Year)))),nrow=boots,ncol=length(levels(XYData$ID.Year))))
colnames(full.t)<-levels(XYData$ID.Year)
full.p <- full.t
full.obs.dist <- full.t

SummaryStats <- function(){}
summstat <- as.data.frame(cbind(levels(XYData$ID.Year)))
summstat$mean.p <- NA
summstat$mean.t <- NA
summstat$mean.obs.dist <- NA
summstat$mean.null.dist <- NA
summstat$pval.ttest.p <- NA
summstat$pval.ttest.t <- NA
summstat$pval.ttest.lc <- NA
summstat$pval.ttest.uc <- NA



for(corctr in 1:dim(full.t)[2]){
	corchoice <- colnames(full.t)[corctr]


ExtractDataByCoral.Year <- function(){}

#corchoice <- 'coral 1A.2008'
# Find outline of the coral, which will be input vertices of the polygon
outline <- XYData[XYData$ID.Year==corchoice & XYData$Item=='coral outline',]

# ID bites and borers
bites <- XYData[XYData$ID.Year==corchoice & XYData$Item=='bite',]
borers <- XYData[XYData$ID.Year==corchoice & XYData$Item=='macroborer',]

# Proceed only if (1) there's a coral outline, (2) there are bites
if(min(c(dim(outline)[1],dim(bites)[1]))>1){
	# and there is at least one bore hole
	if(dim(borers)[1]>0){
	
# Check to be sure all the bite/borer data are within this outline (a "1" means "TRUE")
point.in.polygon(borers$x,borers$y,outline$x,outline$y)
point.in.polygon(bites$x,bites$y,outline$x,outline$y)

MeasureBiteDistance <- function(){}
# For each bite, measure distance to nearest macroborer location:
# Create a matrix of observed distances between each bite and borer
obs.dist <- as.data.frame(c(1:dim(bites)[1]))
for(borecount in 1:dim(borers)[1]){
	obs.dist <- cbind(obs.dist,spDistsN1(cbind(bites$x,bites$y),c(borers$x[borecount],borers$y[borecount])))
}
colnames(obs.dist)[2:(1+dim(borers)[1])]<-c(paste('borer',1:dim(borers)[1],sep=' '))

# Choose the minimum distance
obs.dist$min.dist <- NA
for(bitecount in 1:dim(bites)[1]){
	obs.dist$min.dist[bitecount] <- min(obs.dist[bitecount,2:(dim(borers)[1]+1)])
}
hist(obs.dist$min.dist)

# Rescale that distance with key?



Bootstrap.NullDistrib.T.Test <- function(){}
boots <- boots;
tstats <- rep(NA,boots)
pvals <- tstats
for(bootct in 1:boots){
	
GenerateNullDistrib <- function(){}
# Randomly place the same number of points on the coral head as bites
bitenum <- dim(bites)[1]
null.bites <- as.data.frame(cbind(c(1:bitenum),rep(NA,bitenum),rep(NA,bitenum)))
colnames(null.bites)<-c('NullBite','x','y')

bitecount <- 1
while(bitecount < bitenum+1){
	null.bites$x[bitecount]<- runif(1,min=min(outline$x),max=max(outline$x))
	null.bites$y[bitecount]<- runif(1,min=min(outline$y),max=max(outline$y))
	if(point.in.polygon(null.bites$x[bitecount],null.bites$y[bitecount],outline$x,outline$y)==1){
		bitecount <- bitecount+1
	}
}

# Check to be sure all these random bites are within the outline
# point.in.polygon(null.bites$x,null.bites$y,outline$x,outline$y)
# par(mar=c(4,4,1,1),mfrow=c(1,1))
# plot(outline$x,outline$y)
# points(null.bites$x,null.bites$y,col='red')
# points(bites$x,bites$y,col='blue')
# points(borers$x,borers$y,col='green')

# Calculate distance of these bites to nearest macroborer location
null.dist <- as.data.frame(c(1:dim(null.bites)[1]))
for(borecount in 1:dim(borers)[1]){
	null.dist <- cbind(null.dist,spDistsN1(cbind(null.bites$x,null.bites$y),c(borers$x[borecount],borers$y[borecount])))
}
colnames(null.dist)[2:(1+dim(borers)[1])]<-c(paste('borer',1:dim(borers)[1],sep=' '))

# Choose the minimum distance
null.dist$min.dist <- NA
for(bitecount in 1:dim(bites)[1]){
	null.dist$min.dist[bitecount] <- min(null.dist[bitecount,2:(1+dim(borers)[1])])
}

# # quartz()
# par(mar=c(4,4,1,1),mfrow=c(2,1))
# hist(obs.dist$min.dist)
# hist(null.dist$min.dist)  ## NEED TO FIGURE OUT HOW TO PLOT THESE OVERLAPPING

# Perform t-test on distributions
tstats[bootct]<-t.test(null.dist$min.dist,obs.dist$min.dist)$statistic
pvals[bootct]<-t.test(null.dist$min.dist,obs.dist$min.dist)$p.value

full.obs.dist[bootct,corctr] <- mean(null.dist$min.dist)

} # End Bootstrapping Loop

# hist(tstats)
# hist(pvals)

# Write Data
full.t[,corctr]<- tstats
full.p[,corctr]<- pvals
summstat$mean.p[corctr] <- mean(pvals)
summstat$mean.t[corctr] <- mean(tstats)
summstat$mean.obs.dist[corctr] <- mean(obs.dist$min.dist)
summstat$mean.null.dist[corctr] <- mean(full.obs.dist[,corctr])

ttesthold <- t.test(pvals,mu=.05,alternative='less')
summstat$pval.ttest.p[corctr] <- ttesthold$p.value
summstat$pval.ttest.t[corctr] <- ttesthold$statistic
summstat$pval.ttest.lc[corctr] <- ttesthold$conf.int[1]
summstat$pval.ttest.uc[corctr] <- ttesthold$conf.int[2]

} # End If/Then test of whether there is at least one borer hole
} # End If/Then test of whether coral exists & has bites
} # End looping through coral ID/Year combos

hist(full.p,na.rm=TRUE)
head(full.p)
head(full.t)

SaveFiles <- function(){}
#write.csv(summstat,'Data/SummarySpatialStatistics.csv')

hist(summstat$pval.ttest.p)
tail(summstat)
summstat[!is.na(summstat$mean.p),]
dim(summstat[!is.na(summstat$mean.p),])
##ENDHOLLYADDED





# HOLLY TROUBLESHOOTING MISCELLANY
corctr
par(mar=c(4,4,1,1),mfrow=c(1,1))
plot(outline$x,outline$y)
points(null.bites$x,null.bites$y,col='red')
points(bites$x,bites$y,col='blue')
points(borers$x,borers$y,col='green')


