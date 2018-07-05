#### Logistic regression testing the probability of a parrotfish bite based on MBB density ###
#### Script by Nyssa Silbiger #####
### Edited on 7/5 2018 #########
###############################


# clear the workspace
rm(list=ls())

## load libraries ##
library(tidyverse)
library(lme4)
library(lmerTest)

## load data ###############
### This is the time series from Mallory Rice on bites and bore holes per coral 
TSData<-read.csv('Data/TimeSeries_Photo_Data_7.3.2018.csv')
## read in coral cover data
coverdata<-read.csv('Data/knb-lter-mcr.4_1_20151209.csv')
##read in fish data
FData<-read.csv('Data/MCR_LTER_Annual_Fish_Survey_20160509.csv')


### analysis #############
# normalize bore hole and bite data per cm2
TSData$bore.cm2<-TSData$No.Macroborers/TSData$Surface.area.cm2
TSData$bites.cm2<-TSData$No.Bites/TSData$Surface.area.cm2

# label which corals have a macroborer
yes<-which(TSData$bore.cm2>0)
# calculate a histogram of bites on a coral with and without a macroborer 
a<-hist(TSData$bites.cm2[yes])
b<-hist(TSData$bites.cm2[-yes])
# plot the histogram
plot(a$mids, a$density, 'b')
points(b$mids, b$density, 'b', col = 'red')



# calculate mean porites cover
porites.means<- coverdata %>%
  group_by(Taxonomy...Substrate...Functional.Group, Site, Habitat,Date)%>%
  summarise(.,means = mean(Percent.Cover))%>% 
  filter(Taxonomy...Substrate...Functional.Group=='Porites spp. Massive' & Habitat=='Fringing')%>%
  separate(col = Date,into = 'Year',sep =  "-", extra = 'drop')

# add a 0 or 1 for corals that do and don't have bites
TSData$bite<-ifelse(TSData$bites.cm2>0,1,0)

# run a logisitic regression for probability of getting bitten with density of macroborers
mod1<-glm(bite~bore.cm2, data=TSData, family = 'binomial')
mod2<-glmer(bite~bore.cm2+(1|Year:Site), data=TSData, family = 'binomial')

# plot the best fit line (need to update to plot the effect size from MEM)
png('Output/Logistic plot.png', width = 500, height = 500)
plot(TSData$bore.cm2, TSData$bite, xlab = expression(paste('Density of borers (counts/cm'^2,')')),
     ylab = 'Probability of coral being bitten', cex.lab=1.5, cex.axis=1.5)
curve(predict(mod1,data.frame(bore.cm2=x),type="resp"),add=TRUE, col = 'red', lwd=2) # draws a curve based on prediction from logistic regression model
points(TSData$bore.cm2,fitted(mod1),pch=20, col = 'red') # optional: you could skip this draws an invisible set of points of body size 
dev.off()

# total number of scarides per year and site
Scaridae.counts<- FData %>%
  filter(Family=='Scaridae' & Habitat=='FR')%>% # only put scarides on fringing reef (might want to only include greater than 150mm)
    group_by(Year,Site)%>%
  summarise(.,sum.fish = sum(Count)) # total number of fish

Scaridae.counts$Site<-as.factor(Scaridae.counts$Site)
levels(Scaridae.counts$Site)<-c('LTER1','LTER2','LTER3','LTER4','LTER5','LTER6') # make the levels the same in this dataset

# total number of bites per year and site
TotalBites<-TSData%>%
group_by(Year,Site)%>%
  summarise(.,mean.scars = mean(bites.cm2), mean.bore = mean(bore.cm2) )
TotalBites$Year<-as.factor(TotalBites$Year)

FishBites<-left_join(TotalBites, Scaridae.counts)

plot(FishBites$sum.fish, FishBites$mean.scars)
