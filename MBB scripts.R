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
library(boot)
library(grDevices)
library(jtools)

## load data ###############
### This is the time series from Mallory Rice on bites and bore holes per coral 
TSData<-read.csv('Data/TimeSeries_Photo_Data_7.3.2018.csv')
## read in coral cover data
coverdata<-read.csv('Data/knb-lter-mcr.4_1_20151209.csv')
##read in fish data
FData<-read.csv('Data/MCR_LTER_Annual_Fish_Survey_20160509.csv')
# read in the CHN Data
NutData<-read.csv('Data/Macroalgal CHN.csv')


### analysis #############
# Filter the TS data to include only LTER1, 3, and 4 since the other sites have very low sample sizes of corals
TSData<- TSData %>%
  filter(Site == 'LTER1'| Site == 'LTER3' | Site=="LTER4")

# convert surface area from cm2 to m2
TSData$Surface.area.m2<-TSData$Surface.area.cm2*1e-4

# normalize bore hole and bite data per cm2
TSData$bore.cm2<-TSData$No.Macroborers/TSData$Surface.area.cm2
TSData$bites.cm2<-TSData$No.Bites/TSData$Surface.area.cm2

#normalized to per m2
TSData$bore.m2<-TSData$No.Macroborers/TSData$Surface.area.m2
TSData$bites.m2<-TSData$No.Bites/TSData$Surface.area.m2

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
  filter(Site == 'LTER 1'| Site == 'LTER 3' | Site=="LTER 4" & Habitat =='Fringing') %>%
  filter(Taxonomy...Substrate...Functional.Group=='Porites spp. Massive' | Taxonomy...Substrate...Functional.Group=='Porites')%>%
  separate(col = Date,into = 'Year',sep =  "-", extra = 'drop') %>%
  group_by(Taxonomy...Substrate...Functional.Group, Site, Habitat,Year)%>%
  summarise(.,means = mean(Percent.Cover)) 
  
  
# add a 0 or 1 for corals that do and don't have bites
TSData$bite<-ifelse(TSData$bites.cm2>0,1,0)

# run a logisitic regression for probability of getting bitten with density of macroborers
#mod1<-glm(bite~bore.m2, data=TSData, family = 'binomial')
# scaled values and mixed effect model to account for site by year variance
mod2<-glmer(bite~scale(x = bore.m2, scale = TRUE)+(1|Year:Site), data=TSData, family = 'binomial')
summary(mod2)

#get robust standard errors and unscaled data for the plot
mod3<-glmer(bite~ bore.m2+(1|Year:Site), data=TSData, family = 'binomial')
SE<-effect_plot(mod3, pred = bore.m2, robust = TRUE, interval = TRUE)

#calculate the odds ratio
# back transform the scaled coefs
beta1<-fixef(mod2)[2]/attributes(scale(TSData$bore.m2,center=F))$"scaled:scale"
beta0<-fixef(mod2)[1]/attributes(scale(TSData$bore.m2,center=F))$"scaled:scale"

odds<-exp((beta0 +beta1*100) - (beta0 +beta1*99))
## for every increase in borer density of 1/m2 there is a .00163% increase in the odds of being bitten
# convert to per cm2 for easier interpretation
odds.cm2<-10000*(odds-1)
## for every increase in borer density of 1/cm2 there is a 163% increase in the odds of being bitten 

## Plot the logistic regression
newdat<-data.frame(x=seq(0,17000,length=20))
#mygrey color
grey2<-adjustcolor( "grey", alpha.f = 0.2)

pdf('Output/Logistic plot2.pdf', 6,6,useDingbats = FALSE)
par(mar=c(5.1,5.3,4.1,2.1))
plot(TSData$bore.m2, TSData$bite, xlab = expression(paste('Density of borers (counts/m'^2,')')),
     ylab = 'Probability of coral scar', cex.lab=2, cex.axis=1.5, col = 'grey', pch = 19, cex = 0.5)
#curve(predict(mod2,data.frame(bore.cm2=x),type="resp",re.form=NA),add=TRUE, col = 'black', lwd=2) # draws a curve based on prediction from logistic regression model
#lines(newdat$x, y,col = 'black', lwd=2)
lines(SE$data$bore.m2, SE$data$bite, lwd=2) # prediction
lines(SE$data$bore.m2, SE$data$ymin) # SE lines
lines(SE$data$bore.m2, SE$data$ymax) # SE lines
# fill in with grey polygon
polygon(c(SE$data$bore.m2,rev(SE$data$bore.m2)),c(SE$data$ymax,rev(SE$data$ymin)),col=grey2, border = NA)
dev.off()


# total number of scarides per year and site
Scaridae.counts<- FData %>%
  filter(Site == 1| Site == 3 | Site==4) %>%
  filter(Family=='Scaridae' & Habitat=='FR')%>% # only put scarides on fringing reef (might want to only include greater than 150mm)
    group_by(Year,Site)%>%
  summarise(.,sum.fish = sum(Count)) # total number of fish

Scaridae.counts$Site<-as.factor(Scaridae.counts$Site)
levels(Scaridae.counts$Site)<-c('LTER1','LTER3','LTER4') # make the levels the same in this dataset

# total number of bites per year and site
TotalBites<-TSData%>%
group_by(Year,Site)%>%
  summarise(.,mean.scars = mean(bites.cm2), mean.bore = mean(bore.cm2) )
TotalBites$Year<-as.factor(TotalBites$Year)
# join with the fish data
FishBites<-left_join(TotalBites, Scaridae.counts)
# plot the relationship
plot(FishBites$sum.fish, FishBites$mean.scars)


#### Analysis of bioeroder density and nutrients ####

## make a high vs low nutrient plot

# TSData$Nut<-factor(ifelse(TSData$Site=='LTER1'|TSData$Site=='LTER2'|TSData$Site=='LTER3', 'High','Low'))
# # calculate means by nutrients
# BitesNuts<-TSData%>%
#   group_by(Nut)%>%
#   summarise(.,mean.scars = mean(bites.cm2), mean.bore = mean(bore.cm2), se.scars = sd(bites.cm2)/sqrt(n()), se.bore = sd(bore.cm2)/sqrt(n()) )
# 
# # make a barplot with error bars
# b<-barplot(BitesNuts$mean.bore, names.arg = c('High Nutrient Sites','Low Nutrient Sites'), ylim=c(0,0.06), ylab = expression('mean macroborers cm'^-2), col = 'darkgrey')
# segments(b,BitesNuts$mean.bore+BitesNuts$se.bore, b,BitesNuts$mean.bore-BitesNuts$se.bore, col = 'black')
# #model
# modNut.bore<-lmer(bore.cm2~ Nut +(1|Year), data = TSData)
# anova(modNut.bore)

# # do the same for the scarids
# Scaridae.counts$Nut<-factor(ifelse(Scaridae.counts$Site=='LTER1'|Scaridae.counts$Site=='LTER2'|Scaridae.counts$Site=='LTER3', 'High','Low'))
# 
# ScaridaeNuts<-Scaridae.counts%>%
#   group_by(Nut)%>%
#   summarise(.,mean.scars = mean(sum.fish), se.scars = sd(sum.fish)/sqrt(n()) )
# 
# c<-barplot(ScaridaeNuts$mean.scars, names.arg = c('High Nutrient Sites','Low Nutrient Sites'), ylim=c(0,180),ylab = expression('mean total fish counts'^-2), col = 'darkgrey')
# segments(c,ScaridaeNuts$mean.scars+ScaridaeNuts$se.scars, c,ScaridaeNuts$mean.scars-ScaridaeNuts$se.scars, col = 'black')


#Filter out the 3 sites we are using and the turbanaria %N data
Turb<-NutData %>%
  filter(Habitat =='Fringe', Genus =='Turbinaria', 
         Site == 'LTER 1'| Site == 'LTER 3' | Site=="LTER 4") %>%
         group_by(Site, Year) %>%
         summarise(N = mean(N, na.rm=T)) 

#Calculate the mean bioroder densities per site and year
bore<-TSData %>%
  filter(Site == 'LTER1'| Site == 'LTER3' | Site=="LTER4")%>%
  group_by(Site, Year) %>%
  summarise(bore = mean(bore.cm2),  bore.se = sd(bore.cm2)/sqrt(n()))
#rename the levels
levels(bore$Site)<-c("LTER 1","LTER 2","LTER 3","LTER 4","LTER 5","LTER 6")

# bring together the data frames
Bore.algae<-left_join(bore, Turb)

# bring in the rapid data points for 2016 (These were not available from the LTER data set)
Bore.algae$N[Bore.algae$Site=='LTER 1' & Bore.algae$Year==2016] = 0.61
Bore.algae$N[Bore.algae$Site=='LTER 4' & Bore.algae$Year==2016] = 0.59

# remove the missing values for the analysis
Bore.algae<-Bore.algae[complete.cases(Bore.algae),]
# run a model with a polynomial
N.mod<-lm(bore~N, Bore.algae)
#results
anova(N.mod)
summary(N.mod)

pdf(file = 'Output/BorervsN.pdf', width = 6, height = 6, useDingbats = FALSE)
par(mar=c(5.1,5.3,4.1,2.1))
# make a plot of tissue N versus density of borers
plot(Bore.algae$N, Bore.algae$bore, cex.lab = 2, cex.axis = 1.5, cex = 1.5,ylim = c(0,1100),
     pch = 19, xlab = '% Tissue N', ylab =expression(paste('Mean density of borers (counts/m'^2,')')) )
#generate new data for the fit
xx <- seq(0,1, length=50)
pred<-predict(N.mod, data.frame(N=xx), se.fit = TRUE)
# plot the predictions
polygon(c(xx,rev(xx)),c(pred$fit+pred$se.fit,rev(pred$fit-pred$se.fit)),col=grey2)
lines(xx, predict(N.mod, data.frame(N=xx)), lwd = 2)
lines(xx, pred$fit+pred$se.fit, lty=2)
lines(xx, pred$fit-pred$se.fit, lty=2)
segments(Bore.algae$N, Bore.algae$bore+Bore.algae$bore.se,
         Bore.algae$N, Bore.algae$bore-Bore.algae$bore.se)
dev.off()

# check for normality of residuals
qqnorm(resid(N.mod))
qqline(resid(N.mod))
# looks good!
