#### Logistic regression testing the probability of a parrotfish bite based on MBB density ###
#### Script by Nyssa Silbiger #####
### Edited on 11/16 2018 #########
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
library(car)

## load data ###############
### This is the time series from Mallory Rice on bites and bore holes per coral 
#TSData<-read.csv('Data/TimeSeries_Photo_Data_7.3.2018.csv')
TSData<-read.csv('Data/TimeSeries_Photo_Data_10.25.2018.csv')

## read in coral cover data
coverdata<-read.csv('Data/knb-lter-mcr.4_1_20151209.csv')
##read in fish data
FData<-read.csv('Data/MCR_LTER_Annual_Fish_Survey_20160509.csv')
# read in the CHN Data
NutData<-read.csv('Data/Macroalgal CHN.csv')
#read in groundtruthing data
TruthData<-read.csv('Data/groundtruth.csv')


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
TSData$bore.m2.scaled<-as.numeric(scale(x = TSData$bore.m2, scale = TRUE))
mod2<-glmer(bite~bore.m2.scaled+(1|Year:Site), data=TSData, family = 'binomial')
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
plot(TSData$bore.m2, TSData$bite, xlab = expression(paste('Density of '*italic(Lithophaga)*' (counts per m'^2,')')),
     ylab = 'Probability of parrotfish scar', cex.lab=1.5, cex.axis=1.5, col = 'grey', pch = 19, cex = 0.5)
#curve(predict(mod2,data.frame(bore.cm2=x),type="resp",re.form=NA),add=TRUE, col = 'black', lwd=2) # draws a curve based on prediction from logistic regression model
#lines(newdat$x, y,col = 'black', lwd=2)
lines(SE$data$bore.m2, SE$data$bite, lwd=2) # prediction
lines(SE$data$bore.m2, SE$data$ymin) # SE lines
lines(SE$data$bore.m2, SE$data$ymax) # SE lines
# fill in with grey polygon
polygon(c(SE$data$bore.m2,rev(SE$data$bore.m2)),c(SE$data$ymax,rev(SE$data$ymin)),col=grey2, border = NA)
dev.off()

## assumptions of homoscedasticity
boxplot(resid(mod3)~TSData$Site)
boxplot(resid(mod3)~TSData$Year)

## Raw scaridae data
Scaridae<- FData %>%
  filter(Site == 1| Site == 3 | Site==4) %>%
  filter(Family=='Scaridae' & Habitat=='FR') # only put scarides on fringing reef (might want to only include greater than 150mm)

Scaridae$Site<-as.factor(Scaridae$Site)
levels(Scaridae$Site)<-c('LTER1','LTER3','LTER4') # make the levels the same in this dataset



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
  summarise(bore = mean(bore.m2),  bore.se = sd(bore.m2)/sqrt(n()))
#rename the levels
levels(bore$Site)<-c("LTER 1","LTER 2","LTER 3","LTER 4","LTER 5","LTER 6")

# bring together the data frames
Bore.algae<-left_join(bore, Turb)

# bring in the rapid data points for 2016 from Tom (These were not available from the long-term LTER data set). 
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
par(mar=c(5.1,8.3,4.1,2.1))
# make a plot of tissue N versus density of borers
plot(Bore.algae$N, Bore.algae$bore, cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,ylim = c(0,1100),
     pch = 19, xlab = '% Tissue N', ylab =NA)
mtext(text=expression(paste('Mean density of', italic(' Lithophaga'))), side = 2, line = 4, cex = 1.5 )
mtext(text=expression(paste(' (counts per m'^{2},')')), side = 2, line = 2.2, cex = 1.5 )

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

## TWo way ANOVAS for lithophaga densities, parrotfish biomass, and percent cover of massive porites by site and year
#Lithophaga
MBB.mod<-lm(bore.m2^(1/4)~Site*Year , data = TSData) # need to 4th root transform to meet assumptions
anova(MBB.mod)
summary(MBB.mod)

#TSData$bore.m21<-(TSData$bore.m2+1)
#MBB.mod<-glm(bore.m21~Site*Year, family  = Gamma(link = 'log'), data = TSData) # need to 4th root transform to meet assumptions
qqnorm(resid(MBB.mod))
qqline(resid(MBB.mod))
hist(resid(MBB.mod))
# Not good... very zero-inflated and not normal

#scars
Scar.mod<-lm(log(bites.m2+1)~Site*Year , data = TSData) 
anova(Scar.mod)
summary(Scar.mod)

qqnorm(resid(Scar.mod))
qqline(resid(Scar.mod))
hist(resid(Scar.mod))
# Not good... very zero-inflated and not normal

# run a hurdle model (Zero-inflated Gamma) because too many zeros and non-normal
TSData$bore<-ifelse(TSData$bore.cm2>0,1,0) # make a column of 1's and 0's 

MBB.mod1 <- glm(bore ~ Site*Year, data = TSData, family = binomial(link = logit)) # the logisitic part (0/1)
MBB.mod2 <- glm(bore.m2 ~ Site*Year, data = subset(TSData, bore == 1), family = Gamma(link = log)) # The gamma part (>0)

#qqnorm plots
qqnorm(resid(MBB.mod2))
qqline(resid(MBB.mod2))
hist(resid(MBB.mod2))
#much better!

#binomial coefficient
bin_coef <- plogis(coef(MBB.mod1)[[1]])
(plogis(confint(MBB.mod1)))
Anova(MBB.mod1, type = 3)
summary(MBB.mod1)

#gamma coefficient
gamma_coef <- exp(coef(MBB.mod2)[[1]])
(exp(confint(MBB.mod2)))
Anova(MBB.mod2, type = 3)
summary(MBB.mod2)

#scars
Scar.mod1 <- glm(bite ~ Site*Year, data = TSData, family = binomial(link = logit)) # the logisitic part (0/1)
Scar.mod2 <- glm(bites.m2 ~ Site*Year, data = subset(TSData, bite == 1), family = Gamma(link = log)) # The gamma part (>0)
#check the residuals
qqnorm(resid(Scar.mod2))
qqline(resid(Scar.mod2))
hist(resid(Scar.mod2))
#binomial
Anova(Scar.mod1)
summary(Scar.mod1)
#gamma
Anova(Scar.mod2)
summary(Scar.mod2)

# Fish counts
Scaridae.mod <-lm(log(Count+1)~Site*Year, data = Scaridae)
anova(Scaridae.mod)
summary(Scaridae.mod)

qqnorm(resid(Scaridae.mod))
qqline(resid(Scaridae.mod))
hist(resid(Scaridae.mod))

# Porites cover
#Convert Porites SA to percent cover
TSData$PoritesCover<-100*(TSData$Surface.area.m2/0.25)  # quadrats were 0.25 m2
Porites.mod<-lm(PoritesCover^(1/4) ~Site*Year , data = TSData)
anova(Porites.mod)
summary(Porites.mod)

qqnorm(resid(Porites.mod))
qqline(resid(Porites.mod))
hist(resid(Porites.mod))


### Ground truth the 2D to 3D counts
# simple linear regression between the 2D and 3D counts
# scars
mod.truth.scars<-lm(log(dcount_scar+1)~log(pcount_scar+1), data = TruthData)
hist(resid(mod.truth.scars))  
qqnorm(resid(mod.truth.scars))
qqline(resid(mod.truth.scars))
anova(mod.truth.scars)  
summary(mod.truth.scars)
# p<0.001, R2 = 0.58, slope = 0.59, 2D slightly over estimated (intercept = 1.69)

#Lithophaga
mod.truth.bore<-lm(log(dcount_mbb+1)~log(pcount_mbb+1), data = TruthData)
hist(resid(mod.truth.bore))  
qqnorm(resid(mod.truth.bore))
qqline(resid(mod.truth.bore))
anova(mod.truth.bore)  
summary(mod.truth.bore)
# p<0.001, R2 = 0.58, slope = 1.17, 2D slightly under estimated (intercept -0.19)

par(mfrow=c(1,2))
#plot the parrotfish scars
plot(log(TruthData$pcount_scar+1), log(TruthData$dcount_scar+1), pch = 19, xlab = 'log(Scar photo counts +1)', ylab = 'log(Scar diver counts +1)')
abline(1,1, lty = 2)
pred<-predict(mod.truth.scars, se.fit = TRUE)
ind<-order(TruthData$pcount_scar)
x<-log(TruthData$pcount_scar[ind]+1)
lines(x, pred$fit[ind])
lines(x, pred$fit[ind]+pred$se.fit[ind])
lines(x, pred$fit[ind]-pred$se.fit[ind])
polygon(c(x,rev(x)),c(pred$fit[ind]+pred$se.fit[ind],rev(pred$fit[ind]-pred$se.fit[ind])),col=grey2, border = NA)

# plot the lithophagids
plot(log(TruthData$pcount_mbb+1), log(TruthData$dcount_mbb+1), pch = 19, xlab = 'log(Lithophaga photo counts +1)', ylab = 'log(Lithophaga diver counts +1)')
abline(1,1, lty = 2)
pred<-predict(mod.truth.bore, se.fit = TRUE)
ind<-order(TruthData$pcount_mbb)
x<-log(TruthData$pcount_mbb[ind]+1)
lines(x, pred$fit[ind])
lines(x, pred$fit[ind]+pred$se.fit[ind])
lines(x, pred$fit[ind]-pred$se.fit[ind])
polygon(c(x,rev(x)),c(pred$fit[ind]+pred$se.fit[ind],rev(pred$fit[ind]-pred$se.fit[ind])),col=grey2, border = NA)

