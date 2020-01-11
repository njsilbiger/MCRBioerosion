#### Logistic regression testing the probability of a parrotfish bite based on MBB density ###
#### Script by Nyssa Silbiger #####
### Edited on 11/13/19 #########
###############################


# clear the workspace
rm(list=ls())

## load libraries ##
library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(grDevices)
library(RColorBrewer)
library(cowplot)
library(car)
library(scales)
library(MuMIn)


## load data ###############
### This is the time series from Mallory Rice on bites and bore holes per coral 
#TSData<-read.csv('Data/TimeSeries_Photo_Data_7.3.2018.csv')
TSData<-read.csv('Data/TimeSeries_Photo_Data_Jan2019_final.csv')
##read in fish data
FData<-read.csv('Data/MCR_LTER_Annual_Fish_Survey_20180612.csv')
# read in the CHN Data
NutData<-read.csv('Data/MCR_LTER_Macroalgal_CHN_2005_to_2014_20151031b.csv')
#NutData<-read.csv('Data/MCR_LTER_Macroalgal_CHN_2005_to_2014_20151031.csv')

#read in groundtruthing data
TruthData<-read.csv('Data/groundtruth.csv')
# read in the coral cover data from coralnet analysis
CoralData<-read.csv('Data/CoralNet_benthic_cover_data.csv')


### analysis #############
# Filter the TS data to include only LTER1, 3, and 4 since the other sites have very low sample sizes of corals
TSData<- TSData %>%
  filter(Site == 'LTER1'| Site == 'LTER3' | Site=="LTER4")

# convert surface area from cm2 to m2
#TSData$Surface.area.m2<-TSData$Surface.area.cm2*1e-4 old
TSData$Surface.area.m2<-TSData$Ruler.surface.area..cm2.*1e-04

# normalize bore hole and bite data per cm2
TSData$bore.cm2<-TSData$No.Macroborers/TSData$Ruler.surface.area..cm2.
TSData$bites.cm2<-TSData$No.Bites/TSData$Ruler.surface.area..cm2.

#normalized to per m2
TSData$bore.m2<-TSData$No.Macroborers/TSData$Surface.area.m2
TSData$bites.m2<-TSData$No.Bites/TSData$Surface.area.m2

# label which corals have a macroborer
yes<-which(TSData$bore.cm2>0)

# remove the one crazy outlier that is an order of magnitude higher than everything
remove<-which(TSData$bore.m2>30000)
TSData<-TSData[-remove,]

# add a 0 or 1 for corals that do and don't have bites
TSData$bite<-ifelse(TSData$bites.cm2>0,1,0)

# run a logisitic regression for probability of getting bitten with density of macroborers
# scaled values and mixed effect model to account for site by year variance
TSData$bore.cm2.scaled<-as.numeric(scale(x = TSData$bore.cm2, scale = TRUE))
mod2<-glmer(bite~bore.cm2.scaled+(1|Year:Site), data=TSData, family = 'binomial')
summary(mod2)

#calculate the odds ratio
# back transform the scaled coefs
beta1<-fixef(mod2)[2]/attributes(scale(TSData$bore.cm2,center=F))$"scaled:scale"
beta0<-fixef(mod2)[1]/attributes(scale(TSData$bore.cm2,center=F))$"scaled:scale"

#Calculate odds of going from 1 to 2 per cm squared
odds<-exp((beta0 +beta1*2) - (beta0 +beta1*1))
# odds from 0.5 to 1
odds<-exp((beta0 +beta1*1) - (beta0 +beta1*0.5))

## for every increase in borer density of 1/cm2(or 10000/m2) there is a 4.78 increase in the odds of being bitten
# converted to per cm2 for easier interpretation

## Plot the logistic regression
# #mygrey color
 grey2<-adjustcolor( "grey", alpha.f = 0.2)

# Logistic plot following suggestions from https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
mod3<-glm(bite~ bore.cm2, data=TSData, family = 'binomial')

## some data to predict at: 1000 values over the range of leafHeight
ndata <- with(TSData, data_frame(bore.cm2 = seq(min(bore.cm2), max(bore.cm2),
                                               length = 1000)))
## add the fitted values by predicting from the model for the new data
ndata <-  add_column(ndata, fit = predict(mod3, newdata = ndata, type = 'response'))
TSData$bite_yn<-as.factor(ifelse(TSData$bite==1, 'Yes','No')) # for the rug in the plot
## plot it
plt <- ggplot(ndata, aes(x = bore.cm2, y = fit)) +
  geom_line() +
  geom_rug(aes(y = bite, colour = bite_yn), data = TSData) +
  scale_colour_discrete(name = 'Scar') +
  labs(x = expression(paste('Density of '*italic(Lithophaga)*' (counts per cm'^2,')')), 
       y = 'Probability of parrotfish scar')+
  theme(text = element_text(size=16))
plt

## grad the inverse link function
ilink <- family(mod3)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod3, ndata, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + 2*(se_link)),
                right_lwr = ilink(fit_link - 2*(se_link)))
## show
logisiticplot<-plt + geom_ribbon(data = ndata,
                  aes(ymin = right_lwr, ymax = right_upr),
                  alpha = 0.1)
  

## assumptions of homoscedasticity
boxplot(resid(mod3)~TSData$Site)
boxplot(resid(mod3)~TSData$Year)

## add analysis of density of borers versus scars when both borers and bites present
not0<-which(TSData$bore.cm2>0 & TSData$bites.cm2>0) # remove the 0s
bore.density.mod<-lmer(log(bites.cm2)~log(bore.cm2)+
                         (1|Year:Site), data = TSData[not0,])
anova(bore.density.mod)
summary(bore.density.mod)
qqnorm(resid(bore.density.mod))
qqline(resid(bore.density.mod))
# resquared
 r.squaredGLMM(bore.density.mod)

density.plot<-ggplot(TSData[not0,], aes(x = bore.cm2, y = bites.cm2))+
  geom_point()+
  scale_x_continuous(trans='log',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x))
                    )+
 #   breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  scale_y_continuous(trans='log',breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
#    breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  geom_smooth(method='lm',formula=y~x, color = 'black')+
 # coord_trans(x="log", y="log")+
 labs(x = expression(paste('Density of '*italic(Lithophaga)*' (counts per cm'^2,')')), 
  y = expression(paste('Density of parrotfish scar (counts per cm'^2,')')))+
  #theme_bw()+
  theme(text = element_text(size=16))

# plot the logistic and density plots next to eachother
borevsscar.plot<-plot_grid(logisiticplot, density.plot, labels = c("A", "B"))
ggsave(plot = borevsscar.plot, filename = 'Output/Figure3.pdf', 
       device = 'pdf', width = 10, height = 5)

## Raw scaridae data #####
Scaridae<- FData %>%
  filter(Site == 1| Site == 3 | Site==4) %>%
  filter(Habitat=='FR') %>% # only put scarides on fringing reef 
  mutate(Total_Length = Total_Length / 10) %>%
  filter(Total_Length>10)%>% # only include the fish >10cm because those are the ones that bite
  filter(Taxonomy == 'Chlorurus microrhinos' |
           Taxonomy == 'Chlorurus sordidus'|Taxonomy == 'Scarus frenatus'|Taxonomy == 'Scarus ghobban') %>% # known corallivores
  filter(Year == 2006 |Year == 2008| Year == 2010|Year == 2011|Year == 2013|Year == 2016)
Scaridae$Site<-as.factor(Scaridae$Site)
levels(Scaridae$Site)<-c('LTER1','LTER3','LTER4') # make the levels the same in this dataset

# total number of scarides per year and site
Scaridae.counts<- Scaridae %>%
  group_by(Year,Site)%>%
  summarise(.,sum.fish = sum(Count)) # total number of fish

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
  separate(col = N, into = "N", sep = "%")%>% # remove the % sign
  mutate(N = as.numeric(N))%>%
  group_by(Site, Year) %>%
  dplyr::summarise(N.mean = mean(N, na.rm=TRUE), se.N = sd(N)/sqrt(n())) %>%
  rename(N = N.mean)
        
#Calculate the mean bioroder densities per site and year
bore<-TSData %>%
  filter(Site == 'LTER1'| Site == 'LTER3' | Site=="LTER4")%>%
  group_by(Site, Year) %>%
  summarise(bore = mean(bore.cm2),  bore.se = sd(bore.cm2)/sqrt(n()))
#rename the levels
levels(bore$Site)<-c("LTER 1","LTER 2","LTER 3","LTER 4","LTER 5","LTER 6")

# bring together the data frames
Bore.algae<-left_join(bore, Turb)

# remove the missing values for the analysis
Bore.algae<-Bore.algae[complete.cases(Bore.algae),]
Bore.algae$Year<-as.factor(Bore.algae$Year)
# run a linear model 
N.mod<-lm(bore~N, Bore.algae)
#results
anova(N.mod)
summary(N.mod)

# pdf(file = 'Output/Figure2.pdf', width = 6, height = 6, useDingbats = FALSE)
# par(mar=c(5.1,8.3,4.1,2.1))
# j_brew_colors <- brewer.pal(n = 5, name = "Set2") # custom colors
# # make a plot of tissue N versus density of borers
# plot(Bore.algae$N, Bore.algae$bore, cex.lab = 1.5, cex.axis = 1.5, cex = 2,ylim = c(0,.20),
#      pch = as.numeric(Bore.algae$Site)+12, xlab = '% Tissue N', ylab =NA,col=j_brew_colors[as.factor(Bore.algae$Year)]  )
# #col =c(11,13,12,11,13,11,13,12)
# #j_brew_colors[c(1,2,3,1,2,1,2,3)]
# mtext(text=expression(paste('Mean density of', italic(' Lithophaga'))), side = 2, line = 4, cex = 1.5 )
# mtext(text=expression(paste(' (counts per cm'^{2},')')), side = 2, line = 2.2, cex = 1.5 )
# 
# #generate new data for the fit
# xx <- seq(0,2, length=50)
# pred<-predict(N.mod, data.frame(N=xx), se.fit = TRUE, interval="confidence",
#               level = 0.95)
# # plot the predictions
# #polygon(c(xx,rev(xx)),c(pred$fit+pred$se.fit,rev(pred$fit-pred$se.fit)),col=grey2)
# polygon(c(xx,rev(xx)),c(pred$fit[,2],rev(pred$fit[,3])),col=grey2)
# lines(xx, predict(N.mod, data.frame(N=xx)), lwd = 2)
# lines(xx, pred$fit[,2], lty=2)
# lines(xx, pred$fit[,3], lty=2)
# 
# #lines(xx, pred$fit+pred$se.fit, lty=2)
# #lines(xx, pred$fit-pred$se.fit, lty=2)
# # y error
# segments(Bore.algae$N, Bore.algae$bore+Bore.algae$bore.se,
#          Bore.algae$N, Bore.algae$bore-Bore.algae$bore.se)
# # x error
# segments(Bore.algae$N+Bore.algae$se.N, Bore.algae$bore,
#         Bore.algae$N-Bore.algae$se.N, Bore.algae$bore)
# 
# legend('topleft',c('Sites','LTER 1', 'LTER 3', 'LTER 4'), 
#        pch = c(19,13,15,16), bty='n', col = c('white','black','black','black'),
#        text.font = c(2,1,1,1))
# legend('bottomright',c('2008','2010','2011', '2013', '2016'), pch = c(19), col = j_brew_colors[c(1,2,3,4,5)], bty='n')
# dev.off()

## same plot but with ggplot
ggplot(Bore.algae, aes(N,bore))+
  geom_point(aes(color = Year, shape = Site, size = 1))+
  geom_smooth(method = "lm", color = 'black')+
  geom_errorbar(aes(x = N, ymin = bore - bore.se, ymax = bore+bore.se))+
  geom_errorbarh(aes(xmin = N-se.N, xmax = N+se.N, y = bore))+
  xlab('% Tissue N')+
  ylab(expression(paste('Mean density of', italic(' Lithophaga'),' (counts per cm'^{2},')')))+
  scale_shape_manual(values=c(13, 15, 16))+
  scale_color_brewer(palette="Set2")+
  guides(size=FALSE)+
  ggsave(filename = 'Output/Figure2.pdf', width = 6)

# check for normality of residuals
qqnorm(resid(N.mod))
qqline(resid(N.mod))
# looks good!

## TWo way ANOVAS for lithophaga densities, parrotfish biomass, and percent cover of massive porites by site and year
#Lithophaga
MBB.mod<-lm(bore.cm2^(1/4)~Site*Year , data = TSData) # need to 4th root transform to meet assumptions
anova(MBB.mod)
summary(MBB.mod)

qqnorm(resid(MBB.mod))
qqline(resid(MBB.mod))
hist(resid(MBB.mod))
# Not good... very zero-inflated and not normal

#scars
Scar.mod<-lm(log(bites.cm2+1)~Site*Year , data = TSData) 
anova(Scar.mod)
summary(Scar.mod)

qqnorm(resid(Scar.mod))
qqline(resid(Scar.mod))
hist(resid(Scar.mod))
# Not good... very zero-inflated and not normal

# run a hurdle model (Zero-inflated Gamma) because too many zeros and non-normal
TSData$bore<-ifelse(TSData$bore.cm2>0,1,0) # make a column of 1's and 0's 

MBB.mod1 <- glm(bore ~ Site*Year, data = TSData, family = binomial(link = logit)) # the logisitic part (0/1)
MBB.mod2 <- glm(bore.cm2 ~ Site*Year, data = subset(TSData, bore == 1), family = Gamma(link = log)) # The gamma part (>0)

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
Scar.mod2 <- glm(bites.cm2 ~ Site*Year, data = subset(TSData, bite == 1), family = Gamma(link = log)) # The gamma part (>0)
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

Scaridae.sum<-Scaridae %>% # sum across the transects for analysis
  group_by(Site, Year, Transect) %>%
  summarise(sum.count = sum(Count))

Scaridae.mod <-lm(log(sum.count+1)~Site*Year, data = Scaridae.sum)
anova(Scaridae.mod)
summary(Scaridae.mod)

qqnorm(resid(Scaridae.mod))
qqline(resid(Scaridae.mod))
hist(resid(Scaridae.mod))

# Porites cover
#Convert Porites SA to percent cover
Coralsummary <- CoralData %>% 
  group_by(site, year, pole) %>% 
  filter(year == 2006 |year == 2008| year == 2010|year == 2011|year == 2013|year == 2016)%>%
  summarise(mean = mean(massive, na.rm = TRUE)) # transect level average

#TSData$PoritesCover<-100*(TSData$Surface.area.m2/0.25)  # quadrats were 0.25 m2
Porites.mod<-lm(mean~site*year , data = Coralsummary)
#Porites.mod<-lm(PoritesCover^(1/4) ~Site*Year , data = TSData)
anova(Porites.mod)
summary(Porites.mod)

qqnorm(resid(Porites.mod))
qqline(resid(Porites.mod))
hist(resid(Porites.mod))


### Ground truth the 2D to 3D counts #####
# simple linear regression between the 2D and 3D counts
# scars
mod.truth.scars<-lm(log(dcount_scar+1)~log(pcount_scar+1), data = TruthData)
hist(resid(mod.truth.scars))  
qqnorm(resid(mod.truth.scars))
qqline(resid(mod.truth.scars))
anova(mod.truth.scars)  
summary(mod.truth.scars)

#Lithophaga
mod.truth.bore<-lm(log(dcount_mbb+1)~log(pcount_mbb+1), data = TruthData)
hist(resid(mod.truth.bore))  
qqnorm(resid(mod.truth.bore))
qqline(resid(mod.truth.bore))
anova(mod.truth.bore)  
summary(mod.truth.bore)

# plot the results
pdf(file = 'Output/SupplementalFig2.pdf', width = 10, height = 5, useDingbats = FALSE)
par(mfrow=c(1,2))
par(mar=c(5.1,6.3,4.1,2.1))
#plot the parrotfish scars
plot(log(TruthData$pcount_scar+1), log(TruthData$dcount_scar+1), xlim = c(0,6), cex.lab = 1.5, cex.axis = 1.5,
     ylim = c(0,6), pch = 19, xlab = 'log(Scar photo counts +1)', ylab = 'log(Scar diver counts +1)')
abline(0,1, lty = 2)

pred<-predict(mod.truth.scars, se.fit = TRUE)
ind<-order(TruthData$pcount_scar)
x<-log(TruthData$pcount_scar[ind]+1)
lines(x, pred$fit[ind])
lines(x, pred$fit[ind]+pred$se.fit[ind])
lines(x, pred$fit[ind]-pred$se.fit[ind])
polygon(c(x,rev(x)),c(pred$fit[ind]+pred$se.fit[ind],rev(pred$fit[ind]-pred$se.fit[ind])),col=grey2, border = NA)
legend(-0.75,6,'a)', bty='n', cex = 1.5)
# plot the lithophagids
plot(log(TruthData$pcount_mbb+1), log(TruthData$dcount_mbb+1), pch = 19, xlim = c(0,5), ylim = c(0,5),
     xlab = expression(paste('log(', italic('Lithophaga'),  ' photo counts +1)')), 
     ylab = expression(paste('log(', italic('Lithophaga'), ' diver counts +1)')), cex.lab = 1.5, cex.axis = 1.5)
abline(0,1, lty = 2)
pred<-predict(mod.truth.bore, se.fit = TRUE)
ind<-order(TruthData$pcount_mbb)
x<-log(TruthData$pcount_mbb[ind]+1)
lines(x, pred$fit[ind])
lines(x, pred$fit[ind]+pred$se.fit[ind])
lines(x, pred$fit[ind]-pred$se.fit[ind])
polygon(c(x,rev(x)),c(pred$fit[ind]+pred$se.fit[ind],rev(pred$fit[ind]-pred$se.fit[ind])),col=grey2, border = NA)
legend(-0.75,5,'b)', bty='n', cex = 1.5)
dev.off()

# means
TSData %>%
  group_by(Site) %>%
  summarise(bore.mean = mean(bore.cm2),
            bore.SE = sd(bore.cm2)/sqrt(n()),
            bite.mean = mean(bites.cm2),
            bite.SE = sd(bites.cm2)/sqrt(n()))


##supplemental plot of data over time
#scarids
#sumamrize means and SE
Scaridae.mean <- Scaridae %>%
  group_by(Site, Year, Transect) %>%
  summarise(sum.count = sum(Count))%>%
  group_by(Site,Year)%>%
  summarise(mean = mean(sum.count), se = sd(sum.count)/sqrt(n()))

# plot of Lithophaga densities across site and time
dodge<-position_dodge(width=0.5) # this offsets the points so they don't overlap

#plot
scaridae.plot <- ggplot(Scaridae.mean, aes(x = Year, y = mean, colour = Site, group = Site, fill = Site))+
  geom_line()+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width=0))+
  geom_point(size = 3, shape = 21)+
  xlab('Year')+
  ylab(expression(Parrotfish~"per 250"~{m}^2*~""))+
  theme_classic(base_size=12)+
  scale_fill_manual(values = c('black', 'gray47', 'white'))+
  scale_colour_manual(values = c("black", "gray47", "black"))+
  guides(colour = FALSE, fill = FALSE)
scaridae.plot

# lithophagids
#summarize means and SE
mbb.mean <- TSData %>%
  group_by(Site, Year)%>%
  summarise(mean.bore = mean(bore.cm2),se.bore = sd(bore.cm2)/sqrt(n()),
            mean.bite = mean(bites.cm2), se.bite = sd(bites.cm2)/sqrt(n()),
            mean.SA = mean(Surface.area.m2), se.SA = sd(Surface.area.m2)/sqrt(n()))

# plot
lithophaga.plot <- ggplot(mbb.mean, aes(x = Year, y = mean.bore, colour = Site, group = Site, fill = Site))+
  geom_line()+
  geom_errorbar(aes(ymin = mean.bore - se.bore, ymax = mean.bore + se.bore, width=0))+
  geom_point(size = 3, shape = 21)+
  xlab('Year')+
  ylab(expression(italic(Lithophaga)~"per"~{cm}^2*~""))+
  theme_classic(base_size=12)+
  theme(legend.justification=c(0,0), legend.position=c(0.6,0.55))+
  #theme(legend.background = element_rect(colour="black", size=0.5, linetype="solid"))+
  scale_fill_manual(values = c('black', 'gray47', 'white'))+
  scale_colour_manual(values = c("black", "gray47", "black"))
lithophaga.plot

# plot of percent cover of massive Porites over site and time
#
Coralsummary <- CoralData %>% 
  group_by(site, year) %>% 
  filter(year == 2006 |year == 2008| year == 2010|year == 2011|year == 2013|year == 2016)%>%
  summarise(mean = mean(massive, na.rm = TRUE),
            se = sd(massive)/sqrt(n()))%>%
  mutate(year = as.factor(year))

SA.plot <- ggplot(Coralsummary, aes(x = year, y = mean, group = site, colour = site, fill = site))+
  geom_line(position=position_dodge(width = 0.6))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width=0), 
                position=position_dodge(width = 0.6))+
  geom_point(size = 3, shape = 21, position=position_dodge(width = 0.6))+
  scale_fill_manual(name = "Site", values = c('black', 'gray47', 'white'))+
  scale_colour_manual(name = "Site", values = c("black", "black", "black"))+
  xlab('Year')+
  ylab(expression("Massive"~italic(Porites)~" (%)"))+
  ylim(0,20)+
  theme_classic(base_size=12)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))+
  theme(legend.position = c(0.8, 0.8))+
  guides(colour = FALSE, fill = FALSE)

#SA.plot <- ggplot(mbb.mean, aes(x = Year, y = mean.SA, colour = Site, group = Site, fill = Site))+
 # geom_line()+
  #geom_errorbar(aes(ymin = mean.SA - se.SA, ymax = mean.SA + se.SA, width=0))+
  #geom_point(size = 3, shape = 21)+
  #xlab('Year')+
  #ylab(expression("Massive"~italic(Porites)~"Cover (%)"))+
  #theme_classic(base_size=12)+
  #scale_fill_manual(values = c('black', 'gray47', 'white'))+
  #scale_colour_manual(values = c("black", "gray47", "black"))+
  #guides(colour = FALSE, fill = FALSE)
#SA.plot

# Scars
scar.plot <- ggplot(mbb.mean, aes(x = Year, y = mean.bite, colour = Site, group = Site, fill = Site))+
  geom_line()+
  geom_errorbar(aes(ymin = mean.bite - se.bite, ymax = mean.bite + se.bite, width=0))+
  geom_point(size = 3, shape = 21)+
  xlab('Year')+
  ylab(expression("Parrotfish scars"~"per"~{cm}^2*~""))+
  theme_classic(base_size=12)+
  scale_fill_manual(values = c('black', 'gray47', 'white'))+
  scale_colour_manual(values = c("black", "gray47", "black"))+
  guides(colour = FALSE, fill = FALSE)
scar.plot

supplemental <- plot_grid(lithophaga.plot, SA.plot, scar.plot, scaridae.plot, labels = c("A","B", "C", "D"), ncol = 2, align = 'v')
supplemental

ggsave("supplementalFig1.pdf", supplemental, path = "Output/", width = 6, height = 6, units = "in")

## look at the relationship between parrotfish densities and bite scars
 Scaridae.mean$Year<-as.integer(as.character(Scaridae.mean$Year)) # convert to integer
 all<-left_join(Scaridae.mean, mbb.mean) # join together 
 density_scars<-lm(all$mean.bite~all$mean) # look at the relationship between bites and scars
 anova(density_scars)

