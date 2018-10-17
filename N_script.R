library(tidyverse)

NutData<-read.csv('Data/Macroalgal CHN.csv')

Turb<-NutData %>%
  filter(Habitat =='Fringe', Genus =='Turbinaria', 
         Site == 'LTER 1'| Site == 'LTER 3' | Site=="LTER 4")

Sarg<-NutData %>%
  filter(Habitat =='Fringe', Genus =='Sargassum', 
         Site == 'LTER 1'| Site == 'LTER 3' | Site=="LTER 4") %>%
      group_by(Site,Year) %>%
        summarise(N.sarg = mean(N, na.rm=T))



#mean and SD
mean(Turb$N, na.rm=T)
sd(Turb$N, na.rm=T)

Turb.mean<-Turb %>%
  group_by(Site) %>%
  summarise(N = mean(N, na.rm=T))


bore<-TSData %>%
  filter(Site == 'LTER1'| Site == 'LTER3' | Site=="LTER4")%>%
  group_by(Site, Year) %>%
  summarise(bore = mean(bore.cm2),  bore.se = sd(bore.cm2)/sqrt(n()))

Turb.YS<-Turb %>%
  filter(Site == 'LTER 1'| Site == 'LTER 3' | Site=="LTER 4")%>%
  group_by(Site, Year) %>%
  summarise(N = mean(N, na.rm=T),  N.se = sd(N, na.rm=T)/sqrt(n())) 
  
#rename the levels
levels(bore$Site)<-c("LTER 1","LTER 2","LTER 3","LTER 4","LTER 5","LTER 6")

Bore.alge<-left_join(bore, Turb.YS)



barplot(bore$bore)

plot(Turb.mean$N, bore$bore[c(1,3,4)])

bore$Nuts<-ifelse(bore$Site=='LTER4', 'low','high')

a<-lm(bore$bore~bore$Nuts)


