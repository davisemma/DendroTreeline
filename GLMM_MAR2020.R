##THIS CODE: Is a GLMM to model tree establishment using climate variables ---------------------------
##LAST UPDATE: March 16, 2020
##In this version, run w full dataset (bins up to 2010)

rm(list=ls())
packs <- c('car','ggplot2','ggthemes','MuMIn','dplyr', 'reshape2', 'readr', 'tidyr', 'lme4')
lapply(packs, library, character.only = TRUE)
detach("package:dplyr", unload=TRUE)
library("dplyr")

#Function for testing overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data")

samples <- read.csv("REGION_Samples.csv")

#Giving bin dats to est. years 
breaks<- seq(1590, 2015, by=5) #Create breaks
bin<-as.array((cut(samples$EstYear, breaks, right=FALSE))) #give each est year a bin
bins<-as.numeric(substring(bin,2,5)) #Record the bin value
samples<-(cbind(samples,bins)) #Combine

#Format data - Set levels and select species and time
dat <- samples %>%
  mutate(SiteID = factor(samples$Site, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM')),
        Maturity = factor(samples$Maturity, levels = c('SEE', 'SAP', 'MAT', 'KRU')),
        Species = factor(samples$Species, levels=c('ABLA',  'PIEN',  'LALY',  'ASPEN', 'PIAL',  'MISL'))) %>%
  dplyr::filter(bins >= 1900 & bins < 2015) %>%
  dplyr::filter(Species != 'ASPEN') %>%
  dplyr::filter(Species != 'MISL') %>%
  group_by(SiteID, bins, Species) %>%
  summarise(Count = dplyr::n())

#Creating formatted dataframe and merging w frequencies
bins <- unique(dat$bins)

bins
####
full.samp <- as.data.frame(cbind(SiteID = rep(c('WIL', 'HIL','GSP', 'FTB','FTA', 'HWA', 'SSG','PMR', 'HUM'), each = 23*4), bins = rep(bins, 4*9), Species = rep(rep(c('ABLA', 'PIEN', 'LALY', 'PIAL'), each = 23))))

sp.dat <- merge(dat, full.samp, all.y=TRUE)
sp.dat[is.na(sp.dat)] <- 0

####

#Bringing in climate data
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Climate WNA/Seasonal/")
FTA<-read.csv('5yrFTA.csv')[1:23,]
FTB<-read.csv('5yrFTB.csv')[1:23,]
HUM<-read.csv('5yrHUM.csv')[1:23,]
HWA<-read.csv('5yrHWA.csv')[1:23,]
PMR<-read.csv('5yrPMR.csv')[1:23,]
SSG<-read.csv('5yrSSG.csv')[1:23,]
HIL<-read.csv('5yrHIL.csv')[1:23,]
WIL<-read.csv('5yrWIL.csv')[1:23,]
GSP<-read.csv('5yrGSP.csv')[1:23,]

#Read in establishment data (for formatting)
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data")
estdata <- read_csv("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data/REGION5yrFreqBins.csv")

dfs<-list('WIL'=WIL,'HIL'=HIL, 'GSP'=GSP, 'FTA'=FTA, 'FTB'=FTB,'HUM'=HUM,'HWA'=HWA, 'PMR'=PMR, 'SSG'=SSG)

Tave_wt<- cbind(aveBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  Tave_wt<-x$Tave_wt[1:23])))
Tave_wt<-melt(Tave_wt, 'aveBin', value.name = "Tave_wt")

Tave_sp<- cbind(aveBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  Tave_sp<-x$Tave_sp[1:23])))
Tave_sp<-melt(Tave_sp, 'aveBin', value.name = "Tave_sp")
          
Tave_sm<- cbind(aveBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  Tave_sm<-x$Tave_sm[1:23])))
Tave_sm<-melt(Tave_sm, 'aveBin', value.name = "Tave_sm")

Tave_at<- cbind(aveBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  Tave_at<-x$Tave_at[1:23])))
Tave_at<-melt(Tave_at, 'aveBin', value.name = "Tave_at")

PAS_wt<- cbind(MinBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  PAS_wt<-x$PAS_wt[1:23])))
PAS_wt<-melt(PAS_wt, 'MinBin', value.name = "PAS_wt")

PAS_sp<- cbind(MinBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  PAS_sp<-x$PAS_sp[1:23])))
PAS_sp<-melt(PAS_sp, 'MinBin', value.name = "PAS_sp")

PPT_sm<- cbind(MinBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  PPT_sm<-x$PPT_sm[1:23])))
PPT_sm<-melt(PPT_sm, 'MinBin', value.name = "PPT_sm")

DD5_sm<- cbind(MinBin= seq(1900,2010, 5),as.data.frame(lapply(dfs, function(x) 
  DD5_sm<-x$DD5_sm[1:23])))
DD5_sm<-melt(DD5_sm, 'MinBin', value.name = "DD5_sm")

##
mydata <- estdata %>%
  filter(MinBin >= 1900 & MinBin <2015)%>%
  gather(Site, Frequency, WIL:SSG) %>%
  mutate(Slope = ifelse (Site=='FTA', 28,ifelse (Site=='FTB', 16,ifelse (Site=='HWA', 23,ifelse (Site=='SSG',27, ifelse (Site=='PMR', 27,ifelse (Site=='HIL', 10,ifelse (Site=='WIL', 18,ifelse (Site=='GSP', 28,ifelse (Site=='HUM', 25, NA))))))))),
         Aspect = ifelse (Site=='FTA', 107,ifelse (Site=='FTB', 273, ifelse (Site=='HWA', 177,ifelse (Site=='SSG', 160,ifelse (Site=='PMR', 225, ifelse (Site=='HIL', 45,ifelse (Site=='WIL', 208, ifelse (Site=='GSP', 197,ifelse (Site=='HUM', 25, NA))))))))),
         Latitude = ifelse (Site=='FTA', 50.8295,ifelse (Site=='FTB', 50.8359,ifelse (Site=='HWA', 50.571,ifelse (Site=='SSG',50.235, ifelse (Site=='PMR', 50.226,ifelse (Site=='HIL', 52.195,ifelse (Site=='WIL', 52.226,ifelse (Site=='GSP', 51.2, ifelse (Site=='HUM',50.200, NA))))))))),
         Site = factor(Site, levels = c("WIL","HIL","GSP","FTA","FTB","HUM","HWA","PMR","SSG")))

freq.dat<-as.data.frame(cbind(mydata,
                            Tave_wt=c(Tave_wt$Tave_wt),
                            Tave_sp=c(Tave_sp$Tave_sp),
                            Tave_sm=c(Tave_sm$Tave_sm),
                            Tave_at=c(Tave_at$Tave_at),
                            PAS_wt=c(PAS_wt$PAS_wt),
                            PAS_sp=c(PAS_sp$PAS_sp),
                            PPT_sm=c(PPT_sm$PPT_sm),
                            DD5_sm=c(DD5_sm$DD5_sm)))%>%
  mutate(SiteID = Site,
         bins = MinBin)

#Creating columns for HLI
freq.dat$HLI2<-((1-(cos(freq.dat$Aspect)-45))/2)
freq.dat$fAspect<-(pi/180)*(180 - abs(freq.dat$Aspect - 225))
freq.dat$fSlope<-(pi/180)*freq.dat$Slope
freq.dat$fLatitude<-(pi/180)*freq.dat$Latitude

#Heat load index based on McCune(2002)
freq.dat$HLI<- exp(-1.467+(1.582*cos(freq.dat$fLatitude)*cos(freq.dat$fSlope))-(1.5* cos(freq.dat$fAspect)*sin(freq.dat$fSlope)*sin(freq.dat$fLatitude))-(0.262*sin(freq.dat$fLatitude)*sin(freq.dat$fSlope))+(0.607*sin(freq.dat$fSlope)*sin(freq.dat$fAspect)))

myvars <- c("Tave_wt", "Tave_sp", "Tave_sm", "Tave_at", "PAS_wt", "PAS_sp", "PPT_sm", "DD5_sm", "SiteID","bins","HLI")
sel.dat <- freq.dat[myvars]

#MERGE clim dat w species level freq dat
mod.dat <- merge(sel.dat, sp.dat)

#Reading in subsample WILHIL Data
div.dat <- read_csv("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data/WILHIL_GLMMsampledat.csv")

test <- mod.dat %>%
  filter(SiteID == 'WIL' | SiteID == 'HIL') 
  
test.merge <- merge(div.dat, test, by = c('SiteID', 'bins'), all.y = TRUE) %>%
  mutate(Scaled = Count/Divis,
         Count = Scaled)

myvars2 <- c("Tave_wt", "Tave_sp", "Tave_sm", "Tave_at", "PAS_wt", "PAS_sp", "PPT_sm", "DD5_sm", "SiteID","bins","HLI","Species", "Count")
test.merge2 <-test.merge[myvars2] %>%
  mutate(SiteID = ifelse(SiteID == 'HIL', 'HIL2', ifelse(SiteID=='WIL', 'WIL2', NA)))

mod.dat2 <- rbind(mod.dat, test.merge2) %>%
  mutate(Count = round(Count)) %>%
  filter(SiteID != 'WIL' & SiteID != 'HIL')
mod.dat2[,3:11] <- scale(mod.dat2[,3:11])


#ALL SPECIES TOGETHER
mod.test <-glmer(Count ~ Tave_sm + PAS_wt + PAS_sp + PPT_sm  + Tave_sm*PPT_sm + (1|Species) + (1|bins) + (1|SiteID), family=poisson, data=filter(mod.dat2))

summary(mod.test)

abla.dat <- filter(mod.dat2, Species == 'ABLA')
pien.dat <- filter(mod.dat2, Species == 'PIEN')
laly.dat <- filter(mod.dat2, SiteID == 'FTA' | SiteID == 'FTB' | SiteID == 'GSP') %>%
  filter(., Species == 'LALY')
pial.dat <- filter(mod.dat2, SiteID == 'HWA' | SiteID == 'SSG' | SiteID == 'PMR') %>%
  filter(., Species == 'PIAL')

#MODEL OPTIONS

#ABLA MODEL
amod1<-glmer(Count ~  Tave_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod2<-glmer(Count ~  PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod3<-glmer(Count ~  PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod4<-glmer(Count ~  Tave_sm + PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod5<-glmer(Count ~  Tave_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod6<-glmer(Count ~  PPT_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod7<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod8<-glmer(Count ~  Tave_sm + PPT_sm  + Tave_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod10<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod9<-glmer(Count ~  Tave_sm  + PAS_wt + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod11<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))

amod12<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'ABLA'))


anova(amod1, amod2, amod3, amod4, amod5, amod6, amod7, amod8, amod9, amod10, amod11, amod12)

summary(amod12)

qqnorm(resid(amod12))
qqline(resid(amod12))
hist(resid(amod12))
shapiro.test(resid(amod12))
plot(fitted(amod12),resid(amod12))
acf(resid(amod12))


#PIEN MODEL
pmod1<-glmer(Count ~  Tave_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod2<-glmer(Count ~  PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod3<-glmer(Count ~  PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod4<-glmer(Count ~  Tave_sm + PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod5<-glmer(Count ~  Tave_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod6<-glmer(Count ~  PPT_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod7<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod8<-glmer(Count ~  Tave_sm + PPT_sm  + Tave_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod10<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod9<-glmer(Count ~  Tave_sm  + PAS_wt + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod11<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

pmod12<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + Tave_sm*PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=filter(mod.dat2, Species == 'PIEN'))

anova(pmod1, pmod2, pmod3, pmod4, pmod5, pmod6, pmod7, pmod8, pmod9, pmod10, pmod11, pmod12)
summary(pmod12)

##TESTING TESTING
pamod1<-glmer(Count ~  Tave_sm + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod2<-glmer(Count ~  PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod3<-glmer(Count ~  PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod4<-glmer(Count ~  Tave_sm + PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod5<-glmer(Count ~  Tave_sm + PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod6<-glmer(Count ~  PPT_sm + PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod7<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod8<-glmer(Count ~  Tave_sm + PPT_sm  + Tave_sm*PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod10<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod9<-glmer(Count ~  Tave_sm  + PAS_wt + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod11<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

pamod12<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=pial.dat)

anova(pamod1, pamod2, pamod3, pamod4, pamod5, pamod6, pamod7, pamod8, pamod9, pamod10, pamod11)

summary(pamod1)

##
lmod1<-glmer(Count ~  Tave_sm + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod2<-glmer(Count ~  PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod3<-glmer(Count ~  PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod4<-glmer(Count ~  Tave_sm + PPT_sm + (1|SiteID) + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod5<-glmer(Count ~  Tave_sm + PAS_wt + (1|SiteID) + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod6<-glmer(Count ~  PPT_sm + PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

#lmod7<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod8<-glmer(Count ~  Tave_sm + PPT_sm  + Tave_sm*PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

#lmod10<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod9<-glmer(Count ~  Tave_sm  + PAS_wt + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

lmod11<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)
 
#lmod12<-glmer(Count ~  Tave_sm + PPT_sm + PAS_wt + Tave_sm*PPT_sm + Tave_sm*PAS_wt + SiteID + (1|bins), family=poisson(link='log'), data=laly.dat)

anova(lmod1, lmod2, lmod3, lmod4, lmod5, lmod6, lmod8, lmod9, lmod11)

summary(lmod9)


#Check for dispersion (deviance/df.resid), should be less than 1, definitely not greater than 5
summary(mod.test)
809.1/153

qqnorm(resid(mod.test))
qqline(resid(mod.test))
hist(resid(mod.test))
shapiro.test(resid(mod.test))
plot(fitted(mod.test),resid(mod.test))
acf(resid(mod.test))

###############

r.squaredGLMM(abla.mod)
x <- filter(mod.dat2, SiteID !='WIL' & SiteID !='HIL')
hist(x$Count, breaks = 50)

a.mod1<-glmer(Count ~ PAS_wt + Tmin_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson, data=filter(mod.dat2, Species == 'ABLA' & SiteID != 'WIL' & SiteID != 'HIL'))
summary(a.mod1)

a.mod2<-glmer(Count ~ Tmin_sm + Tmin_at + PAS_wt + PAS_sp + PPT_sm  + Tmin_sm*PPT_sm + (1|SiteID) + (1|bins), family=poisson, data=filter(mod.dat2, Species == 'ABLA' & SiteID != 'WIL' & SiteID != 'HIL'))
summary(a.mod2)

###
exp.dat <- pien.dat %>%
  mutate(Tmin_smPAS_wt = Tmin_sm*PAS_wt,
         Tmin_smPPT_sm = Tmin_sm*PPT_sm)

plot(exp.dat$Tmin_smPAS_wt,exp.dat$Count)
plot(exp.dat$PAS_wt,exp.dat$Count)
plot(exp.dat$Tmin_sm,exp.dat$Count)
