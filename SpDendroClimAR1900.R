#Within and between analysis of species and sites - 1900-2000
rm(list=ls())
library(dplR)
library(bootRes)
library(dplyr)

setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Chronology Development/Tucson/3 - Species Levels")
myFiles <- list.files(pattern = "*.txt") 
myNames <- myFiles %>%
  substr(.,1,nchar(.)-7)
dat<-as.list(lapply(myFiles, read.rwl))
names(dat)<-myNames

#ADJUST END YEARS!
WilHil <- list('WIL_ABLA' = dat$WIL_ABLA, 'HIL_PIEN' = dat$HIL_PIEN,  'WIL_PIEN' = dat$WIL_PIEN) #Series to age adjust -2 years 
WilHil2 <- list('HIL_ABLA' = dat$HIL_ABLA) #Series to age adjust -3 years 

yr2<-lapply(WilHil, function(x){
  row.names(x)<-as.numeric(row.names(x))-2 #subtracts form end yr
  x[1:nrow(x)-1,]}) #Eliminates last row, incompletes

yr3<-lapply(WilHil2, function(x){
  row.names(x)<-as.numeric(row.names(x))-3
  x[1:nrow(x)-1,]})

dat$WIL_ABLA <- yr2$WIL_ABLA
dat$WIL_PIEN <- yr2$WIL_PIEN
dat$HIL_PIEN <- yr2$HIL_PIEN
dat$HIL_ABLA <- yr3$HIL_ABLA

#summary.stats <- lapply(dat, rwl.report)

det.chron <- lapply(dat, detrend, make.plot=FALSE, method=c("ModNegExp")) #Detrending series
chr.fin <- lapply(det.chron, chron, prefix = det.chron$x, biweight = TRUE, prewhiten = TRUE) #Creating chronology from detrended series; auto-refressive prewhittening 

yrs<-as.character(seq(1900, 1999, by=1))
chr.yrs<- lapply(chr.fin, function(x) {
  subset(x, rownames(x) %in% yrs, select=c(res, samp.depth))
}) %>%
  lapply(., function(x) {
    na.omit(x)
  })

###EPS: 0.85 is accepted threshold; search for EPS among detrended series *as opposed to final chronology*
# eps.cut <- 0.85 # EPS Cutoff
# eps.test <- lapply(det.chron, rwi.stats.running, window.length=20) %>%
#   lapply(., function(x) x[x$eps > 0.85, c("start.year")]) %>%
#   lapply(., function(x) x[c(1)]) %>%
#   as.data.frame(.) %>%
#   gather(., Chrono, EPSyear, FTA_ABLA:WIL_PIEN, factor_key=TRUE)

eps.1900 <- lapply(det.chron, rwi.stats.running, window.length=20) %>%
  lapply(., function(x) x[x$start.year > 1900,]) %>%
  lapply(., function(x) x[c(1,15)]) %>%
  do.call(rbind, . ) %>%
  mutate(Chron = row.names(.)) %>%
  mutate(Chron = substr(Chron, 1, 8))

###eps.test give the years where EPS starts being >0.85

#Read in climate data
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Chronology Development/Climate Data")

WILsea<- read.csv("WILsea.csv")
HILsea<- read.csv("HILsea.csv")
GSPsea<- read.csv("GSPsea.csv")
FTAsea<- read.csv("FTAsea.csv")
FTBsea<- read.csv("FTBsea.csv")
HUMsea<- read.csv("HUMsea.csv")
HWAsea<- read.csv("HWAsea.csv")
PMRsea<- read.csv("PMRsea.csv")
SSGsea<- read.csv("SSGsea.csv")

dev.off() #Prep graphics for plotting
names <- names(WILsea)[-1]

###
WIL.ABLA <- bootRes::dcc(chr.yrs$WIL_ABLA, WILsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "WIL",
         Species = "ABLA",
         Variable = names)
#dcplot(WIL.ABLA)

WIL.PIEN <- bootRes::dcc(chr.yrs$WIL_PIEN, WILsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "WIL",
         Species = "PIEN",
         Variable = names)
# dcplot(WIL.PIEN)

HIL.ABLA <- bootRes::dcc(chr.yrs$HIL_ABLA, HILsea, vnames=c("SEA"), method="corr", start = 1, end =12) %>%
  mutate(Site = "HIL",
         Species = "ABLA",
         Variable = names)
# dcplot(HIL.ABLA)

HIL.PIEN <- bootRes::dcc(chr.yrs$HIL_PIEN, HILsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "HIL",
         Species = "PIEN",
         Variable = names)
# dcplot(HIL.PIEN)

GSP.LALY <- bootRes::dcc(chr.yrs$GSP_LALY, GSPsea,  vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "GSP",
         Species = "LALY",
         Variable = names)
# dcplot(GSP.LALY)

GSP.ABLA <- bootRes::dcc(chr.yrs$GSP_ABLA, GSPsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "GSP",
         Species = "ABLA",
         Variable = names)
# dcplot(GSP.ABLA)

GSP.PIEN <- bootRes::dcc(chr.yrs$GSP_PIEN, GSPsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "GSP",
         Species = "PIEN",
         Variable = names)
# dcplot(GSP.PIEN)

FTA.ABLA <- bootRes::dcc(chr.yrs$FTA_ABLA, FTAsea,  vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "FTA",
         Species = "ABLA",
         Variable = names)
# dcplot(FTA.ABLA)

FTA.PIEN <- bootRes::dcc(chr.yrs$FTA_PIEN, FTAsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "FTA",
         Species = "PIEN",
         Variable = names)
# dcplot(FTA.PIEN)

FTA.LALY <- bootRes::dcc(chr.yrs$FTA_LALY, FTAsea,  vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "FTA",
         Species = "LALY",
         Variable = names)
# dcplot(FTA.LALY)

FTB.ABLA <- bootRes::dcc(chr.yrs$FTB_ABLA, FTBsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "FTB",
         Species = "ABLA",
         Variable = names)
# dcplot(FTB.ABLA)

FTB.LALY <- bootRes::dcc(chr.yrs$FTB_LALY, FTBsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "FTB",
         Species = "LALY",
         Variable = names)
# dcplot(FTB.LALY)

FTB.PIEN <- bootRes::dcc(chr.yrs$FTB_PIEN, FTBsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
mutate(Site = "FTB",
       Species = "PIEN",
       Variable = names)
#dcplot(FTB.PIEN) #Can't run, not long enough

HUM.ABLA <- bootRes::dcc(chr.yrs$HUM_ABLA, HUMsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "HUM",
         Species = "ABLA",
         Variable = names)
# dcplot(HUM.ABLA) #CORRECT start should be 21

HUM.PIEN <- bootRes::dcc(chr.yrs$HUM_PIEN, HUMsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "HUM",
         Species = "PIEN",
         Variable = names)
# dcplot(HUM.PIEN) #CORRECT start should be 27

HWA.ABLA <- bootRes::dcc(chr.yrs$HWA_ABLA, HWAsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "HWA",
         Species = "ABLA",
         Variable = names)
# dcplot(HWA.ABLA) #CORRECT start should be 20

HWA.PIEN <- bootRes::dcc(chr.yrs$HWA_PIEN, HWAsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
mutate(Site = "HWA",
       Species = "PIEN",
       Variable = names)
#dcplot(HWA.PIEN) #CORRECT is 16:33; can't run

HWA.PIAL <- bootRes::dcc(chr.yrs$HWA_PIAL, HWAsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "HWA",
         Species = "PIAL",
         Variable = names) #CORRECT should be 42:69
# dcplot(HWA.PIAL)

PMR.ABLA <- bootRes::dcc(chr.yrs$PMR_ABLA, PMRsea,  vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "PMR",
         Species = "ABLA",
         Variable = names)
# dcplot(PMR.ABLA)

PMR.PIAL <- bootRes::dcc(chr.yrs$PMR_PIAL, PMRsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "PMR",
         Species = "PIAL",
         Variable = names)
# dcplot(PMR.PIAL)

PMR.PIEN <- bootRes::dcc(chr.yrs$PMR_PIEN, PMRsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "PMR",
         Species = "PIEN",
         Variable = names)
# dcplot(PMR.PIEN)

SSG.ABLA <- bootRes::dcc(chr.yrs$SSG_ABLA, SSGsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "SSG",
         Species = "ABLA",
         Variable = names) #Should be start 73 
# dcplot(SSG.ABLA)

SSG.PIAL <- bootRes::dcc(chr.yrs$SSG_PIAL, SSGsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "SSG",
         Species = "PIAL",
         Variable = names) #Should be start 63
# dcplot(SSG.PIAL)

SSG.PIEN <-bootRes::dcc(chr.yrs$SSG_PIEN, SSGsea, vnames=c('SEA'), method="corr", start = 1, end =12) %>%
  mutate(Site = "SSG",
         Species = "PIEN",
         Variable = names)
# dcplot(SSG.PIEN)

##MARCH 2020: Change analysis to Tave and fix proportion calculations 
##################
sp.list <- list("WIL.ABLA" = WIL.ABLA, "WIL.PIEN" = WIL.PIEN, "HIL.ABLA" = HIL.ABLA, "HIL.PIEN" = HIL.PIEN, "GSP.LALY"= GSP.LALY, "GSP.ABLA" = GSP.ABLA, "GSP.PIEN" = GSP.PIEN, "FTA.ABLA" = FTA.ABLA, "FTA.PIEN" = FTA.PIEN, "FTA.LALY" = FTA.LALY, "FTB.ABLA" = FTB.ABLA, "FTB.LALY" = FTB.LALY, "FTB.PIEN" = FTB.PIEN, "HUM.ABLA" = HUM.ABLA, "HUM.PIEN" = HUM.PIEN, "HWA.ABLA" = HWA.ABLA, "HWA.PIAL" = HWA.PIAL, "HWA.PIEN" = HWA.PIEN, "PMR.ABLA" = PMR.ABLA, "PMR.PIAL" = PMR.PIAL, "PMR.PIEN" = PMR.PIEN, "SSG.ABLA" = SSG.ABLA, "SSG.PIAL"= SSG.PIAL, "SSG.PIEN" = SSG.PIEN) #List the clim.dcc

var.order <- WIL.ABLA[,7]

dat <- Reduce(rbind, sp.list) %>%
  mutate(Sig = as.numeric(ifelse(significant == TRUE, 1, 0))) %>%
  mutate(Variable = factor(Variable, levels = (var.order)),
         Direction = ifelse(coef > 0, "Positive", "Negative")) %>%
  group_by(Variable, Site, Species, Direction) %>%
  dplyr::summarise(Freq = sum(Sig)) %>%
  dplyr::mutate(Prop = (if_else(Species == "ABLA", Freq/9, if_else(Species == "PIEN", Freq/9, if_else(Species == "LALY", Freq/3, Freq/3)))),
                PltFreq = ifelse(Direction == c("Negative"), Freq*-1, Freq*1)) %>%
  ungroup() %>%
  mutate(Variable = factor(Variable, levels = c("PPT_wt","PPT_sp","PPT_sm","PPT_at","Tmax_wt", "Tmax_sp", "Tmax_sm","Tmax_at","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at")))


a<-rep(c('bold'), 4) #Bolds for differentiating vars
b<-rep(c('plain'), 4)

theme_emma <- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12, family = "Times"),
    axis.text.x = element_text(angle = 45, hjust = 1, face=c(a,b)),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size = 12, face = "bold")
  )}


mypal <- c("#4B4B4B","#A2A5A5", "#717171", '#E2E2E2')

sig.plot <-ggplot(data=dat, aes(x=Variable, y=PltFreq/24, fill = Site))+
  geom_bar(stat="identity", width=0.5, color = "black", size = 0.3)+
  facet_wrap(~Species)+
  #ylim(-0.5, 0.5)+
  #scale_fill_manual(values = mypal)+
  xlab('Seasonal climate variables')+
  ylab("Proportion of chronologies")+
  # ggtitle("Proportion of total chronologies significantly correlated \nto seasonal climate variables") +
  theme_bw()+
  theme_emma()
  #theme(legend.title=element_blank(), legend.justification=c(0,0), legend.position=c(.01,.6))

sig.plot

#ggsave(plot=sig.plot, "TRWClim2.pdf", device = "pdf", width = 4, height = 4, units = c("in"), dpi= 600)

mypal <- c("#4B4B4B","#A2A5A5",'#E2E2E2', "#717171")
a<-rep(c('bold'), 4) #Bolds for differentiating vars
b<-rep(c('plain'), 4)



sig.plot <-ggplot(data=dat.sum, aes(x=Variable, y=CountDir, fill = Dir))+
  geom_bar(stat="identity", width=0.5, color = "black", size = 0.3)+
  xlab('Seasonal climate variables')+
  ylab("n Sites x direction of correlation")+
  scale_fill_manual(values = mypal, name="Direction of Correlation")+
  guides(fill=FALSE)+
  theme_bw()+
  ylim(-4, 4)+
  theme_emma()

sig.plot
