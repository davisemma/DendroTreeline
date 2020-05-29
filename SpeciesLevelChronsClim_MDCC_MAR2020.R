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
HWA2 <- list('HWA_ABLA' = dat$HWA_ABLA, 'HWA_PIEN' = dat$HWA_PIEN, 'HWA_PIAL' = dat$HWA_PIAL)

yr2<-lapply(WilHil, function(x){
  row.names(x)<-as.numeric(row.names(x))-2 #subtracts form end yr
  x[1:nrow(x)-1,]}) #Eliminates last row, incompletes

yr3<-lapply(WilHil2, function(x){
  row.names(x)<-as.numeric(row.names(x))-3
  x[1:nrow(x)-1,]})

yr1<-lapply(HWA2, function(x){
  row.names(x)<-as.numeric(row.names(x))-1
  x[1:nrow(x)-1,]})

#Replace with adjusted data
dat$WIL_ABLA <- yr2$WIL_ABLA
dat$WIL_PIEN <- yr2$WIL_PIEN
dat$HIL_PIEN <- yr2$HIL_PIEN
dat$HIL_ABLA <- yr3$HIL_ABLA
dat$HWA_ABLA <- yr1$HWA_ABLA
dat$HWA_PIEN <- yr1$HWA_PIEN
dat$HWA_PIAL <- yr1$HWA_PIAL

#summary.stats <- lapply(dat, rwl.report)

det.chron <- lapply(dat, detrend, make.plot=FALSE, method=c("ModNegExp")) #Detrending series
chr.fin <- lapply(det.chron, chron, prefix = det.chron$x, biweight = TRUE, prewhiten = TRUE) #Creating chronology from detrended series; auto-refressive prewhittening 

yrs<-as.character(seq(1900, 2015, by=1))
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
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Chronology Development/Climate Data/Chrono Climate Data MAR 2020")

WILsea<- read.csv("WILsea_2020.csv")
HILsea<- read.csv("HILsea_2020.csv")
GSPsea<- read.csv("GSPsea_2020.csv")
FTAsea<- read.csv("FTAsea_2020.csv")
FTBsea<- read.csv("FTBsea_2020.csv")
HUMsea<- read.csv("HUMsea_2020.csv")
HWAsea<- read.csv("HWAsea_2020.csv")
PMRsea<- read.csv("PMRsea_2020.csv")
SSGsea<- read.csv("SSGsea_2020.csv")

dev.off() #Prep graphics for plotting
names <- names(WILsea)[-1]

#MDCC analyses

WIL.ABLA <- bootRes::mdcc(chr.yrs$WIL_ABLA, WILsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

WIL.PIEN <- bootRes::mdcc(chr.yrs$WIL_PIEN, WILsea,
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

HIL.ABLA <- bootRes::mdcc(chr.yrs$HIL_ABLA, HILsea,  
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

HIL.PIEN <- bootRes::mdcc(chr.yrs$HIL_PIEN, HILsea,
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

GSP.LALY <- bootRes::mdcc(chr.yrs$GSP_LALY, GSPsea, 
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

GSP.ABLA <- bootRes::mdcc(chr.yrs$GSP_ABLA, GSPsea, 
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

GSP.PIEN <- bootRes::mdcc(chr.yrs$GSP_PIEN, GSPsea, 
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

FTA.ABLA <- bootRes::mdcc(chr.yrs$FTA_ABLA, FTAsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

FTA.PIEN <- bootRes::mdcc(chr.yrs$FTA_PIEN, FTAsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

FTA.LALY <- bootRes::mdcc(chr.yrs$FTA_LALY, FTAsea, 
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

FTB.ABLA <- bootRes::mdcc(chr.yrs$FTB_ABLA, FTBsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

FTB.LALY <- bootRes::mdcc(chr.yrs$FTB_LALY, FTBsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

FTB.PIEN <- bootRes::mdcc(chr.yrs$FTB_PIEN, FTBsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)
#dcplot(FTB.PIEN) #Can't run, not long enough

HUM.ABLA <- bootRes::mdcc(chr.yrs$HUM_ABLA, HUMsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)
# dcplot(HUM.ABLA) #CORRECT start should be 21

HUM.PIEN <- bootRes::mdcc(chr.yrs$HUM_PIEN, HUMsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

HWA.ABLA <- bootRes::mdcc(chr.yrs$HWA_ABLA, HWAsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

HWA.PIEN <- bootRes::mdcc(chr.yrs$HWA_PIEN, HWAsea, 
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

HWA.PIAL <- bootRes::mdcc(chr.yrs$HWA_PIAL, HWAsea,  
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

PMR.ABLA <- bootRes::mdcc(chr.yrs$PMR_ABLA, PMRsea,  
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

PMR.PIAL <- bootRes::mdcc(chr.yrs$PMR_PIAL, PMRsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

PMR.PIEN <- bootRes::mdcc(chr.yrs$PMR_PIEN, PMRsea,  
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

SSG.ABLA <- bootRes::mdcc(chr.yrs$SSG_ABLA, SSGsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

SSG.PIAL <- bootRes::mdcc(chr.yrs$SSG_PIAL, SSGsea, 
                          vnames= NULL, method="corr", 
                          win.size = 30, start = 1, end = 8, 
                          win.offset =5)

SSG.PIEN <-bootRes::mdcc(chr.yrs$SSG_PIEN, SSGsea,  
                         vnames= NULL, method="corr", 
                         win.size = 30, start = 1, end = 8, 
                         win.offset =5)

 
##################
sp.list <- list("WIL.ABLA" = WIL.ABLA$coef, "WIL.PIEN" = WIL.PIEN$coef, "HIL.ABLA" = HIL.ABLA$coef, 
                "HIL.PIEN" = HIL.PIEN$coef, "GSP.LALY"= GSP.LALY$coef, "GSP.ABLA" = GSP.ABLA$coef,
                "GSP.PIEN" = GSP.PIEN$coef, "FTA.ABLA" = FTA.ABLA$coef, "FTA.PIEN" = FTA.PIEN$coef,
                "FTA.LALY" = FTA.LALY$coef, "FTB.ABLA" = FTB.ABLA$coef, "FTB.LALY" = FTB.LALY$coef,
                "FTB.PIEN" = FTB.PIEN$coef, "HUM.ABLA" = HUM.ABLA$coef, "HUM.PIEN" = HUM.PIEN$coef,
                "HWA.ABLA" = HWA.ABLA$coef, "HWA.PIAL" = HWA.PIAL$coef, "HWA.PIEN" = HWA.PIEN$coef,
                "PMR.ABLA" = PMR.ABLA$coef, "PMR.PIAL" = PMR.PIAL$coef, "PMR.PIEN" = PMR.PIEN$coef,
                "SSG.ABLA" = SSG.ABLA$coef, "SSG.PIAL" = SSG.PIAL$coef, "SSG.PIEN" = SSG.PIEN$coef) #List the clim.dcc

var.order <- WIL.ABLA[,7]



###
theme_emma <- function(){  
  theme(
    axis.text = element_text(size =8),
    text = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=8, face="bold"),
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )}

vars <- names(WILsea[2:9])

test <- lapply(sp.list, mutate, Variable = vars) %>%
  lapply(., melt, id = 'Variable') %>%
  lapply(., setNames, nm = c('Variable', 'Years', 'r'))

test.melt <- melt(test, id = c('Variable', 'Years', 'r')) %>%
  mutate(Site = substr(L1, 1, 3),
         Species = substr(L1, 5,8))

ggplot(filter(test.melt, Site == 'GSP'), aes(x = Years, y = Variable))+
  geom_raster(aes(fill = r))+
  scale_fill_gradient2(low = "blue", high = "red", 
                       midpoint = 0, 
                       limit = c(-1, 1))+
  facet_wrap(~Species, scales = 'free')+
  theme_few()+
  theme_emma()

#ggsave(plot=sig.plot, "~/Desktop/DendroChrons.pdf", device = "pdf", width = 8.5, height = 4, units = c("in"), dpi= 600)
