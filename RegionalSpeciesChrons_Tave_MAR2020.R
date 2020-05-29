#THIS CODE: Produces regional species chronologies and conducts a 
#dendroclim analysis using seasonal climate data

rm(list=ls())
pkgs <- c("dplR", "bootRes", "dplyr", "ggsci", "ggthemes", "scales")
lapply(pkgs, require, character.only=TRUE)

#Reading in text files
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Chronology Development/Tucson/3 - Species Levels")
#Read only .txt extentions
myFiles <- list.files(pattern = "*.txt")  
#Get file names without extention
myNames <- myFiles %>%
  substr(.,1,nchar(.)-7)
#Read in as 'dat' and rename
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

#Summary statistics for each .rwl
#summary.stats <- lapply(dat, rwl.report)

###CHRONOLOGY BUILDING
#Detrending series
det.chron <- lapply(dat, detrend, make.plot=FALSE, method=c("ModNegExp"))

#Creating chronology from detrended series; auto-refressive prewhittening 
chr.fin <- lapply(det.chron, chron, prefix = det.chron$x, biweight = TRUE, prewhiten = TRUE)

#REGIONAL ABLA CHRONOLOGY
abla<-combine.rwl(chr.fin$FTA_ABLA, chr.fin$FTB_ABLA) %>%
  combine.rwl(., chr.fin$GSP_ABLA) %>%
  combine.rwl(., chr.fin$HIL_ABLA) %>%
  combine.rwl(., chr.fin$HUM_ABLA[14:71,]) %>%
  combine.rwl(., chr.fin$HWA_ABLA) %>%
  combine.rwl(., chr.fin$PMR_ABLA) %>%
  combine.rwl(., chr.fin$SSG_ABLA) %>%
  combine.rwl(., chr.fin$WIL_ABLA) 

#abla1<-abla[143:242, c(1,4,7,10,13,16,19,22,25)]
abla1<-abla[143:258, c(2,5,8,11,14,17,20,23,26)] #RESID CHRONOLOGY
colnames(abla1)<-c('FTA_ABLA', 'FTB_ABLA', "GSP_ABLA", "HIL_ABLA", "HUM_ABLA", "HWA_ABLA", "PMR_ABLA", "SSG_ABLA", "WIL_ABLA")
rwl.report(abla1)
abla.chrn <- chron(abla1, prefix = 'ABL', biweight = TRUE, prewhiten = FALSE)

#PIEN
pien<-combine.rwl(chr.fin$FTA_PIEN, chr.fin$FTB_PIEN) %>%
  combine.rwl(., chr.fin$GSP_PIEN) %>%
  combine.rwl(., chr.fin$HIL_PIEN) %>%
  combine.rwl(., chr.fin$HUM_PIEN[16:70,]) %>%
  combine.rwl(., chr.fin$HWA_PIEN) %>%
  combine.rwl(., chr.fin$PMR_PIEN) %>%
  combine.rwl(., chr.fin$SSG_PIEN) %>%
  combine.rwl(., chr.fin$WIL_PIEN) 

# pien1<-pien[244:343, c(1,4,7,10,13,16,19,22,25)]
pien1<-pien[244:359, c(2,5,8,11,14,17,20,23,26)]
colnames(pien1)<-c('FTA_PIEN', 'FTB_PIEN', "GSP_PIEN", "HIL_PIEN", "HUM_PIEN", "HWA_PIEN", "PMR_PIEN", "SSG_PIEN", "WIL_PIEN")
rwl.report(pien1)
pien.chrn <- chron(pien1, prefix = 'PIE', biweight = TRUE, prewhiten = FALSE)

#LALY
laly<-combine.rwl(chr.fin$FTA_LALY, chr.fin$FTB_LALY) %>%
  combine.rwl(., chr.fin$GSP_LALY)

#laly1<-laly[173:272, c(1,4,7)]
laly1<-laly[173:288, c(2,5,8)]
colnames(laly1)<-c('FTA_LALY', 'FTB_LALY', "GSP_LALY")
rwl.report(laly1)
laly.chrn <- chron(laly1, prefix = 'LAL', biweight = TRUE, prewhiten = FALSE)

#PIAL
pial<-combine.rwl(chr.fin$PMR_PIAL, chr.fin$HWA_PIAL) %>%
  combine.rwl(., chr.fin$SSG_PIAL)

# pial1<-pial[221:292, c(1,4,7)]
pial1<-pial[221:304, c(2,5,8)]
colnames(pial1)<-c('PMR_PIAL', 'HWA_PIAL', "SSG_PIAL")
rwl.report(pial1)
pial.chrn <- chron(pial1, prefix = 'PIA', biweight = TRUE, prewhiten = FALSE)

###

#PLOTS --------
#7 x 3
names(abla.chrn) = c("Chron", "nSamples")
abla.chrn$Species = c("ABLA")
names(pien.chrn) = c("Chron", "nSamples")
pien.chrn$Species = c("PIEN")
names(laly.chrn) = c("Chron", "nSamples")
laly.chrn$Species = c("LALY")
names(pial.chrn) = c("Chron", "nSamples")
pial.chrn$Species = c("PIAL")

plt.dat<-rbind(abla.chrn, pien.chrn, pial.chrn, laly.chrn) 
my_greys <- c('grey10', 'grey39', 'gray48', 'gray78')
my_greys <- c("#E2E2E2", "#717171", "#A5A2A2", "#4B4B4B")

theme_emma <- function(){  
  theme(
    axis.text = element_text(size =8),
    text = element_text(size=10),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=8, face="bold"),
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )}

growth <- ggplot(plt.dat, aes(x = rownames(plt.dat), y=Chron, group=Species))+
  #geom_hline(yintercept = 0.85, color = '#A2A5A5', size = 0.3)+
  geom_hline(yintercept = 1, colour = 'dark grey', size=0.1)+
  geom_path(size = 0.4) +
  geom_smooth(method= 'loess', colour = 'red', size = 0.2, span = 0.4, se= FALSE)+
  scale_color_manual(values = "#717171")+
  scale_x_discrete(breaks = seq(1900,2015, by = 20))+
  facet_wrap(~Species, nrow = 1)+
  ylab('Regional residual chronology')+
  xlab('Year')+
  theme_few()+
  theme_emma2()
growth
#ggsave(plot=growth, "~/Desktop/RegChronGrowth.pdf", device = "pdf", width = 8.5, height = 3, units = c("in"), dpi= 600)  

x <- runCor(abla.chrn$ABLstd, Tave_reg$Tave_sm, n = 20)
plot(seq(1900, 2015, by = 1),x, type = 'l')


y <- runCor(pien.chrn$PIEstd, Tave_reg$Tave_sm, n = 20)
plot(seq(1900, 2015, by = 1),y, type = 'l')


z <- runCor(laly.chrn$LALstd, Tave_laly$Tave_sm, n = 20)
plot(seq(1900, 2015, by = 1),z, type = 'l')

g <- runCor(pial.chrn$PIAstd, Tave_pial$Tave_sm[29:112], n = 20)
plot(seq(1928, 2011, by = 1),g, type = 'l')


#####NEXT: Try dendroclim analysis w these regional chrons .... check old JUNE 2019 Code
#DENDROCLIM ----
Tave_reg<- read.csv("Tave_reg2020.csv")[2:6]
PPT_reg<- read.csv("PPT_reg2020.csv")[2:6]
reg_ave <- as.data.frame(cbind(Tave_reg, PPT_reg[2:5], Tave_reg[2:5]))
names <- names(reg_ave)[2:9]

Tave_laly<- read.csv("Tave_laly2020.csv")[2:6]
PPT_laly<- read.csv("PPT_laly2020.csv")[2:6]
laly_ave <- as.data.frame(cbind(Tave_laly, PPT_laly[2:5], Tave_laly[2:5]))

Tave_pial<- read.csv("Tave_pial2020.csv")[2:6]
PPT_pial<- read.csv("PPT_pial2020.csv")[2:6]
pial_ave <- as.data.frame(cbind(Tave_pial, PPT_pial[2:5], Tave_pial[2:5]))


#DENDROCLIM
ABLA.dcc <- dcc(abla.chrn, reg_ave, vnames=c('SEA'), method="corr", start = 1, end = 8) %>%
  mutate(Variable = names,
         Species = "ABLA")

PIEN.dcc <- dcc(pien.chrn, reg_ave, vnames=c('SEA'), method="corr", start = 1, end = 8) %>%
  mutate(Variable = names,
         Species = "PIEN")

LALY.dcc <- dcc(laly.chrn, laly_ave, vnames=c('SEA'), method="corr", start = 1, end =8) %>%
  mutate(Variable = names,
         Species = "LALY")

PIAL.dcc <- dcc(pial.chrn, pial_ave, vnames=c('SEA'), method="corr", start = 1, end =8) %>%
  mutate(Variable = names,
         Species = "PIAL")

names(abla.chrn) = c("Chron", "nSamples")
abla.chrn$Species = c("ABLA")
names(pien.chrn) = c("Chron", "nSamples")
pien.chrn$Species = c("PIEN")
names(laly.chrn) = c("Chron", "nSamples")
laly.chrn$Species = c("LALY")
names(pial.chrn) = c("Chron", "nSamples")
pial.chrn$Species = c("PIAL")

dcdat <- rbind(ABLA.dcc, LALY.dcc, PIAL.dcc, PIEN.dcc) %>%
  mutate(Significant = significant) %>%
  mutate(Variable = factor(Variable, levels = c("Tave_wt","Tave_sp", "Tave_sm","Tave_at",
                                                "PPT_wt","PPT_sp","PPT_sm","PPT_at")))

a<-rep(c('bold'), 4) #Bolds for differentiating vars
b<-rep(c('plain'), 4)

theme_emma2<- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1, face=c(a,b)),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size = 12, face = "bold")
  )}

mypal <- c("#A2A5A5", "#4B4B4B")
mypal <- c("#E1E1E1", "#717171")


sig.plot <-ggplot(data=dcdat, aes(x=Variable, y=coef, fill = Significant))+
  geom_bar(stat="identity", width=0.5, color = "black", size = 0.2)+
  geom_errorbar(aes(ymin = ci.lower, ymax=ci.upper), size = 0.2, width = 0, color = 'black')+
  facet_wrap(~Species, nrow = 1)+
  scale_fill_manual(values = mypal)+
  theme_bw()+
  theme_emma()+
  ylim(-0.7, +0.7)+
  geom_hline(yintercept = 0, size=0.2)+
  theme(legend.position = "none")
sig.plot  

#ggsave(plot=sig.plot, "DendroClim2.pdf", device = "pdf", width = 8.5, height = 3, units = c("in"), dpi= 600)

#9 x 3







