rm(list=ls())
library("car")
library('dplyr')
library('ggplot2')
library('reshape2')
library('psych')
library('ggthemes')

setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data")

#READ IN AND FORMAT SAMPLES from WIL and HIL
#Assuming all WIL and HIL dated samples are all samples
samples <- read.csv("REGION_Samples.csv") %>%
  filter(., Maturity != 'KRU', StudyArea != 'KC', StudyArea != 'GSP') %>%
  mutate(SiteID = factor(SiteID, levels = c('WIL','HIL')),
         Maturity = factor(Maturity, levels = c('SEE', 'SAP', 'MAT')))

#bins (1900 = 1900, 01, 02, 03, 04)
Bin<-as.numeric((substring(cut(samples$EstYear, seq(1590, 2015, by=5), right=FALSE),2,5)))

samples<-as.data.frame(cbind(samples,Bin))

mean.age<- samples %>%
  group_by(SiteID, Elev, Maturity) %>%
  summarise(Mean = mean(EstYear),
            Perc25 = quantile(EstYear, c(.25), na.rm = TRUE),
            Perc75 = quantile(EstYear, c(.75), na.rm = TRUE))

work.data<-merge(mean.age, samples, by=c('SiteID', 'Elev', 'Maturity'), all.y=TRUE)

MeanBin <- as.numeric((substring(cut(work.data$Mean, seq(1590, 2015, by=5), right=FALSE),2,5)))
Perc25Bin <-as.numeric((substring(cut(work.data$Perc25, seq(1590, 2015, by=5), right=FALSE),2,5)))
Perc75Bin <- as.numeric((substring(cut(work.data$Perc75, seq(1590, 2015, by=5), right=FALSE),2,5)))

work.data<-as.data.frame(cbind(work.data,MeanBin, Perc25Bin, Perc75Bin))

#Density calculations adjusted for lower seedling sample size
actual_ages<- work.data %>%
  mutate(AgeBin=Bin)%>% 
  group_by(SiteID,Maturity, AgeBin)%>%
  summarise(NumTrees=n())%>%
  mutate(CumSum=cumsum(NumTrees),
         Density=ifelse(SiteID == "WIL" & Maturity == "SEE", CumSum/0.585, ifelse(SiteID == "HIL" & Maturity == "SEE", CumSum/0.405, CumSum/0.825))) %>%
  group_by(SiteID, AgeBin) %>%
  summarise(NumTrees = sum(NumTrees),
            CumSum = sum(CumSum),
            Density = sum(Density)) %>%
  mutate(Max=max(Density),
         Year=AgeBin,
         Diff=Density-lag(Density, default=first(Density)),
         PropChange=Diff/Max) %>%
  filter(., Year > 1899)

mean_ages <- work.data %>%
  mutate(AgeBin=MeanBin)%>% 
  group_by(SiteID,Maturity, AgeBin)%>%
  summarise(NumTrees=n())%>%
  mutate(CumSum=cumsum(NumTrees),
         Density=ifelse(SiteID == "WIL" & Maturity == "SEE", CumSum/0.585, ifelse(SiteID == "HIL" & Maturity == "SEE", CumSum/0.405, CumSum/0.825))) %>%
  group_by(SiteID, AgeBin) %>%
  summarise(NumTrees = sum(NumTrees),
            CumSum = sum(CumSum),
            Density = sum(Density)) %>%
  mutate(Max=max(Density),
         Year=AgeBin,
         Diff=Density-lag(Density, default=first(Density)),
         PropChange=Diff/Max) %>%
  filter(., Year > 1899)

perc25_ages <- work.data %>%
  mutate(AgeBin=Perc25Bin)%>% 
  group_by(SiteID,Maturity, AgeBin)%>%
  summarise(NumTrees=n())%>%
  mutate(CumSum=cumsum(NumTrees),
         Density=ifelse(SiteID == "WIL" & Maturity == "SEE", CumSum/0.585, ifelse(SiteID == "HIL" & Maturity == "SEE", CumSum/0.405, CumSum/0.825))) %>%
  group_by(SiteID, AgeBin) %>%
  summarise(NumTrees = sum(NumTrees),
            CumSum = sum(CumSum),
            Density = sum(Density)) %>%
  mutate(Max=max(Density),
         Year=AgeBin,
         Diff=Density-lag(Density, default=first(Density)),
         PropChange=Diff/Max) %>%
  filter(., Year > 1899)

perc75_ages <- work.data %>%
  mutate(AgeBin=Perc75Bin)%>% 
  group_by(SiteID,Maturity, AgeBin)%>%
  summarise(NumTrees=n())%>%
  mutate(CumSum=cumsum(NumTrees),
         Density=ifelse(SiteID == "WIL" & Maturity == "SEE", CumSum/0.585, ifelse(SiteID == "HIL" & Maturity == "SEE", CumSum/0.405, CumSum/0.825))) %>%
  group_by(SiteID, AgeBin) %>%
  summarise(NumTrees = sum(NumTrees),
            CumSum = sum(CumSum),
            Density = sum(Density)) %>%
  mutate(Max=max(Density),
         Year=AgeBin,
         Diff=Density-lag(Density, default=first(Density)),
         PropChange=Diff/Max) %>%
  filter(., Year > 1899)


theme_emma2 <- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )}


all_plot<-ggplot()+
  geom_point(data = actual_ages,
             aes(x=Year, y=Density), size=1)+
  geom_path(data = actual_ages,
            aes(x=Year, y=Density),color = "#717171")+
  
  geom_point(data = mean_ages,
             aes(x=Year, y=Density), size=1)+
  geom_path(data = mean_ages,
            aes(x=Year, y=Density),color = "red")+
  
  geom_point(data = perc25_ages,
             aes(x=Year, y=Density), size=1)+
  geom_path(data = perc25_ages,
            aes(x=Year, y=Density),color = "blue")+
  
  geom_point(data = perc75_ages,
             aes(x=Year, y=Density), size=1)+
  geom_path(data = perc75_ages,
            aes(x=Year, y=Density),color = "green")+
  #breaks=elevbreaks)+
  facet_wrap(~SiteID, scales="free_y", ncol=3, nrow=3)+
  xlab('Year')+
  ylab('Tree density (trees/ha)\n')+
  scale_x_continuous(breaks = seq(1900, 2000, by = 20), limits=c(1900, 2000))+
  theme_bw()+
  theme_emma2()

all_plot

#BREAK POINT ANALYSIS
#INSTALL BREAKPOINT PACKAGE - see how breakpoint year differs between three methods
#ACTUAL
WILa <- as.data.frame(actual_ages[actual_ages$SiteID == 'WIL', ])
HILa <- as.data.frame(actual_ages[actual_ages$SiteID == 'HIL', ])

test_a <- list("WIL" = as.data.frame(WILa[,5]), "HIL" = as.data.frame(HILa[,5]))
test2_a <- lapply(test_a, CE.Normal.Mean, Nmax = 3, h = 7)

WILbp_a <- data.frame(BP.LOC = test2_a[[1]]$BP.Loc) %>%
  mutate(Site = 'WIL')
HILbp_a <- data.frame(BP.LOC = test2_a[[2]]$BP.Loc) %>%
  mutate(Site = 'HIL')

bp.dat_a <- rbind(WILbp_a, HILbp_a)

WIL_a <- as.data.frame(actual_ages[actual_ages$SiteID == 'WIL', ]) [c(13),]
HIL_a <- as.data.frame(actual_ages[actual_ages$SiteID == 'HIL', ]) [c(12),]
bp.yrs.a <- rbind(WIL_a, HIL_a)

#MEAN
WILmean <- as.data.frame(mean_ages[mean_ages$SiteID == 'WIL', ])
HILmean <- as.data.frame(mean_ages[mean_ages$SiteID == 'HIL', ])

test_mean<- list("WIL" = as.data.frame(WILmean[,5]), "HIL" = as.data.frame(HILmean[,5]))
test2_mean <- lapply(test_mean, CE.Normal.Mean, Nmax = 3, h = 7)

WILbp_mean <- data.frame(BP.LOC = test2_mean[[1]]$BP.Loc) %>%
  mutate(Site = 'WIL')
HILbp_mean <- data.frame(BP.LOC = test2_mean[[2]]$BP.Loc) %>%
  mutate(Site = 'HIL')

bp.dat_mean <- rbind(WILbp_mean, HILbp_mean)

WIL_mean <- as.data.frame(mean_ages[mean_ages$SiteID == 'WIL', ]) [c(11),]
HIL_mean <- as.data.frame(mean_ages[mean_ages$SiteID == 'HIL', ]) [c(11),]
bp.yrs.mean <- rbind(WIL_mean, HIL_mean)

#PERC25
WILperc75 <- as.data.frame(perc75_ages[perc75_ages$SiteID == 'WIL', ])
HILperc75 <- as.data.frame(perc75_ages[perc75_ages$SiteID == 'HIL', ])

test_perc75<- list("WIL" = as.data.frame(WILperc75[,5]), "HIL" = as.data.frame(HILperc75[,5]))
test2_perc75 <- lapply(test_perc75, CE.Normal.Mean, Nmax = 3, h = 7)

WILbp_perc75 <- data.frame(BP.LOC = test2_perc75[[1]]$BP.Loc) %>%
  mutate(Site = 'WIL')
HILbp_perc75 <- data.frame(BP.LOC = test2_perc75[[2]]$BP.Loc) %>%
  mutate(Site = 'HIL')

bp.dat_perc75 <- rbind(WILbp_perc75, HILbp_perc75)

WIL_perc75 <- as.data.frame(perc75_ages[perc75_ages$SiteID == 'WIL', ]) [c(11),]
HIL_perc75 <- as.data.frame(perc75_ages[perc75_ages$SiteID == 'HIL', ]) [c(11),]
bp.yrs.perc75 <- rbind(WIL_perc75, HIL_perc75)
