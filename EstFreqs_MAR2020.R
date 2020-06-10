#THIS CODE: Creates models and plots of establishment frequencies and residuals in 5yr bins
#Establishment frequency calculations 
#June 10, 2020
rm(list=ls())
library('car')
library('dplyr')
library('ggplot2')
library('RColorBrewer')
library('ggthemes')

setwd("~/DendroTreeline")

samples <- read.csv("REGION_Samples.csv")

#Set levels order for graphing
samples$SiteID <- factor(samples$Site, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM'))
samples$Maturity <- factor(samples$Maturity, levels = c('SEE', 'SAP', 'MAT', 'KRU'))
samples$Species<-factor(samples$Species, levels=c('ABLA',  'PIEN',  'LALY',  'ASPEN', 'PIAL',  'MISL'))
samples$Species<-as.factor(samples$Species)
samp.sel <- samples %>%
  filter(., Species !="MISL")%>%
  filter(., Species !="ASPEN") %>%
  mutate(Species = factor(Species, levels = c('ABLA', 'LALY', 'PIAL', 'PIEN'))) %>%
  dplyr::mutate(Vline = ifelse(SiteID == "WIL", 1960, ifelse(SiteID == "HIL", 1955, ifelse(SiteID== "GSP", 1975, ifelse(SiteID == "FTA", 1965, ifelse(SiteID == "FTB", 1965, ifelse(SiteID == "PMR", 1955, ifelse(SiteID == "SSG", 1970, NA))))))))


#PLOT BY SITE
mypal <- c("#4B4B4B","#A2A5A5",'#E2E2E2', "#717171")
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

breaks<- seq(1590, 2015, by=5) #Create breaks
bin<-as.array((cut(samp.sel$EstYear, breaks, right=FALSE))) #give each est year a bin
bins<-as.numeric(substring(bin,2,5)) #Record the bin value
samp.sel<-(cbind(samp.sel,bins)) #Combine

#write.csv(samp.sel, 'REGIONAL_Samples_FIN.csv')


site.plt <- ggplot(samp.sel,
       aes(bins))+
  geom_histogram(binwidth=5, aes(fill = Species), color = "black", size = 0.1)+
  facet_wrap(~ SiteID, nrow=3, ncol=3, scales = 'free')+
  scale_fill_manual(values=mypal)+
  xlab("Year") + 
  ylab("Establishment frequency")+
  scale_x_continuous(breaks = seq(1800, 2050, by = 50), limits=c(1795, 2020))+
  theme_few()+ #Changed from theme_bw()
  theme_emma()+
  theme(legend.position="bottom")+
  geom_vline(aes(xintercept = Vline), linetype="longdash", col='red')

site.plt

#ggsave(plot=site.plt, "SpEstFreq_March2020b.pdf", device = "pdf", width = 6.5, height = 5, units = c("in"), dpi= 600)


#Set levels order for graphing
samples$SiteID <- factor(samples$Site, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM'))
samples$Maturity <- factor(samples$Maturity, levels = c('SEE', 'SAP', 'MAT', 'KRU'))
samples$Species<-factor(samples$Species, levels=c('ABLA',  'PIEN',  'LALY',  'ASPEN', 'PIAL',  'MISL'))
samples$Species<-as.factor(samples$Species)
samp.sel <- samples %>%
  filter(., Species !="MISL")%>%
  filter(., Species !="ASPEN") %>%
  mutate(Species = factor(Species, levels = c('ABLA', 'LALY', 'PIAL', 'PIEN'))) %>%
  dplyr::mutate(Vline = ifelse(SiteID == "WIL", 1960, ifelse(SiteID == "HIL", 1955, ifelse(SiteID== "GSP", 1975, ifelse(SiteID == "FTA", 1965, ifelse(SiteID == "FTB", 1965, ifelse(SiteID == "PMR", 1955, ifelse(SiteID == "SSG", 1970, NA))))))))

###PLOT BY SPECIES
mypal <- c("#4B4B4B","#A2A5A5",'#E2E2E2', "#717171")
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

breaks<- seq(1590, 2015, by=5) #Create breaks
bin<-as.array((cut(samp.sel$EstYear, breaks, right=FALSE))) #give each est year a bin
bins<-as.numeric(substring(bin,2,5)) #Record the bin value
samp.sel<-(cbind(samp.sel,bins)) #Combine

#Plot of establishment by species x site
sp.site <- ggplot(samp.sel,
                  aes(Elev))+
  geom_histogram(binwidth=5, aes(fill = Species), color = "black", size = 0.1)+
  facet_wrap(SiteID ~ Species, scales = 'free')+
  scale_fill_manual(values=mypal)+
  xlab("Year") + 
  ylab("Establishment frequency")+
  scale_x_continuous(breaks = seq(-130, 130, by = 10), limits=c(-130, 130))+
  theme_few()+ #Changed from theme_bw()
  theme_emma()+
  theme(legend.position="bottom")
  #geom_vline(aes(xintercept = Vline), linetype="longdash", col='red')

sp.site

