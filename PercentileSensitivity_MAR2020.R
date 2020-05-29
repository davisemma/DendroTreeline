#THIS CODE: Generates analysis of treeline advance (15th percentile establishment with elevation) 
#(MAR 23/20)-------------------------------------------

rm(list=ls())
library('car')
library('dplyr')
library('ggplot2')
library('ggthemes')
library('ggsci')

setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Regional Treering /Data")

samples <- read.csv("REGION_SamplesJAN2019.csv") %>%
  select(-"X") %>%
  rename(Elev = Elev.y) %>%
  rename(Zone = Elev.x)

breaks<- seq(1590, 2015, by=5) #Create breaks
bin<-as.array((cut(samples$EstYear, breaks, right=FALSE))) #give each est year a bin
bins<-as.numeric(substring(bin,2,5)) #Record the bin value

#Combine, set levels for graphing
samples<-(cbind(samples,bins)) %>%
  mutate(SiteID = factor(SiteID, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM')),
         Maturity = factor(Maturity, levels = c('SEE', 'SAP', 'MAT', 'KRU')),
         Species = factor(samples$Species, levels=c('ABLA',  'PIEN',  'LALY',  'ASPEN', 'PIAL',  'MISL')))

#Colour pallette
mypal <- c("#E64B35B2" ,"#4DBBD5B2", "#00A087B2","#3C5488B2" ,"#F39B7FB2", "#8491B4B2","#91D1C2B2", "#DEBD7EB2", "#7E6148B2")

##CONSIDER:
mypal2 = pal_jco(alpha = 0.9)(9)
pal_jco
mypal2
show_col(mypal2)

show_col(mypal)
#Define labels
xlab <- c("\nEstablishment Year (5 yr bins)")
ylab <- c("\nProportional Tree Establishment")
matbreaks <- c('SEE','SAP', 'MAT')
sitebreaks <- c('WIL', 'HIL','GSP', 'FTB','FTA', 'HWA', 'SSG','PMR', 'HUM')
spbreaks <- c('ABLA',  'PIEN',  'LALY',  'ASPEN', 'PIAL',  'MISL')
yrbreaks <- seq(1580, 2020, by=20)


#SUMMARY
data<- samples %>%
  group_by(SiteID, Maturity) %>%
  summarise(nSamples = dplyr::n())

#15% establishment quantiles
theme_emma <- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=10),
    plot.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )}

#Standardizes elevations based on sites so they can be compared directly
std.samples <- samples %>%
  group_by(SiteID) %>%
  mutate(Std.Elev = as.numeric(scale(Elev)))

std.quantiles<- std.samples %>%
  group_by(SiteID, Std.Elev) %>%
  summarise(Perc = quantile(EstYear, c(.15), na.rm=TRUE)) %>%
  mutate(Std.Elev = as.numeric(Std.Elev, levels=Std.Elev))

quantiles<- samples %>%
  group_by(SiteID, Elev) %>%
  summarise(Perc = quantile(EstYear, c(.15), na.rm=TRUE)) %>%
  mutate(Elev = as.numeric(Elev, levels=Elev))

quant.est <- ggplot(std.quantiles,
                    aes(x=Perc, y=Std.Elev, group=1, color = SiteID))+
  # geom_path(color = "grey")+
  geom_smooth(method=lm, color='grey', se=FALSE, size = .5)+
  geom_point(size = 2)+
  scale_color_manual(values = mypal2)+
  facet_wrap(~ SiteID, nrow=3, ncol=3, scales="free_x")+ 
  xlab('\n Establishment year of 15th percentile')+
  ylab('Standardized elevation\n')+
  theme_bw()+
  theme(legend.position="none")

quant.est
#(Former green shade #456934)

#ggsave(plot=quant.est, "15thPerc_2019.pdf", device = "pdf", width = 6.75, height = 5, units = c("in"), dpi= 600)

#Percent established after 1900, after 1950
estperc<- na.omit(samples) %>%
  mutate(Old = ifelse(EstYear <= 1900, 1, 0),
         Young = ifelse(EstYear > 1900, 1, 0),
         MidCent = ifelse(EstYear > 1950, 1,0),
         Samples = 1)%>%
  group_by(SiteID)%>%
  summarise(Pre1900 = sum(Old),
            Post1900 = sum(Young),
            Post1950 = sum(MidCent),
            N = sum(Samples),
            Perc1900 = Post1900/N,
            Perc1950 = Post1950/N)

Av1900 <- mean(estperc$Perc1900)
SD1900 <- sd(estperc$Perc1900)
Av1950 <- mean(estperc$Perc1950)
SD1950 <- sd(estperc$Perc1950)


#15% quantiles, all sites
TLdist <- samples %>%
  group_by(SiteID, Zone) %>%
  summarise(Perc = quantile(EstYear, c(.15), na.rm=TRUE)) %>%
  mutate(Zone = as.numeric(Zone, levels=Zone))

all.site <- ggplot(TLdist, 
       aes(x = Perc, y = Zone,
           colour = SiteID, group=1))+
  geom_smooth(method = lm, col = 'grey', linetype = 2, se = FALSE, size = 0.5)+
  geom_smooth(aes(group=SiteID), method=lm, se = FALSE, size = 0.5)+
  geom_point(size = 2)+
  scale_colour_manual(values = mypal2, name = "Site ID")+
  xlab('Establishment year of oldest 15th percentile')+
  ylab('Distance from elevational midpoint (m)')+
  # xlim(1800, 2000)+
  ylim(-150, 150)+
  theme_bw()+
  guides(shape = FALSE, fill = FALSE, colour = FALSE)

all.site
#ggsave(plot=all.site, "AllSiteEstJCO_2019.pdf", device = "pdf", width = 4.5, height = 3.25, units = c("in"), dpi= 600)

#Building linear models
quantiles_test<- samples %>%
  group_by(SiteID, Elev) %>%
  summarise(Perc10 = quantile(EstYear, c(.10), na.rm=TRUE),
            Perc15 = quantile(EstYear, c(.15), na.rm=TRUE),
            Perc20 = quantile(EstYear, c(.20), na.rm=TRUE)) %>%
  mutate(Elev = as.numeric(Elev, levels=Elev),
         PercDiff10 = (abs(Perc10-Perc15)/Perc15)*100,
         PercDiff20 = (abs(Perc20-Perc15)/Perc15)*100)


quant_site15 = split(quantiles15, quantiles15$SiteID)

quant_site_map15 = map(quant_site15, function(x){
  lm(Elev ~ Perc15,
     data = x)
})
quant_site_map15

coefs_15 <- quant_site_map15 %>% map( ~ coefficients(.x))

#median
quantilesmed<- samples %>%
  group_by(SiteID, Elev) %>%
  summarise(Med = median(EstYear, na.rm = TRUE)) %>%
  mutate(Elev = as.numeric(Elev, levels=Elev))
quant_site_med = split(quantilesmed, quantilesmed$SiteID)

quant_site_map_med = map(quant_site_med, function(x){
  lm(Elev ~ Med,
     data = x)
})
quant_site_map_med

coefs_med <- quant_site_map_med %>% map( ~ coefficients(.x))

#25h percentile
quantiles25<- samples %>%
  group_by(SiteID, Elev) %>%
  summarise(Perc25 = quantile(EstYear, c(.25), na.rm=TRUE)) %>%
  mutate(Elev = as.numeric(Elev, levels=Elev))
quant_site25 = split(quantiles25, quantiles25$SiteID)

test <- as.data.frame(cbind(quantiles$SiteID, quantiles$Perc, quantiles25$Perc25, quantiles10$Perc10))




quant_site_map25 = map(quant_site25, function(x){
  lm(Elev ~ Perc25,
     data = x)
})
quant_site_map25

coefs_25 <- quant_site_map25 %>% map( ~ coefficients(.x))


test <- cbind(as.data.frame(unlist(coefs_med)),
             as.data.frame(unlist(coefs_15)),
             as.data.frame(unlist(coefs_25))) %>%
  mutate(Vars = row.names(.))
names(test) <- c('Median', 'Perc15', 'Perc25', 'Vars')

test2 <- test %>%
  mutate(PercDiff_Median = (Median - Perc15)/Perc15*100,
         PercDiff_25 = (Perc25 - Perc15)/Perc15*100)

test3 <- test2[c(2,4,6,8,10,12,14,16,18),]


modWIL<-lm(Elev ~ Perc, SiteID=='WIL', data=quantiles)
summary(modWIL)

modHIL<-lm(Elev ~ Perc, SiteID=='HIL', data=quantiles)
summary(modHIL)

modGSP<-lm(Elev ~ Perc, SiteID=='GSP', data=quantiles)
summary(modGSP)

modFTA<-lm(Elev ~ Perc, SiteID=='FTA', data=quantiles)
summary(modFTA)

modFTB<-lm(Elev ~ Perc, SiteID=='FTB', data=quantiles)
summary(modFTB)

modHUM<-lm(Elev ~ Perc, SiteID=='HUM', data=quantiles)
summary(modHUM)

modHWA<-lm(Elev ~ Perc, SiteID=='HWA', data=quantiles)
summary(modHWA)

modPMR<-lm(Elev ~ Perc, SiteID=='PMR', data=quantiles)
summary(modPMR)

modSSG<-lm(Elev ~ Perc, SiteID=='SSG', data=quantiles)
summary(modSSG)

#R2 and p-value from linear regression
mod.sum <- rbind(as.data.frame(cbind(Site = c('WIL'), r2 = summary(modWIL)$r.squared, p = summary(modWIL)$coefficients[2,4])),as.data.frame(cbind(Site = c('HIL'), r2 = summary(modHIL)$r.squared, p = summary(modHIL)$coefficients[2,4])),as.data.frame(cbind(Site = c('GSP'), r2 = summary(modGSP)$r.squared, p = summary(modGSP)$coefficients[2,4])),as.data.frame(cbind(Site = c('FTA'), r2 = summary(modFTA)$r.squared, p = summary(modFTA)$coefficients[2,4])),as.data.frame(cbind(Site = c('FTB'), r2 = summary(modFTB)$r.squared, p = summary(modFTB)$coefficients[2,4])),as.data.frame(cbind(Site = c('HUM'), r2 = summary(modHUM)$r.squared, p = summary(modHUM)$coefficients[2,4])),as.data.frame(cbind(Site = c('HWA'), r2 = summary(modHWA)$r.squared, p = summary(modHWA)$coefficients[2,4])),as.data.frame(cbind(Site = c('PMR'), r2 = summary(modPMR)$r.squared, p = summary(modPMR)$coefficients[2,4])),as.data.frame(cbind(Site = c('SSG'), r2 = summary(modSSG)$r.squared, p = summary(modSSG)$coefficients[2,4])))

#write.csv(mod.sum, "EstRegSumJAN2019.csv")

#Movement rates from slope of regression line
rates<-as.data.frame(rbind('WIL'=coef(modWIL), 'HIL'=coef(modHIL), 'GSP'=coef(modGSP),'FTA'=coef(modFTA),'FTB'=coef(modFTB),
                           'HUM'=coef(modHUM),'HWA'=coef(modHWA),
                           'PMR'=coef(modPMR),'SSG'=coef(modSSG)))
colnames(rates)<-c('Intercept', 'Slope')

#write.csv(rates, 'advratesJAN2019.csv')

#t-test if oldest 15th percentile of ALL TREES is older below than above treeline
tdat<-samples %>%
  mutate(Pos = ifelse(Zone<0, 'B','A')) %>%
  group_by(SiteID, Pos) %>%
  summarise(Perc = quantile(EstYear, c(.15), na.rm=TRUE)) %>%
  dcast(SiteID~Pos)

t.test(tdat$A, tdat$B, alternative=c("greater"), paired=TRUE)
mean(tdat$A)
mean(tdat$B)

tlong <-melt(tdat, id = "SiteID") %>%
  mutate(SiteID = factor(SiteID, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM')))

theme_emma <- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")
  )}

#Supple figure showing establishment years: Above below treeline
est.yrs <- ggplot(tlong, 
                  aes(x = SiteID, y = value,
                      fill = SiteID, shape = variable, alpha = variable))+
  geom_point(size=4)+
  scale_shape_manual(values = c(24, 25), name = "Position relative to midpoint", labels=c("Above", "Below"))+
  theme_few()+
  theme_emma()+
  scale_alpha_discrete(range = c(.5, 1))+
  scale_fill_manual(values = mypal2)+
  ylab("Est. date of oldest 15th percentile of trees")+
  xlab("Site")+
  ylim(1800, 2000)+
  guides(fill = FALSE, alpha = FALSE, shape = FALSE)+
  theme(legend.position="bottom")

est.yrs

#ggsave(plot=est.yrs, "~/Desktop/EstYrs.pdf", device = "pdf", width = 5, height = 4.5, units = c("in"), dpi= 600)

