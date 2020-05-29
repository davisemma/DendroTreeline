#THIS CODE: Estimates changes in mean sumemr temperatures and compares to observed changes in treeline position
#MARCH 2020 ----------------------------------------
rm(list=ls())
setwd("~/Desktop/PhD/Research Data/Dendro Analysis/Climate WNA/Mean Growing Season/")

WIL<-read_csv('WIL_AveGS.csv') %>% .[1:100,] %>% mutate(SM_Ave = rowMeans(.[,3:5]))
HIL<-read_csv('HIL_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
GSP<-read_csv('GSP_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
FTA<-read_csv('FTA_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
FTB<-read_csv('FTB_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
HUM<-read_csv('HUM_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
HWA<-read_csv('HWA_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
PMR<-read_csv('PMR_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))
SSG<-read_csv('SSG_AveGS.csv') %>% .[1:100,]%>% mutate(SM_Ave = rowMeans(.[,3:5]))


#Mean SUMMER
SM_Ave<-as.data.frame(cbind(Year=c(FTA$Year) ,WIL=c(WIL$SM_Ave), HIL=c(HIL$SM_Ave), GSP=c(GSP$SM_Ave), FTA=c(FTA$SM_Ave), FTB=c(FTB$SM_Ave), HUM=c(HUM$SM_Ave), HWA=c(HWA$SM_Ave), PMR=c(PMR$SM_Ave), SSG=c(SSG$SM_Ave)))

WILlm <- lm(SM_Ave$WIL ~ SM_Ave$Year)
HILlm <- lm(SM_Ave$HIL ~ SM_Ave$Year)
GSPlm <- lm(SM_Ave$GSP ~ SM_Ave$Year)
FTAlm <- lm(SM_Ave$FTA ~ SM_Ave$Year)
FTBlm <- lm(SM_Ave$FTB ~ SM_Ave$Year)
HUMlm <- lm(SM_Ave$HUM ~ SM_Ave$Year)
HWAlm <- lm(SM_Ave$HWA ~ SM_Ave$Year)
PMRlm <- lm(SM_Ave$PMR ~ SM_Ave$Year)
SSGlm <- lm(SM_Ave$SSG ~ SM_Ave$Year)

#Added range of 4 to 8 degrees to show insensitivity ---
exp.obs<-as.data.frame(rbind(coefficients(WILlm),coefficients(HILlm),coefficients(GSPlm),coefficients(FTAlm),coefficients(FTBlm),coefficients(HWAlm),coefficients(SSGlm),coefficients(PMRlm),coefficients(HUMlm))) %>%
  mutate(Site = c("WIL", "HIL", "GSP", "FTA", "FTB", "HWA", "SSG","PMR", "HUM"),
         SiteID = factor(Site, levels = c('WIL', 'HIL', 'GSP', 'FTA','FTB', 'HWA', 'SSG','PMR', 'HUM'))) %>% 
  rename(., Intercept = `(Intercept)`, Slope = `SM_Ave$Year`) %>%
  mutate(MPY = (Slope*1000)/6,
         MPY_Hi = (Slope*1000)/4,
         MPY_Lo = (Slope*1000)/8, 
         EstMov = c(1.995245902, 0.308164764, 0.371605703, 0.676354241, 0.379501445,1.142358753, 1.83386617, 0.230821958, 0.537934365),
         Perc = (EstMov/MPY) *100,
         Var = "SM_Ave")

temp.sum<-exp.obs %>%
  mutate(Outcome = ifelse(Perc > 100, "Exceed", "Lag"),
         Mov100 = EstMov * 100,
         Clim100 = MPY * 100,
         Clim100_Hi = MPY_Hi * 100,
         Clim100_Lo = MPY_Lo * 100,)

#Mean and standard dev of mean summer T
x <- temp.sum$Slope * 100
mean(x)
sd(x)

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

mypal <- c("#E64B35B2" ,"#4DBBD5B2", "#00A087B2","#3C5488B2" ,"#F39B7FB2", "#8491B4B2","#91D1C2B2", "#DEBD7EB2", "#7E6148B2")

mypal2 = pal_jco(alpha = 0.9)(9)

adv.lag <- ggplot(temp.sum, aes(x = Clim100, y = Mov100, fill = SiteID, shape = Var))+
  geom_errorbarh(aes(xmin = Clim100_Lo, xmax = Clim100_Hi, color = SiteID))+
  scale_color_manual(values = mypal2)+
  geom_jitter(size = 4, width = 4)+
  geom_abline(intercept = 0, slope = 1, color = "#717171", linetype = 2)+
  xlim(0,300)+
  ylim(0,300)+
  scale_shape_manual(values = c(21))+
  scale_fill_manual(values = mypal2)+
  theme_few()+
  theme_emma()+
  ylab("Estimated 20th century advance (m)")+
  xlab("Expected 20th century advance (m)")+
  guides(shape = FALSE, fill = FALSE, color = FALSE)

adv.lag
#CC submission - 4 x 3.7
#ggsave(plot=adv.lag, "~/Desktop/AdvLag2020.pdf", device = "pdf", width = 4, height = 3.7, units = c("in"), dpi= 600)

#write.csv(dat.sum, 'SignifSum.csv')
#write.csv(temp.sum, 'ClimChSum2.csv')

###END HERE 
mean(temp.sum$EstMov)
sd(temp.sum$EstMov)
median(temp.sum$Mov100)
median(temp.sum$Clim100)
median(temp.sum$Perc)
mean(temp.sum$Perc)


sm.dat<- as.data.frame(cbind(WIL[,1], WIL[,7], HIL[,7], GSP[,7], FTA[,7], FTB[,7], HWA[,7], HUM[,7], PMR[,7], SSG[,7]))
sm.dat$Means <- rowMeans(gs.dat[,2:10])

plot(sm.dat[,1], sm.dat[,11], type = 'l')

