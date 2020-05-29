#THIS CODE: Produces species level chronologies by site, finds EPS range, calculates loess smoother (20yr), and plots chrons for each site
rm(list=ls())
pkgs <- c("dplR", "bootRes", "dplyr", "ggsci", "ggthemes", "scales")
lapply(pkgs, require, character.only=TRUE)

#Reading in text files
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

#Selecting only years 1900-2015
yrs<-as.character(seq(1900, 2015, by=1))
chr.yrs<- lapply(chr.fin, function(x) {
  subset(x, rownames(x) %in% yrs, select=c(res, samp.depth))
}) %>%
  lapply(., function(x) {
    na.omit(x)
  })

###EPS: 0.85 is accepted threshold; search for EPS among detrended series *as opposed to final chronology*
eps.test <- lapply(det.chron, rwi.stats.running, window.length=20) %>%
  lapply(., function(x) x[x$eps > 0.85, c("start.year")]) %>%
  lapply(., function(x) x[c(1)]) %>%
  as.data.frame(.) %>%
  gather(., Chrono, EPSyear, FTA_ABLA:WIL_PIEN, factor_key=TRUE)

#EPS for full time span
eps.1900 <- lapply(det.chron, rwi.stats.running, window.length=20) %>%
  lapply(., function(x) x[x$start.year > 1900,]) %>%
  lapply(., function(x) x[c(1,15)]) %>%
  do.call(rbind, . ) %>%
  mutate(Chron = row.names(.)) %>%
  mutate(Chron = substr(Chron, 1, 8))

###
chrons <- ls(chr.yrs, sorted=FALSE)
WIL.ABLA <-(chr.yrs$WIL_ABLA) %>%
  mutate(Chrono = chrons[2],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

WIL.PIEN <-  (chr.yrs$WIL_PIEN)%>%
  mutate(Chrono = chrons[1],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HIL.ABLA <-  (chr.yrs$HIL_ABLA)%>%
  mutate(Chrono = chrons[15],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HIL.PIEN <-  (chr.yrs$HIL_PIEN)%>%
  mutate(Chrono = chrons[14],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

GSP.LALY <-  (chr.yrs$GSP_LALY)%>%
  mutate(Chrono = chrons[17],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

GSP.ABLA <-  (chr.yrs$GSP_ABLA)%>%
  mutate(Chrono = chrons[18],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

GSP.PIEN <-  (chr.yrs$GSP_PIEN)%>%
  mutate(Chrono = chrons[16],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTA.ABLA <-  (chr.yrs$FTA_ABLA)%>%
  mutate(Chrono = chrons[24],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTA.PIEN <-  (chr.yrs$FTA_PIEN)%>%
  mutate(Chrono = chrons[22],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTA.LALY <-  (chr.yrs$FTA_LALY)%>%
  mutate(Chrono = chrons[23],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTB.ABLA <-  (chr.yrs$FTB_ABLA)%>%
  mutate(Chrono = chrons[21],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTB.LALY <-  (chr.yrs$FTB_LALY)%>%
  mutate(Chrono = chrons[20],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

FTB.PIEN <-  (chr.yrs$FTB_PIEN)%>%
  mutate(Chrono = chrons[19],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HUM.ABLA <-  (chr.yrs$HUM_ABLA)%>%
  mutate(Chrono = chrons[13],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HUM.PIEN <-  (chr.yrs$HUM_PIEN)%>%
  mutate(Chrono = chrons[12],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HWA.ABLA <-  (chr.yrs$HWA_ABLA)%>%
  mutate(Chrono = chrons[11],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HWA.PIEN <-  (chr.yrs$HWA_PIEN)%>%
  mutate(Chrono = chrons[9],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

HWA.PIAL <-  (chr.yrs$HWA_PIAL)%>%
  mutate(Chrono = chrons[10],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

PMR.ABLA <-  (chr.yrs$PMR_ABLA)%>%
  mutate(Chrono = chrons[8],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

PMR.PIAL <-  (chr.yrs$PMR_PIAL)%>%
  mutate(Chrono = chrons[7],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

PMR.PIEN <-  (chr.yrs$PMR_PIEN)%>%
  mutate(Chrono = chrons[6],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

SSG.ABLA <-  (chr.yrs$SSG_ABLA)%>%
  mutate(Chrono = chrons[5],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

SSG.PIAL <-  (chr.yrs$SSG_PIAL)%>%
  mutate(Chrono = chrons[4],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

SSG.PIEN <-  (chr.yrs$SSG_PIEN)%>%
  mutate(Chrono = chrons[3],
         Site = substr(Chrono, 1,3),
         Species = substr(Chrono, 5,8),
         Year = row.names(.),
         Loess = predict(loess(res ~ Year, span=20/nrow(.))))

plt.dat<-rbind(WIL.ABLA,WIL.PIEN,HIL.ABLA,HIL.PIEN,GSP.ABLA,GSP.LALY, GSP.PIEN,FTA.ABLA,FTA.LALY, FTA.PIEN,FTB.ABLA,FTB.LALY,FTB.PIEN,HUM.ABLA,HUM.PIEN,HWA.ABLA, HWA.PIAL, HWA.PIEN,PMR.ABLA,PMR.PIAL,PMR.PIEN,SSG.ABLA,SSG.PIAL, SSG.PIEN) %>%
  rename(., Chron = Chrono) %>%
  mutate(Chron = factor(Chron, levels = c("WIL_ABLA","WIL_PIEN","HIL_ABLA","HIL_PIEN","GSP_ABLA","GSP_LALY", "GSP_PIEN","FTA_ABLA","FTA_LALY", "FTA_PIEN","FTB_ABLA","FTB_LALY","FTB_PIEN","HUM_ABLA","HUM_PIEN","HWA_ABLA", "HWA_PIAL", "HWA_PIEN","PMR_ABLA","PMR_PIAL","PMR_PIEN","SSG_ABLA","SSG_PIAL", "SSG_PIEN")))

my_greys <- c('grey10', 'grey39', 'gray48', 'gray78')
my_greys <- c("#E2E2E2", "#717171", "#A5A2A2", "#4B4B4B")
eps.plt <- eps.1900 %>%
  mutate(Site = substr(Chron, 1, 3),
         Species = substr(Chron, 5, 8),
         Year = as.character(start.year)) %>%
  mutate(Chron = factor(Chron, levels = c("WIL_ABLA","WIL_PIEN","HIL_ABLA","HIL_PIEN","GSP_ABLA","GSP_LALY", "GSP_PIEN","FTA_ABLA","FTA_LALY", "FTA_PIEN","FTB_ABLA","FTB_LALY","FTB_PIEN","HUM_ABLA","HUM_PIEN","HWA_ABLA", "HWA_PIAL", "HWA_PIEN","PMR_ABLA","PMR_PIAL","PMR_PIEN","SSG_ABLA","SSG_PIAL", "SSG_PIEN"))) %>%
  filter(., eps < 3)

theme_emma2 <- function(){  
  theme(
    axis.text = element_text(size =10),
    text = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size = 12, face = "bold")
  )}

abla.plt <- ggplot(filter(plt.dat, Species == 'ABLA'), aes(x = Year, y=res, group=Chron))+
  #geom_area(data = eps.plt, 
  #         aes(x= Year, y = eps),
  #         fill = "#E2E2E2")+
  #geom_hline(yintercept = 0.85, color = '#A2A5A5', size = 0.3)+
  geom_hline(yintercept = 1, colour = 'grey19', size=0.3)+
  geom_path(size = 0.4) +
  geom_line(data=filter(plt.dat, Species == 'ABLA'),
            aes(x=Year, y=Loess),
            size=0.3, color = 'red')+
  scale_color_manual(values = "#717171")+
  scale_x_discrete(breaks = seq(1900,2000, by = 20))+
  scale_y_continuous(breaks = seq(0.5,1.5, by = 0.5))+
  coord_cartesian(ylim=c(0.5, 1.5))+
  facet_wrap(~Chron, ncol = 3)+
  ylab('Residual ring width')+
  theme_bw()+
  theme_emma2()

abla.plt

#ggsave(plot=abla.plt, "ABLA_Chrons_March2020.pdf", device = "pdf", width = 6.5, height = 5, units = c("in"), dpi= 600)

pien.plt <- ggplot(filter(plt.dat, Species == 'PIEN'), aes(x = Year, y=res, group=Chron))+
  #geom_area(data = eps.plt, 
  #         aes(x= Year, y = eps),
  #         fill = "#E2E2E2")+
  #geom_hline(yintercept = 0.85, color = '#A2A5A5', size = 0.3)+
  geom_hline(yintercept = 1, colour = 'grey19', size=0.3)+
  geom_path(size = 0.4) +
  geom_line(data=filter(plt.dat, Species == 'PIEN'),
            aes(x=Year, y=Loess),
            size=0.3, color = 'red')+
  scale_color_manual(values = "#717171")+
  scale_x_discrete(breaks = seq(1900,2000, by = 20))+
  scale_y_continuous(breaks = seq(0.5,1.5, by = 0.5))+
  coord_cartesian(ylim=c(0.5, 1.5))+
  facet_wrap(~Chron, ncol = 3)+
  ylab('Residual ring width')+
  theme_bw()+
  theme_emma2()

pien.plt

#ggsave(plot=pien.plt, "PIEN_Chrons_March2020.pdf", device = "pdf", width = 6.5, height = 5, units = c("in"), dpi= 600)

laly.plt <- ggplot(filter(plt.dat, Species == 'LALY'), aes(x = Year, y=res, group=Chron))+
  #geom_area(data = eps.plt, 
  #         aes(x= Year, y = eps),
  #         fill = "#E2E2E2")+
  #geom_hline(yintercept = 0.85, color = '#A2A5A5', size = 0.3)+
  geom_hline(yintercept = 1, colour = 'grey19', size=0.3)+
  geom_path(size = 0.4) +
  geom_line(data=filter(plt.dat, Species == 'LALY'),
            aes(x=Year, y=Loess),
            size=0.3, color = 'red')+
  scale_color_manual(values = "#717171")+
  scale_x_discrete(breaks = seq(1900,2000, by = 20))+
  scale_y_continuous(breaks = seq(0.5,1.5, by = 0.5))+
  coord_cartesian(ylim=c(0.5, 1.5))+
  facet_wrap(~Chron, ncol = 3)+
  ylab('Residual ring width')+
  theme_bw()+
  theme_emma2()

laly.plt

#ggsave(plot=laly.plt, "LALY_Chrons_March2020.pdf", device = "pdf", width = 6.5, height = 2.25, units = c("in"), dpi= 600)

pial.plt <- ggplot(filter(plt.dat, Species == 'PIAL'), aes(x = Year, y=res, group=Chron))+
  #geom_area(data = eps.plt, 
  #         aes(x= Year, y = eps),
  #         fill = "#E2E2E2")+
  #geom_hline(yintercept = 0.85, color = '#A2A5A5', size = 0.3)+
  geom_hline(yintercept = 1, colour = 'grey19', size=0.3)+
  geom_path(size = 0.4) +
  geom_line(data=filter(plt.dat, Species == 'PIAL'),
            aes(x=Year, y=Loess),
            size=0.3, color = 'red')+
  scale_color_manual(values = "#717171")+
  scale_x_discrete(breaks = seq(1900,2000, by = 20))+
  scale_y_continuous(breaks = seq(0.5,1.5, by = 0.5))+
  coord_cartesian(ylim=c(0.5, 1.5))+
  facet_wrap(~Chron, ncol = 3)+
  ylab('Residual ring width')+
  theme_bw()+
  theme_emma2()

pial.plt

#ggsave(plot=pial.plt, "PIAL_Chrons_March2020.pdf", device = "pdf", width = 6.5, height = 2.25, units = c("in"), dpi= 600)






