#THIS CODE: Creates figure and summary data for changes in climate conditions
#Creates Mean summer temperature figure
rm(list=ls())
library("car")
library('dplyr')
library('ggplot2')
library('RColorBrewer')

source("http://peterhaschke.com/Code/multiplot.R")

data <- read.csv("Region2_Seasonal.csv")

mypal <-c("#428EAD","#FF875F","#22A782","#997965","#E0C06A",'#D19762', "#C48758", "grey")

mypal <-c("blue", "blue4", "gray20")

mypal <- c('darkolivegreen', 'peru', 'steelblue4', 'steelblue4', 'grey10','indianred3')

breaks<-c("Historical", "CNRM_RCP45", "CNRM_RCP85")

#SUMMER TEMPERATURE FIGURE
SMfig<-ggplot(filter(data, Data == 'Historical'),
              aes(x=Year, y=Tave_sm,
                  colour=Data))+
  geom_line(aes(linetype=Data, size=Data))+
  scale_linetype_manual(values=c("solid"),name="", breaks=breaks)+
  geom_smooth(method=loess, col=c('red'), size=0.4)+
  theme_few()+
  scale_colour_manual(values = "black", 
                      name="", 
                      breaks=breaks)+ 
  scale_size_manual(values=c(0.4,0.4,0.4),name="", breaks=breaks)+
  xlab('\nYear')+
  ylab('Regional mean\nsummer temperature (°C)\n')+
  xlim(1900, 2015)
SMfig
#(5 x 3.5)
###EXTRA ANALYSIS
data2<-data[-c(116:120, 206:210),]

AveSM<-ggplot(data2,
             aes(x=Year, y=Tave_sm,
                 colour=Data))+
  geom_line(aes(linetype=Data, size=Data))+
  scale_linetype_manual(values=c("dashed", "dashed", "solid"),name="", breaks=breaks)+
  geom_smooth(method=loess, col=c('indianred'), size=0.5)+
  theme_classic()+
  scale_colour_manual(values = c('dark grey', 'light grey', 'black'), 
                      name="", 
                      breaks=breaks)+ 
  scale_size_manual(values=c(0.4,0.4,0.4),name="", breaks=breaks)+
  xlab('\nYear')+
  ylab('Mean summer temperature (°C)\n')+
  theme(legend.position="bottom")
AveSM

PPTsm<-ggplot(data,
              aes(x=Year, y=PPT_sm,
                  colour=Data))+
  geom_line(aes(linetype=Data, size=Data))+
  scale_linetype_manual(values=c("solid", "solid", "solid"),name="", breaks=breaks)+
  geom_smooth(method=loess, col=c('indianred4'), size=0.5)+
  theme_classic()+
  scale_colour_manual(values = mypal, 
                      name="", 
                      breaks=breaks)+ 
  scale_size_manual(values=c(0.5,0.5,0.5),name="", breaks=breaks)+
  xlab('\nYear')+
  ylab('Summer Precipitation (mm)\n')+
  theme(legend.position="bottom")
PPTsm


PASwt<-ggplot(data,
              aes(x=Year, y=PAS_wt,
                  colour=Data))+
  geom_line(aes(linetype=Data, size=Data))+
  scale_linetype_manual(values=c("solid", "solid", "solid"),name="", breaks=breaks)+
  geom_smooth(method=loess, col=c('indianred4'), size=0.5)+
  theme_classic()+
  scale_colour_manual(values = mypal, 
                      name="", 
                      breaks=breaks)+ 
  scale_size_manual(values=c(0.5,0.5,0.5),name="", breaks=breaks)+
  xlab('\nYear')+
  ylab('Winter precipitation\n')+
  theme(legend.position="bottom")
PASwt


clim.sum <- data %>%
  group_by(Data) %>%
  summarize(AveSM = mean(Tave_sm),
            PASwt = mean(PAS_wt),
            PPTsm = mean(PPT_sm))







