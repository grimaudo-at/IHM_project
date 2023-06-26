library(tidyverse)
library(ggpubr)
#This script is for plotting all of our environmental data that was collected in IHM sites over the study period. 

#Let's start with barometer data:

bar.data <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/barometer_master.csv")
#Barometer master database. 

bar.data$datetime <- as.POSIXct(strptime(bar.data$datetime, "%Y-%m-%d %H:%M:%S",tz='EST')) 
bar.data$date <- as.Date(bar.data$date, format="%Y-%m-%d")
#Formating datetime and date correctly

hg.p <- ggplot(aes(x=datetime, y=pressure_hg, color=site), data=bar.data) +
  geom_line(size=0.6);hg.p
#Plot of barometric pressure over the study period in each site. 

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/barometric_pressure_sites.PNG",hg.p,scale=3,width=12,height=6,units="cm",dpi=600)

bar.temp.p <- ggplot(aes(x=datetime, y=temp, color=site), data=bar.data) +
  geom_line(size=0.6);bar.temp.p
#Plot of tempertaures recorded by barometer data loggers

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/barometer_temperatures_sites.PNG",bar.temp.p,scale=3,width=12,height=6,units="cm",dpi=600)

bar.rh.p <- ggplot(aes(x=datetime, y=rh, color=site), data=bar.data) +
  geom_line(size=0.6);bar.rh.p
#Plot of relative humidity values recorded by barometer data loggers

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/barometer_rh_sites.PNG",bar.rh.p,scale=3,width=12,height=6,units="cm",dpi=600)


#Moving on to psychrometer and HOBO data:

psych.hobo.dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/psychrometer_hobo_data_26JUN2023.csv")
#This dataframe contains data from both psychrometers and ProV2 loggers. I am going to make a 'uniq' column that
#contains a unique ID for each one, using its site, section, and logger_model metadata:

psych.hobo.dat$uniq.id <- paste(psych.hobo.dat$site, psych.hobo.dat$section, psych.hobo.dat$logger_model)
unique(psych.hobo.dat$uniq.id)
#There are 36 unique loggers in this dataset. 

#To create a datetime column, need to combine the date and time columns, which were separated:
psych.hobo.dat$datetime <- paste(psych.hobo.dat$date, psych.hobo.dat$time, sep=" ")
#Then, need to define the data type as datetime:
psych.hobo.dat$datetime <- as.POSIXct(strptime(psych.hobo.dat$datetime, "%Y-%m-%d %H:%M",tz='EST')) 

temp.psych.p <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat)+
  geom_line()+
  facet_wrap(~site);temp.psych.p
#Psychromter/hobo plot of temperatures
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/temp_psych_facet_sites.PNG",temp.psych.p,scale=3,width=12,height=6,units="cm",dpi=600)

#Also want to plot sites individually:
temp.psych.p.bball <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="BLACKBALL",])+
  ggtitle("BLACKBALL")+
  geom_line();temp.psych.p.bball
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/blackball_section_temps.PNG",temp.psych.p.bball,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.cptun <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="CP TUNNEL",])+
  ggtitle("CP TUNNEL")+
  geom_line();temp.psych.p.cptun
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/cptunnel_section_temps.PNG",temp.psych.p.cptun,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.elroy <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="ELROY SPARTA",])+
  ggtitle("ELROY SPARTA")+
  geom_line();temp.psych.p.elroy
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/elroy_section_temps.PNG",temp.psych.p.elroy,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.graphite <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="GRAPHITE MINE",])+
  ggtitle("GRAPHITE MINE")+
  geom_line();temp.psych.p.graphite
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/graphite_section_temps.PNG",temp.psych.p.graphite,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.mead <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="MEAD MINE",])+
  ggtitle("MEAD MINE")+
  scale_y_continuous(limits=c(5,10))+
  geom_line();temp.psych.p.mead
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/mead_section_temps.PNG",temp.psych.p.mead,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.southlake <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="SOUTH LAKE MINE",])+
  ggtitle("SOUTH LAKE MINE")+
  geom_line();temp.psych.p.southlake
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/southlake_section_temps.PNG",temp.psych.p.southlake,scale=3,width=6,height=6,units="cm",dpi=600)

temp.psych.p.zimm <- ggplot(aes(x=datetime, y=temp_1, color=section, linetype=logger_model, group=uniq.id), data=psych.hobo.dat[psych.hobo.dat$site=="ZIMMERMAN",])+
  ggtitle("ZIMMERMAN")+
  geom_line();temp.psych.p.zimm
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/zimmerman_section_temps.PNG",temp.psych.p.zimm,scale=3,width=6,height=6,units="cm",dpi=600)
