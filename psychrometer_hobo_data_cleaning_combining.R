library(tidyverse)

#This script is quick and sloppy and unreliable. 
dat <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/Logger Data/ENV_DAT_MASTER.csv")
dat$date<-as.Date(dat$date, format="%Y-%m-%d")

dat <- dat %>%
  filter(site == 'MAIN GRAPHITE MINE' | site=="MEAD MINE" | site=="CP TUNNEL" | site=="BLACKBALL" | site=="ZIMMERMAN" | site=='SOUTH LAKE MINE' | site=="ELROY SPARTA") %>% 
  filter(date >= "2021-10-01" & date <= "2022-05-01")

dat.graphite<-read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/Logger Data/Raw Data/NY 2022 apr/trimmed csv/Graphite C trimmed.csv")
dat.graphite$state<-"NY"
dat.graphite<-separate(dat.graphite, datetime, into=c('date','time'), sep=" ")
dat.graphite$date <- as.Date(dat.graphite$date, format="%m/%d/%y")
dat.graphite$vpd <- NA
dat<-rbind(dat, dat.graphite)


dat$section[dat$section=="B START"]<-"B"
dat$section[dat$section=="NA"]<-NA
dat$section[dat$section=="FALSE"]<-"F"

dat$logger_model[dat$logger_model=="psy-MX2303"]<-"MX2303"

meta<-read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv") %>%
  mutate(date_deployed = as.Date(date_deployed, format="%m/%d/%y"), date_retrieved=as.Date(date_retrieved, format="%m/%d/%y"))
meta$site[meta$site=="GRAPHITE"] <- 'GRAPHITE MINE'
meta$site[meta$site=="SOUTH LAKE"] <- 'SOUTH LAKE MINE'

dat$date_deployed <- meta$date_deployed[match(dat$site, meta$site)]
dat$date_retrieved <- meta$date_retrieved[match(dat$site, meta$site)]
dat <- filter(dat, date >= (date_deployed+1) & date <= (date_retrieved-1))

dat<-select(dat, "site","state","section","date","time","temp_1","temp_2","rh","vpd","logger_model")

#write.csv(dat, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/psychrometer_hobo_data_26JUN2023.csv", row.names=F)
