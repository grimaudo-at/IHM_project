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

i.dat<-read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_5JUN2023.csv")
sec.summ.bats <- i.dat %>% 
  mutate(c=1) %>% 
  group_by(site, season, section) %>% 
  summarise(num.bats = sum(c))

dat$section[dat$section=="B START"]<-"B"
dat$section[dat$section=="NA"]<-NA
sec.summ.env <- dat %>% 
  mutate(c=1) %>% 
  group_by(site, section) %>% 
  summarise(env.dat.avail = mean(c))
View(sec.summ.env)

sec.summ.bats <- pivot_wider(data=sec.summ.bats, id_cols = c(site, section), names_from = season, values_from = num.bats)

dat$logger_model[dat$logger_model=="psy-MX2303"]<-"MX2303"
#write.csv(dat, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/psychrometer_hobo_data_26JUN2023.csv", row.names=F)
