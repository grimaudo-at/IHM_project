library(tidyverse)

#The parameters that Skylar estimated for deriving r. 
a <- 0.719
g <- 5.262
p <- 0.309
Tmax <- 21.5
Topt <- 14.05

#Equation for deriving r. This is called a Logan-10 function. 
r <- a*(((1+g*exp(-p*Temp))^-1) - exp(-(Tmax-Temp)/(Tmax-Topt)))


t.dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_raw_working.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  filter(behavior=="Torpor" & !is.na(datetime)) %>%
  group_by(site, trans_id, date) %>%
  summarise(temp = mean(temp))
#Reading in raw transmitter data and summarising the daily temperature of each, excluding arousal bouts. 


i.dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_25APR2023.csv") %>%
  filter(season=="hiber_earl") %>%
  mutate(gdL.c = mean_gdL + 0.0000000001, date=as.Date(date, format="%Y-%m-%d")) 
#This contains only the early hibernation data. I added an extremely small constant to make all gd values of 0 a positive number that can be 
#fed to the pathogen growth model. 


t.dat.fst <- t.dat %>%
  group_by(trans_id) %>%
  filter(date==min(date)) %>%
  mutate(gdL = i.dat$gdL.c[match(trans_id, i.dat$trans_id)]) %>%
  select(trans_id, date, gdL)
t.dat$gdL <- left_join(t.dat, t.dat.fst, by=c("trans_id", "date"), keep=NULL)
#Column on which to simulate pathogen growth. The first value for each bat needs to be its load value measured in early hibernation. 
