library(tidyverse)

t.dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_raw_working.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  filter(behavior=="Torpor" & !is.na(datetime)) %>%
  group_by(site, trans_id, date) %>%
  summarise(temp = mean(temp))
#Reading in raw transmitter data and summarising the daily temperature of each, excluding arousal bouts. 


i.dat.all <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_5JUN2023.csv")
  mutate(date = as.Date(date, format="%m/%d/%y"))
#This dataframe contains all infection data from transmitter bats, both early and late sampling events. 

i.dat <- i.dat.all %>% 
  mutate(gdL.c = mean_gdL + 0.0000000001) %>%
  filter(season=="hiber_earl")
#This contains only the early hibernation data. I added an extremely small constant to make all gd values of 0 a positive number that can be 
#fed to the pathogen growth model. r * 0 always = 0, but bats need to become infected at some point (simulate times to infection?). 

t.dat.fst <- t.dat %>%
  group_by(site, trans_id) %>%
  filter(date==min(date)) %>%
  mutate(gdL = i.dat$gdL.c[match(trans_id, i.dat$trans_id)]) %>%
  select(trans_id, date, gdL)
#This dataframe contains each individual's early hibernation load value, if available, with the small constant. 
t.dat <- left_join(t.dat, t.dat.fst, by=c("site"="site", "trans_id"="trans_id", "date"="date")) %>%
  filter(trans_id %in% unique(t.dat.fst$trans_id[!is.na(t.dat.fst$gdL)])) #This line removes from the dataframe all temp/gdL data from bats that did not have early hib load data
#This dataframe now contains each bat's temp data as well as their starting load value.
#The first value for each individual is its early hibernation load value, if available. 


#The parameters that Skylar estimated for deriving r. 
a <- 0.719
g <- 5.262
p <- 0.309
Tmax <- 21.5
Topt <- 14.05

#The below function requires a temperature input, Temp, and then returns the value of r from the temperature performance curve with 
#the above parameters estimated by Skylar. 
rFUN <- function(Temp) {
  r <- a*(1/((1+g*exp(-p*Temp))) - exp(-(Tmax-Temp)/(Tmax-Topt)))
  return(r)
}

h.temps <- as.vector(seq(0.01,21.5,0.01))
r.temps <- rFUN(h.temps)
r.crv <- as.data.frame(cbind(h.temps, r.temps))
r.crv.p <- ggplot(aes(x=h.temps, y=r.temps), data=r.crv)+
  geom_line();r.crv.p
#This is our temperature response curve, but it actually doesn't seem to match Skylar's curve exactly. For example,
#why isn't the greatest r value at 14.05, the Topt value? 

t.dat$r <- rFUN(t.dat$temp)

for(i in 2:nrow(t.dat)) {if(is.na(t.dat[i,5])==FALSE) {t.dat[i,5] <- t.dat[i,5]}
  else{t.dat[i,5] <- t.dat[i-1,5]*((10^(t.dat[i-1,6]))^(1/7))}}
#This function breaks the weekly pathogen growth rate at a given, temperature-dependent r (calculated by Skylar), into the 7th root,
#which represents a daily lambda. It then takes that daily lambda and multiplies it by the day's pathogen load value to 
#simulate the following day's pathogen load value. 

dl<-ggplot(aes(x=date, y=log10(gdL)), data=t.dat)+
  geom_line()+
  facet_wrap(~trans_id);dl
#Plot of all bats' pathogen growth. 

pos.inf <-filter(t.dat, trans_id %in% t.dat.fst$trans_id[t.dat.fst$gdL>0.0000000001])
#Dataframe of only those bats with a detectable infection in early hibernation

dl2<-ggplot()+
  geom_line(aes(x=date, y=log10(gdL)), data=pos.inf)+
  geom_point(aes(x=date, y=log10(true_mean_gdL)), color="blue", data=i.dat.all[i.dat.all$trans_id %in% pos.inf$trans_id,]) +
  facet_wrap(~trans_id);dl2
#Plot of only those bats with detectable early infection. Blue points are their true infection values. 

late.loads <- pos.inf %>%
  group_by(trans_id) %>%
  filter(date==max(date))
hist(log10(late.loads$gdL))
#Histogram of simulated late hibernation pathogen loads. 

i.dat.late <- filter(i.dat.all, season=="hiber_late")
hist(log10(i.dat.late$true_mean_gdL))
