library(tidyverse)

#### Simulation on transmitter dataset ####
t.dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  filter(behavior=="Torpor" & !is.na(datetime)) %>%
  group_by(site, trans_id, date) %>%
  summarise(temp = mean(temp))
#Reading in raw transmitter data and summarising the daily temperature of each, excluding arousal bouts. 


i.dat.all <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_5JUN2023.csv") %>%
  mutate(date = as.Date(date, format="%m/%d/%y"))
#This dataframe contains all infection data from transmitter bats, both early and late sampling events. 

i.dat <- i.dat.all %>% 
 # mutate(gdL.c = mean_gdL + 10^((40-22.04942)/-3.34789)) %>%
  filter(season=="hiber_earl")
i.dat$mean_gdL[i.dat$mean_gdL == 0] <- 10^((40-22.04942)/-3.34789)
#This contains only the early hibernation data. For all bats without detectable infection in early hibernation, the 40ct value was assigned. 

t.dat.fst <- t.dat %>%
  group_by(site, trans_id) %>%
  filter(date==min(date)) %>%
  mutate(gdL = i.dat$mean_gdL[match(trans_id, i.dat$trans_id)]) %>%
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
  r <- a*((1/(1+g*exp(-p*Temp))) - exp(-(Tmax-Temp)/(Tmax-Topt)))
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

#pos.inf <-filter(t.dat, trans_id %in% t.dat.fst$trans_id[t.dat.fst$gdL>10^((40-22.04942)/-3.34789)])
#Dataframe of only those bats with a detectable infection in early hibernation

i.dat.all$true_mean_gdL.c[!is.na(i.dat.all$mean_gdL)]<- i.dat.all$mean_gdL[!is.na(i.dat.all$mean_gdL)] + 10^((40-22.04942)/-3.34789)
dl2<-ggplot()+
  geom_line(aes(x=date, y=log10(gdL)), data=t.dat)+
  geom_point(aes(x=date, y=log10(true_mean_gdL.c)), color="blue", data=i.dat.all[i.dat.all$trans_id %in% t.dat$trans_id,]) +
  facet_wrap(~trans_id);dl2
#Plot of only those bats with detectable early infection. Blue points are their true infection values. 

#Let's just plot those individuals that have transmitter data for the whole study period:

t.dat.lst <- t.dat %>%
  group_by(site, trans_id) %>%
  filter(date==max(date)) %>%
  filter(date>"2022-03-01")
#This dataframe contains the ending date and simulated pathogen load of those bats with transmitter data avaialbe for the 
#entire study period. 

dl3<-ggplot()+
  geom_line(aes(x=date, y=log10(gdL), color=site), data=t.dat[t.dat$trans_id %in% unique(t.dat.lst$trans_id),])+
  geom_point(aes(x=date, y=log10(true_mean_gdL.c)), color="blue", data=i.dat.all[i.dat.all$trans_id %in% t.dat.lst$trans_id,]) +
  facet_wrap(~trans_id);dl3

late.loads <- t.dat %>%
  group_by(trans_id) %>%
  filter(date==max(date))
hist(log10(late.loads$gdL))
#Histogram of simulated late hibernation pathogen loads. 

i.dat.late <- filter(i.dat.all, season=="hiber_late")
hist(log10(i.dat.late$true_mean_gdL))


#### Simulation on early temperature point ####

t.dat2 <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  filter(behavior=="Torpor" & !is.na(datetime)) %>%
  group_by(site, trans_id, date) %>%
  summarise(temp = mean(temp))
#Reading in raw transmitter data and summarising the daily temperature of each, excluding arousal bouts. 
#This is just to get the appropriate length of each transmitter's dataframe. The temperature data will be 
#replaced with each bat's early hibernation point temperature:
t.dat2$temp <- i.dat$temp[match(t.dat2$trans_id, i.dat$trans_id)]

#Again, bring in the early infection data to the t.dat2 dataframe:
t.dat2 <- left_join(t.dat2, t.dat.fst, by=c("site"="site", "trans_id"="trans_id", "date"="date")) %>%
  filter(trans_id %in% unique(t.dat.fst$trans_id[!is.na(t.dat.fst$gdL)])) #This line removes from the dataframe all temp/gdL data from bats that did not have early hib load data

t.dat2$r <- rFUN(t.dat2$temp)
#Calculating r values for each day-temperature

for(i in 2:nrow(t.dat2)) {if(is.na(t.dat2[i,5])==FALSE) {t.dat2[i,5] <- t.dat2[i,5]}
  else{t.dat2[i,5] <- t.dat2[i-1,5]*((10^(t.dat2[i-1,6]))^(1/7))}}
#Running the simulation

#pos.inf2 <-filter(t.dat2, trans_id %in% t.dat.fst$trans_id[t.dat.fst$gdL>10^((40-22.04942)/-3.34789)])
#Dataframe of only those bats with a detectable infection in early hibernation

late.loads2 <- t.dat2 %>%
  group_by(trans_id) %>%
  filter(date==max(date))
#Ending hibernation loads


#### Simulation on mean torpor bout temps ####

t.dat3 <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
  filter(behavior=="Torpor" & !is.na(datetime)) %>%
  group_by(site, trans_id, date) %>%
  summarise(temp = mean(temp))
#Reading in raw transmitter data and summarising the daily temperature of each, excluding arousal bouts. 
#This is just to get the appropriate length of each transmitter's dataframe. The temperature data will be 
#replaced with each bat's mean torpor bout temperature:

torpor.temps <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/arousals_torpors_working.csv")

t.dat3$temp <- torpor.temps$mean.torpor.temp[match(t.dat3$trans_id, torpor.temps$trans_id)]

#Again, bring in the early infection data to the t.dat2 dataframe:
t.dat3 <- left_join(t.dat3, t.dat.fst, by=c("site"="site", "trans_id"="trans_id", "date"="date")) %>%
  filter(trans_id %in% unique(t.dat.fst$trans_id[!is.na(t.dat.fst$gdL)])) #This line removes from the dataframe all temp/gdL data from bats that did not have early hib load data

t.dat3$r <- rFUN(t.dat3$temp)
#Calculating r values for each day-temperature

for(i in 2:nrow(t.dat3)) {if(is.na(t.dat3[i,5])==FALSE) {t.dat3[i,5] <- t.dat3[i,5]}
  else{t.dat3[i,5] <- t.dat3[i-1,5]*((10^(t.dat3[i-1,6]))^(1/7))}}
#Running the simulation

#pos.inf2 <-filter(t.dat2, trans_id %in% t.dat.fst$trans_id[t.dat.fst$gdL>10^((40-22.04942)/-3.34789)])
#Dataframe of only those bats with a detectable infection in early hibernation

late.loads3 <- t.dat3 %>%
  group_by(trans_id) %>%
  filter(date==max(date))
#Ending hibernation loads


#### Plotting simulation output #### 
#Density plot of simulation output and qpcr results:

#First, making data source columns:

late.loads$data = "Model Output\n(Overwinter Microclimate Data)"
late.loads2$data = "Model Output\n(Early Hibernation Data)"
late.loads3$data = "Model Output\n(Mean Torpor Temp)"
i.dat.late$data = "qPCR Data"

late.loads.red <- late.loads %>% ungroup() %>% select(data, site, trans_id, gdL)
late.loads2.red <- late.loads2 %>% ungroup() %>% select(data, site, trans_id, gdL)
late.loads3.red <- late.loads3 %>% ungroup() %>% select(data, site, trans_id, gdL)
i.dat.late2 <- i.dat.late %>% filter(trans_id %in% unique(late.loads$trans_id)) %>% select(data, site, trans_id, true_mean_gdL) 
colnames(i.dat.late2) <- c("data", "site","trans_id", "gdL")

load.compar = rbind(late.loads.red, late.loads2.red, i.dat.late2)

load.compar<-filter(load.compar, site!="MEAD MINE")

load.compar$lgdL <- log10(load.compar$gdL)

load.compar.means <- load.compar %>%
  group_by(data) %>%
  summarise(mean.load = mean(lgdL, na.rm=T))

load.compar$data<-factor(load.compar$data, levels = c("qPCR Data", "Model Output\n(Early Hibernation Data)", "Model Output\n(Overwinter Microclimate Data)"))
load.compar.means$data<-factor(load.compar.means$data, levels = c("qPCR Data", "Model Output\n(Early Hibernation Data)", "Model Output\n(Overwinter Microclimate Data)"))

g.load.compar = ggplot(aes(lgdL, fill=data), data=load.compar) +
  geom_density(alpha=0.8)+
  geom_vline(xintercept=load.compar.means$mean.load, linetype="dashed", color="black", linewidth=1)+
  scale_fill_manual(values=c("firebrick1","deepskyblue4","blueviolet"))+
  ylab("Density")+
  xlab(expression("Late Hibernation Pathogen Load (log-Nanograms of Fungal DNA)"))+
  scale_x_continuous(limits=c(-6,4))+
  scale_y_continuous(limits=c(0,0.5))+
  theme(
    panel.grid.major.y = element_line(color="grey89"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill="white", color="black", linewidth=1),
    legend.key.size = unit(0.9,'cm'),
    legend.background = element_rect(fill="white", color="black"),
    legend.key = element_blank(),
    legend.position = c(0.18,0.85),
    legend.spacing.y = unit(0.3,'cm'),
    legend.title = element_blank(),
    legend.text = element_text(family="Arial", color="black", size=12),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(family="Arial", size=15, color="black"),
    axis.text.x = element_text(family="Arial", size=15, color="black"),
    axis.title = element_text(family="Arial", size=20, color="black"),
    strip.text.x = element_text(family="Arial", size=15, color="black"),
    strip.background = element_rect(color=NULL, fill=NULL)
  )+
  guides(fill=guide_legend(byrow = T));g.load.compar
#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/For ESA 2023/late_lgdL~qpcr_simulation_3.PNG",g.load.compar,scale=3,width=8,height=6,units="cm",dpi=600)


