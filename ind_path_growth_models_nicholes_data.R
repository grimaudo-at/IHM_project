library(tidyverse)
library(lme4)


#### Fitting logistic growth models to culture data ####

culture_dat <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/nichole_model_fitting/nl_merge_culture_temp_1_6_22.csv") %>%
  filter(isolate_id %in% c("WH2020", 'PY2019', 'FS2020', 'BC2021','BM2019','BB2020','SC2020','NR2019','WL2020')) %>%
  mutate(experiment_day = as.integer(as.Date(date_sampled, format="%m/%d/%y") - as.Date(date_start, format="%m/%d/%y")))
#For now, filtering Nichole's raw culture data to just be contemporary isolates. 
#The experiment_day column is the day into the experiment that the sample was taken for qPCR. 

p.raw.contemporary <- ggplot()+
  geom_point(aes(x=experiment_day, y=lgdL, fill=temp), data=culture_dat, shape=21, size=4, alpha=0.8)+
  facet_wrap(~isolate_id);p.raw.contemporary



#### Fitting linear thermal performance curve through sub-11C temps ####

## First pass: simple line through Kate's r estimates
#Kate's fitted r values:
kates_r <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/fitted_r_on_nicholes_data_KATEJUL2023.csv")

#Column with unique strain name: 
kates_r$strain <- paste(kates_r$site, kates_r$collection.year, sep="_")

#Only going to use contemporary strains:
kates_r_contemporary <- filter(kates_r, collection.year>=2019)

#Plot of strain-specific r estimates:
p.r.ests1 <- ggplot(aes(x=temp, y=r, fill=strain, group=strain), data=kates_r_contemporary)+
  geom_point(shape=21, size=4, alpha=0.8);p.r.ests1
#Facet_wrapped by strain:
p.r.ests1.f <- ggplot(aes(x=temp, y=r, fill=strain, group=strain), data=kates_r_contemporary)+
  geom_point(shape=21, size=4, alpha=0.8)+
  facet_wrap(~strain);p.r.ests1.f

#For now, we're just going to try defining the thermal performance of the pathogen as a linear regression
#through the r estimates between the temp range 1 and 11C:

kates_r_contemporary_sub11 <- filter(kates_r_contemporary, temp<=11)
#contains only r's from contemporary strains sub-11 C.

m.linear <- lm(r~temp, data=kates_r_contemporary_sub11);summary(m.linear)
#Linear regression through r estimates. Parameter estimates below:
#Intercept = 0.202363
#Temp Beta = 0.015758

m.linear.temps <- data.frame(temp=seq(0,11,0.1))
m.linear.yhat <- as.data.frame(predict(m.linear, m.linear.temps, se.fit=T, type='response'))
m.linear.yhat.fin <- cbind(m.linear.temps, m.linear.yhat)
#The linear model's predictions

#Plot sub-11C data with linear model predictions:
p.r.ests1.sub11 <- ggplot()+
  geom_point(aes(x=temp, y=r, fill=strain), data=kates_r_contemporary_sub11, shape=21, size=4, alpha=0.8)+
  geom_line(aes(x=temp, y=fit), color="blue", data=m.linear.yhat.fin);p.r.ests1.sub11

#### Preparing iButton dataset for simulations ####
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
#This contains only the early hibernation data. 
#For all bats without detectable infection in early hibernation, the 40ct value was assigned. 
#If we haven't received swab sample data back yet, then this value is left as NA. 

t.dat.fst <- t.dat %>%
  group_by(site, trans_id) %>%
  filter(date==min(date)) %>%
  mutate(gdL = i.dat$mean_gdL[match(trans_id, i.dat$trans_id)]) %>%
  select(trans_id, date, gdL)
#This dataframe contains each individual's early hibernation load value, IF THEIR SWAB DATA WAS RETURNED BY NAU. 
#If they weren't positive for Pd in early hibernation, the 40ct value is used. 

t.dat <- left_join(t.dat, t.dat.fst, by=c("site"="site", "trans_id"="trans_id", "date"="date")) %>%
  filter(trans_id %in% unique(t.dat.fst$trans_id[!is.na(t.dat.fst$gdL)])) #This line removes from the dataframe all temp/gdL data from bats that did not have early hib load data
#This dataframe now contains each bat's temp data as well as their starting load value.
#The first value for each individual is its early hibernation load value, if available. 


## However, nearly every transmitter prematurely stopped recording data between two weeks and nearly a month and a half
# prior to late hibernation sampling. To fill in this missing data, I am going to assume their last recorded daily 
# torpor bout temperature was the temperature they remained at until the end of hibernation. We can tweak this assumption
# later. 

#First, need to construct a list of dataframes that mirrors t.dat but with empty temp and gdL columns and a date column
#that extends from the start to end of the sampling period. 

dfs <- list()
ids <- as.list(unique(t.dat$trans_id))
for (i in 1:length(unique(t.dat$trans_id))) {dfs[i] <- data.frame(date = seq((min(i.dat.all$date[i.dat.all$trans_id==unique(t.dat$trans_id)[i]])),
                                                                              max(i.dat.all$date[i.dat.all$trans_id==unique(t.dat$trans_id)[i]]),1))}
dfs2<-lapply(dfs, function(x) {as.data.frame(x)})
for( i in seq_along(dfs2)){
  dfs2[[i]]$trans_id <- rep(ids[i],nrow(dfs2[[i]]))
  colnames(dfs2[[i]]) <- c('date','trans_id')
}




#### Simulation using linear thermal performance on sub-11 temperatures ####


