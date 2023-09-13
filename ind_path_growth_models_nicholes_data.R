library(tidyverse)
library(lme4)


#### Fitting logistic growth models to culture data ####

culture_dat <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/nichole_model_fitting/nl_merge_culture_temp_1_6_22.csv") %>%
  filter(isolate_id %in% c("WH2020", 'PY2019', 'FS2020', 'BC2021','BM2019','BB2020','SC2020','NR2019','WL2020')) %>%
  mutate(experiment_day = as.integer(as.Date(date_sampled, format="%m/%d/%y") - as.Date(date_start, format="%m/%d/%y")))
#For now, filtering Nichole's raw culture data to just be contemporary isolates. 
#The experiment_day column is the day into the experiment that the sample was taken for qPCR. 

p.raw.contemporary <- ggplot()+
  geom_point(aes(x=experiment_day, y=gdL, fill=temp), data=culture_dat, shape=21, size=4, alpha=0.8)+
  scale_y_continuous(limits=c(0,50))+
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
i.dat.all$mean_gdL[i.dat.all$mean_gdL == 0] <- 10^((40-22.04942)/-3.34789)
i.dat.all$lgdL <- log10(i.dat.all$mean_gdL)
#This dataframe contains all infection data from transmitter bats, both early and late sampling events. 
#For all bats without detectable infection in early hibernation, the 40ct value was assigned. 
#If we haven't received swab sample data back yet, then this value is left as NA. 

i.dat <- filter(i.dat.all, season=="hiber_earl")
#This contains only the early hibernation data. 
i.dat.late <- filter(i.dat.all, season=="hiber_late") #Infection data only from late hibernation
#Late hibernation infection data. 


t.dat.fst <- t.dat %>%
  group_by(site, trans_id) %>%
  filter(date==min(date)) %>%
  mutate(gdL = i.dat$mean_gdL[match(trans_id, i.dat$trans_id)]) %>% 
  dplyr::select(trans_id, date, gdL)
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
ids <- as.list(unique(t.dat$trans_id)) #List of all unique transmitter ID's
for (i in 1:length(unique(t.dat$trans_id))) {dfs[i] <- data.frame(date = seq((min(i.dat.all$date[i.dat.all$trans_id==unique(t.dat$trans_id)[i]])),
                                                                              max(i.dat.all$date[i.dat.all$trans_id==unique(t.dat$trans_id)[i]]),1))}
#^^^List of dataframes containing a sequence of dates for each transmitter ID, from early to late hibernation sampling dates. 
dfs2<-lapply(dfs, function(x) {as.data.frame(x)})
for( i in seq_along(dfs2)){
  dfs2[[i]]$trans_id <- rep(ids[[i]],nrow(dfs2[[i]]))
  colnames(dfs2[[i]]) <- c('date','trans_id')
}
#Matches the transmitter ID's in ids to date dataframes in dfs, matching on list index. 

dfs3 <- lapply(dfs2, function(x) {filter(x, date != min(date) & date != max(date))})
#Removes each dataframe's first and last date because sampling was done on those days and therefore, temperature/behavior data 
#could be inaccurate/misleading. 

#Can now merge all dataframes together into a single dataframe:
t.dat.empty <- bind_rows(dfs3)

t.dat.empty$trans_id <- as.factor(t.dat.empty$trans_id)
t.dat$trans_id <- as.factor(t.dat$trans_id)
#Now, on this t.dat.empty dataframe, I am going to join in the data from t.dat:
t.dat.full <- left_join(t.dat.empty, t.dat, by=c("trans_id","date"))

#For most transmitter datasets, there's a bunch of empty data on the end of the hibernation period when the logger stopped recording.
#I need to fill in these missing temp data WITH THE LAST RECORDED TEMP VALUE. To do so, I'm just going to write a for-loop that replaces NA values 
#with the value from the above cell, as long as the current cell is NA:

for(i in 2:nrow(t.dat.full)) {if(is.na(t.dat.full[i,4]) == 'TRUE') {t.dat.full[i,4] <- t.dat.full[i-1,4]}} #Filling in missing temp data
for(i in 2:nrow(t.dat.full)) {if(is.na(t.dat.full[i,3]) == 'TRUE') {t.dat.full[i,3] <- t.dat.full[i-1,3]}} #Filling in missing site data


#t.dat.full is now ready for the model. 



#### Simulation using linear thermal performance on sub-11 temperatures ####

#The pathogen growth model using linear terms for both temperature-dependent pathogen growth and decay will take the form:
# dP/dt = rP - uP
# Where r = mT + b ; u = cT + d ; T = average daily iButton temperature
# -> dP/dt = (r - u)P
# -> dP/dt = ((mT + b) - (cT + d))P
# Thus, P(t) = P(0)exp(((mT + b) - (cT + d))*t)
# Or, P(t+1) = P(0)exp((mT + b) - (cT + d))

# The values of m and b were estimated from our simple linear model earlier in the script:
b = 0.202363
m = 0.015758

#A first guess at what the values of c and d could be:
d = 0.1
c = 0.033


tp.lines <- ggplot()+
  geom_abline(slope = m, intercept = b, color="red")+
  geom_abline(slope = c, intercept = d, color="blue")+
  scale_x_continuous(limits=c(0,10))+
  scale_y_continuous(limits=c(0,0.7));tp.lines
#Plot of our two linear performance models

#Applying model to empty iButton dataframe, t.dat.full.lmodel:
t.dat.full.lmodel.growth.decay <- t.dat.full #Duplicating t.dat.full dataframe
for(i in 2:nrow(t.dat.full.lmodel.growth.decay)) {if(is.na(t.dat.full.lmodel.growth.decay[i,5])=="TRUE") #running model
  {t.dat.full.lmodel.growth.decay[i,5] <- (t.dat.full.lmodel.growth.decay[i-1,5])*exp(((m * t.dat.full.lmodel.growth.decay[i-1,4]) + b) - ((c*t.dat.full.lmodel.growth.decay[i-1,4])+d))}} #(t.dat.full.lmodel.growth.decay[i-1,5])*exp((((m * t.dat.full.lmodel.growth.decay[i-1,4]) + b) - ((c*t.dat.full.lmodel.growth.decay[i-1,4])+d))/5)}}

#Plot of each bat's loads over time:
p1<- ggplot(aes(x=date, y=log10(gdL), color=site), data=t.dat.full.lmodel.growth.decay) +
  geom_line()+
  facet_wrap(facets = "trans_id", scales='free');p1

#Let's plot density plots of these late hibernation loads and d.lgdL values against measured values.

#Need to construct a dataframe from model output:
end.loads.growth.decay <- t.dat.full.lmodel.growth.decay %>%
  group_by(trans_id) %>%
  filter(date==max(date, na.rm=T)) %>%
  mutate(lgdL.late = log10(gdL), data="Model Output") 
#Final loads predicted by the model. 
end.loads.growth.decay$lgdL.early <- i.dat$lgdL[match(end.loads.growth.decay$trans_id, i.dat$trans_id)]
#Bringin in early lgdL to calculate change in loads. 
end.loads.growth.decay$date.early <- i.dat$date[match(end.loads.growth.decay$trans_id, i.dat$trans_id)]
#Bringing in fall sampling date to calculate sampling time and daily pathogen growth
end.loads.growth.decay$sampling.length <- end.loads.growth.decay$date - end.loads.growth.decay$date.early
#Length of sampling interval
end.loads.growth.decay$d.lgdL <- end.loads.growth.decay$lgdL.late - end.loads.growth.decay$lgdL.early
#Over-winter change in pathogen loads
end.loads.growth.decay$d.lgdL.daily <- end.loads.growth.decay$d.lgdL/as.numeric(end.loads.growth.decay$sampling.length)
#Daily over-winter change in pathogen loads
end.loads.growth.decay <- select(end.loads.growth.decay, trans_id, site, data, lgdL.late, d.lgdL, d.lgdL.daily)
#Reducing dataframe to minimum columns

#And now another dataframe constructed from empirical data only:
i.dat$lgdL.late <- i.dat.late$lgdL[match(i.dat$trans_id, i.dat.late$trans_id)]
#Bringing late hibernation load data to early hibernation dataframe
i.dat$date.late <- i.dat.late$date[match(i.dat$trans_id, i.dat.late$trans_id)]
#Bringing late hibernation date data to early hibernation dataframe
i.dat$sampling.length <- i.dat$date.late - i.dat$date
#Length of sampling period
i.dat$d.lgdL <- i.dat$lgdL.late - i.dat$lgdL
#Over-winter change in pathogen load
i.dat$d.lgdL.daily <- i.dat$d.lgdL/as.numeric(i.dat$sampling.length)
#Daily over-winter change in pathogen load
i.dat.red <- select(i.dat, trans_id, site, lgdL.late, d.lgdL, d.lgdL.daily)
i.dat.red <- filter(i.dat.red, !is.na(d.lgdL.daily) & d.lgdL.daily != 0)
#Reducing dataframe
i.dat.red$data <- "qPCR"
#Data type column for merging
i.dat.red$trans_id <- as.factor(i.dat.red$trans_id)
#Needs to be a factor for merging

end.loads.growth.decay2 <- rbind(end.loads.growth.decay, i.dat.red)
#Dataframe for plotting


#end.loads.growth.decay.density <- ggplot(aes(x=lgdL), data=end.loads.growth.decay) +
 # geom_density(fill="blue", alpha=0.6)+
  #geom_density(aes(x=lgdL), fill="red", alpha=0.6, data=i.dat.all[i.dat.all$season=="hiber_late",])+
  #scale_fill_discrete(labels=c("Model Output", "qPCR Values"));end.loads.growth.decay.density









#I also want a model using just the pathogen growth term, without the decay term, to compare outputs:
t.dat.full.lmodel.growth.only <- t.dat.full #Duplicating dataframe
for(i in 2:nrow(t.dat.full.lmodel.growth.only)) {if(is.na(t.dat.full.lmodel.growth.only[i,5])=="TRUE") #running model
{t.dat.full.lmodel.growth.only[i,5] <- (t.dat.full.lmodel.growth.only[i-1,5])*exp((m * t.dat.full.lmodel.growth.only[i-1,4]) + b)}} # (t.dat.full.lmodel.growth.only[i-1,5])*exp(((m * t.dat.full.lmodel.growth.only[i-1,4]) + b)/5)}}

p2<- ggplot(aes(x=date, y=log10(gdL), color=site), data=t.dat.full.lmodel.growth.only) +
  geom_line()+
  facet_wrap(facets = "trans_id");p2



