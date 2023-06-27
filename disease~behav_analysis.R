library(tidyverse)
library(lme4)
library(effects)
library(emmeans)
library(lmerTest)
library(ggExtra)
library(glmmTMB)
library(betareg)

dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/arousals_torpors_working.csv")
#This dataframe contains a list and summary of every individual's arousal and torpor events. 

dat$start.datetime <- as.POSIXct(strptime(dat$start.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
dat$end.datetime <- as.POSIXct(strptime(dat$end.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
#Formatting date correctly

dat$event.length.days[dat$event.length.days == 0] <- NA
#There are several instances where there is an extremely short bout of "torpor", for less that the sampling interval of the logger.
#These need to be NA when calculating change in torpor bout temperature independent variables

#### Building independent variables of interest #####

#First, going to start simply with the number of arousals, arousal frequency, average torpor bout length, and
#mean torpor bout temperature (already calculated and in the dat dataframe)

#Starting with arousal frequency
ind.summ <- dat %>%
  mutate(c=1) %>%
  filter(behavior=="Arousal") %>%
  group_by(site, trans_id) %>%
  summarise(num.arousals = sum(c))
#This dataframe, ind.summ, stands for "individual summary." It will contain the independent variables for every individual, 
#which will then be linked to disease severity metrics. Currently, it just contains the number of arousals per individual.

#The number of arousals is meaningless without correcting for the length of the sampling period (arousal frequency). 
#In other papers using arousal data, authors remove the first and last torpor bouts from the sampling duration because they could have been 
#artificially shortened due to handling. I will do the same, but there are three bats whose first event is actually an
#arousal, which should be fine to include in the sampling period duration. Additionally, four individuals had arousals
#as their final events, and those should be included in their duration as well. 

first.event <- dat %>%
  group_by(site, trans_id) %>%
  filter(event_num == min(event_num)) %>%
  select(trans_id, behavior)
#Dataframe containing each individual's first event type
last.event <- dat %>%
  group_by(site, trans_id) %>%
  filter(event_num == max(event_num)) %>%
  select(trans_id, behavior)
#Dataframe containing each individual's last event type
first.last <- first.event
colnames(first.last) <- c("site","trans_id","first.event.type")
first.last$last.event.type <- last.event$behavior[match(first.last$trans_id, last.event$trans_id)]
#This dataframe now contains the behavior type of each individual's first and last event. 

tt <- unique(first.last$trans_id[first.last$first.event.type=="Torpor" & first.last$last.event.type=="Torpor"])
#These are the trans_id's of all the individuals whose first and last events were torpor bouts. There are 73 such individuals
ta <- unique(first.last$trans_id[first.last$first.event.type=="Torpor" & first.last$last.event.type=="Arousal"])
#These are the trans_id's of all the individuals whose first event was a torpor, but last was an arousal. There are 4 such individuals
at <- unique(first.last$trans_id[first.last$first.event.type=="Arousal" & first.last$last.event.type=="Torpor"])
#These are the trans_id's of all the individuals whose first event was an arousal, but last was torpor. There are 3 such individuals
#There are no bats that started and ended with an arousal.

tt.durations <- dat %>%
  filter(trans_id %in% tt) %>%
  group_by(trans_id) %>%
  filter(event_num == 2 | event_num == max(event_num)-1) %>%
  mutate(sampling.duration = max(end.datetime) - min(start.datetime)) %>%
  summarise(sampling.duration=mean(sampling.duration))
#This dataframe contains the sampling durations for the tt bats. 

ta.durations <- dat %>%
  filter(trans_id %in% ta) %>%
  group_by(trans_id) %>%
  filter(event_num == 2 | event_num == max(event_num)) %>%
  mutate(sampling.duration = max(end.datetime) - min(start.datetime)) %>%
  summarise(sampling.duration=mean(sampling.duration))
#This dataframe contains the sampling durations for the tt bats. 
  
at.durations <- dat %>%
  filter(trans_id %in% at) %>%
  group_by(trans_id) %>%
  filter(event_num == 1 | event_num == max(event_num)-1) %>%
  mutate(sampling.duration = max(end.datetime) - min(start.datetime)) %>%
  summarise(sampling.duration=mean(sampling.duration))
#This dataframe contains the sampling durations for the at bats. 

durations <- rbind(tt.durations, ta.durations, at.durations)
#Combining the three duration dataframes

ind.summ$sampling.duration.days <- as.numeric(durations$sampling.duration[match(ind.summ$trans_id, durations$trans_id)])
#Bringing the sampling duration data into the individual summary dataframe

ind.summ$arousal.freq.days <- as.numeric(ind.summ$num.arousals/ind.summ$sampling.duration.days)
#This is the arousal frequency or, on average, the proportion of the sampling time that passed until an arousal. 


#Moving on to average torpor bout length.

#Again, if a torpor bout was the first, last, or both the first and last events for a bat, these torpor bouts cannot be 
#used to calculate average lengths because they could have been artificially changed by handling. 
tt.average.torpor.lengths <- dat %>%
  filter(trans_id %in% tt & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num >1 & event_num < max(event_num)) %>%
  summarise(mean.torpor.length.days = mean(event.length.days, na.rm=T))
#This dataframe contains the mean torpor bout length for the tt bats. 

ta.average.torpor.lengths <- dat %>%
  filter(trans_id %in% ta & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num >1) %>%
  summarise(mean.torpor.length.days = mean(event.length.days, na.rm=T))
#This dataframe contains the mean torpor bout length for the ta bats. 

at.average.torpor.lengths <- dat %>%
  filter(trans_id %in% at & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num < max(event_num)) %>%
  summarise(mean.torpor.length.days = mean(event.length.days, na.rm=T))
#This dataframe contains the mean torpor bout length for the at bats. 

mean.torpor.length <- rbind(tt.average.torpor.lengths, ta.average.torpor.lengths, at.average.torpor.lengths)
#This dataframe contains each bat's average torpor bout length. 

ind.summ$mean.torpor.length.days <- mean.torpor.length$mean.torpor.length.days[match(ind.summ$trans_id, mean.torpor.length$trans_id)]
#Looping average torpor bout length data into the individual summary dataframe

#Finally, just need to loop in the torpor bout temp data:
ind.summ$mean.torpor.temp <- dat$mean.torpor.temp[match(ind.summ$trans_id, dat$trans_id)]


# The next independent variables I want to construct are a little bit more complicated. I want to capture changes in 
# roosting temperature after arousal events, whether they were caused by movements or otherwise. To do so, I am first going
# to construct a change in average torpor bout temperature metric, which is just the difference in mean temperature during
# torpor bouts after an arousal event. I will then calculate the weighted average of those values, using the torpor bout
# duration as weights. I will also just construct the unweighted average, because why not? 

# I will also calculate a simple sum of the temperature changes, which should give a sort of deviation
# from the mean temperature of the first torpor bout. 

# Finally, I want to construct a sum of the deviations from the mean of the first torpor bout, corrected by the 
# sampling duration. Basically, if we took the average temperature of the first torpor bout, I want to then calculate
# the difference between that value and the average temperature of each subsequent torpor bout, and then find the 
# duration-corrected average of those values. This metric is cool in that it has memory built-in, so it's like comparing
# what happens when a bat moves or changes temperature to what would happen if it stayed in its first roosting location
# for the entirety of hibernation. 

#First, the change in averaage torpor bout metric:
dat$d.mean.torpor.temp <- NA
#This is the column in which I'll store the change in temperature data. 

for(i in 3:nrow(dat)) {if(dat[i,3] == "Arousal") {dat[i,15] <- NA}
  else{if(dat[i,2] != dat[i-1,2]) {dat[i,15] <- NA}
    else{if(dat[i,2] != dat[i-2,2]) {dat[i,15] <- NA}
      else{if(is.na(dat[i-2,12])=="TRUE") {dat[i,15] <- NA}
        else{dat[i,15] <- dat[i,7] - dat[i-2,7]}}
      }
    }
  }
#This for-loop calculates the change in mean torpor bout temperature following arousal events. Of course, no change in temp
#value is given to the first torpor bout. 

mean.d.torpor.temps.weighted <- dat %>%
  filter(event.length.days > 0) %>%
  group_by(trans_id) %>%
  summarise(mean.d.torpor.temps.weighted = weighted.mean(x=d.mean.torpor.temp, w=event.length.days, na.rm=T))
#This dataframe contains the weighted change in torpor temps for each bat. 


mean.d.torpor.temps.unweighted <- dat %>%
  group_by(trans_id) %>%
  summarise(mean.d.torpor.temps.unweighted = mean(d.mean.torpor.temp, na.rm=T))
#This dataframe contains the unweighted change in torpor temps for each bat. 

ind.summ$mean.d.torpor.temp.weighted <- mean.d.torpor.temps.weighted$mean.d.torpor.temps.weighted[match(ind.summ$trans_id, mean.d.torpor.temps.weighted$trans_id)]
ind.summ$mean.d.torpor.temp.unweighted <- mean.d.torpor.temps.unweighted$mean.d.torpor.temps.unweighted[match(ind.summ$trans_id, mean.d.torpor.temps.unweighted$trans_id)]
#Looping this data into the individual summary dataframe. 


#Now constructing the deviations from first torpor temp metric:

first.torpor.temps <- dat %>%
  filter(behavior == "Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num == min(event_num)) %>%
  summarise(mean.temp.first.torpor = mean(mean.temp))
#This is a dataframe of the mean temperature of each bat's first torpor bout. 

dat$mean.temp.first.torpor <- first.torpor.temps$mean.temp.first.torpor[match(dat$trans_id, first.torpor.temps$trans_id)]
#Bringing the data into the dat dataframe

dat$temp.dev.first.torpor <- NA
dat$temp.dev.first.torpor[dat$behavior=="Torpor"] <- dat$mean.temp[dat$behavior=="Torpor"] - dat$mean.temp.first.torpor[dat$behavior=="Torpor"]
#This column now contains the deviation of each torpor bout's mean temperature from the mean temperature of the first torpor bout for each bat. 
dat$temp.dev.first.torpor[dat$temp.dev.first.torpor==0] <- NA
#There are instances where the deviation value is 0 because it is the first torpor bout, so of course the deviation is 0. 
#Making these NA because I don't want them included in any averages. 

mean.temp.dev.first.torpor.weighted <- dat %>%
  group_by(trans_id) %>%
  filter(!is.na(event.length.days)) %>%
  summarise(mean.temp.dev.first.torpor.weighted = weighted.mean(x=temp.dev.first.torpor, w=event.length.days, na.rm=T))
#This dataframe contains the weighted average (again weighted by torpor bout length) of the deviations of each torpor
#bout's mean temperature from that of the first torpor bout. 
mean.temp.dev.first.torpor.unweighted <- dat %>%
  group_by(trans_id) %>%
  summarise(mean.temp.dev.first.torpor.unweighted = mean(temp.dev.first.torpor, na.rm=T))
#This is the unweighted analog of above

ind.summ$mean.temp.dev.first.torpor.weighted <- NA
ind.summ$mean.temp.dev.first.torpor.weighted <- mean.temp.dev.first.torpor.weighted$mean.temp.dev.first.torpor.weighted[match(ind.summ$trans_id, mean.temp.dev.first.torpor.weighted$trans_id)]
ind.summ$mean.temp.dev.first.torpor.unweighted <- NA
ind.summ$mean.temp.dev.first.torpor.unweighted <- mean.temp.dev.first.torpor.unweighted$mean.temp.dev.first.torpor.unweighted[match(ind.summ$trans_id, mean.temp.dev.first.torpor.unweighted$trans_id)]
#Bringing the weighted and unweighted averages into the individual summary dataframe. 




# I also want to use the ambient temperature data recorded by psychrometers/pro v2 data loggers to construct a metric of an individual's microclimate deviation from 
# available microclimate conditions. I'm going to do this in two ways: first, for each site, I'll identify the warmest/coldest available sections from psych/hobo data and
# calculate a DAILY deviation from those dataset. Basically, for every day a bat's transmitter was recording data, I'll calculate a deviation from the average daily 
# temperature recorded by the psych/hobo in the coldest and warmest sections and then average all of those daily deviation values for each bat. Second, I'll do the same
# thing, but instead of using data from the coldest/warmest section, I'll use psych/hobo data from the section in which the bats were initially captured in early hibernation,
# if available. 

#First, need to read in the raw transmitter dataset and filter away any arousal data: 
raw.trans.data <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv") %>%
  filter(behavior != "Arousal") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d"))

#As well as the psych/hobo dataset: 
env.data <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/psychrometer_hobo_data_26JUN2023.csv") %>%
  mutate(date = as.Date(date, format="%Y-%m-%d"))

#Need daily average, non-arousal temperature values for each individual:
trans.daily.temp.avgs <- raw.trans.data %>% 
  group_by(site, trans_id, date) %>% 
  summarise(temp.mean.daily = mean(temp))
#This dataframe contains daily torpor temperature averages for each bat. 


#Also need to loop in the data on the sections in which each bat was captured in early hibernation from the disease dataframe:
dis.dat.master <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_5JUN2023.csv")
#This is the disease data for bats in this study available as of the 5th of June, 2023. This dataset does not include all infection
#data that will ultimately be available, as Jeff Foster et al. are still processing/correcting samples. If the mean_gdL value
#is NA, the sample has not yet been processed and we haven't yet received data back. 

dis.early <- filter(dis.dat.master, season=="hiber_earl") #This is the disease dataframe reduced to just early hibernation data.
trans.daily.temp.avgs$origin.section <- dis.early$section[match(trans.daily.temp.avgs$trans_id, dis.early$trans_id)]
#Each bat's first section, or the section in which it was captured in early hibernation, is now matched into the average daily temp dataframe

#Now we can calculate the average daily temperatures for each logger. However, first, we need to remove data from a few loggers that didn't record the entire
#sampling period. To do so, I'm first making a new column with a 'uniq.id' for each logger based on the site, section, and logger_model, as no section should
#have more than one logger of the same model:
env.data$uniq.id <- paste(env.data$site, env.data$section, env.data$logger_model)
#Now, I'll summarise each logger's maximum date to identify those that didn't record data for the entire study period:
env.data.ends <- env.data %>% 
  group_by(uniq.id) %>% 
  summarise(end.date=max(date)) %>%
  filter(end.date >= "2022-03-01")
length(unique(env.data$uniq.id)) - length(unique(env.data.ends$uniq.id))
#There are four that prematurely ended recording. If it ended before March, it prematurely ended. The above dataframe contains the loggers that didn't prematurely end. 

env.data <- filter(env.data, uniq.id %in% unique(env.data.ends$uniq.id))
#This filters away all logger data from loggers that prematurely stopped logging. 

env.data.daily.avg <- env.data %>%
  group_by(uniq.id, site, section, logger_model, date) %>% 
  summarise(temp.mean.daily = mean(temp_1))
#This dataframe contains each logger's daily average temperatures. 

#I now need to identify the coldest and warmest sections/loggers in each site: 
overall.avg.logger.temps <- env.data %>% 
  group_by(uniq.id, site) %>% 
  summarise(temp.mean.overall = mean(temp_1), temp.median.overall = median(temp_1)) %>% 
  ungroup() %>% 
  group_by(site) %>%
  mutate(temp.mean.max = max(temp.mean.overall), temp.mean.min = min(temp.mean.overall)) %>% 
  filter(temp.mean.overall == temp.mean.max | temp.mean.overall == temp.mean.min) 
#This dataframe contains the mean and median temperature values for the warmest and coldest section-loggers in each site. 

env.data.daily.avg.warmest <- env.data.daily.avg %>% 
  filter(uniq.id %in% unique(overall.avg.logger.temps$uniq.id[overall.avg.logger.temps$temp.mean.overall == overall.avg.logger.temps$temp.mean.max])) %>%
  mutate(site.date = paste(site, date, sep=" ")) #this column will be used to match this temp data into the transmitter dataframe. 
env.data.daily.avg.coldest <- env.data.daily.avg %>% 
  filter(uniq.id %in% unique(overall.avg.logger.temps$uniq.id[overall.avg.logger.temps$temp.mean.overall == overall.avg.logger.temps$temp.mean.min])) %>%
  mutate(site.date = paste(site, date, sep=" ")) #this column will be used to match this temp data into the transmitter dataframe. 
#These two dataframes now only contain daily average temperature data from those loggers in warmest and coldest sections, respectively. 

#Now we can create two new columns in trans.daily.temp.avgs that contain the daily average temps from the warmest and coldest sections:
trans.daily.temp.avgs$site.date <- paste(trans.daily.temp.avgs$site, trans.daily.temp.avgs$date, sep=" ") #Matching on this column. 
trans.daily.temp.avgs$temp.mean.daily.warmest.section <- env.data.daily.avg.warmest$temp.mean.daily[match(trans.daily.temp.avgs$site.date, env.data.daily.avg.warmest$site.date)]
trans.daily.temp.avgs$temp.mean.daily.coldest.section <- env.data.daily.avg.coldest$temp.mean.daily[match(trans.daily.temp.avgs$site.date, env.data.daily.avg.coldest$site.date)]

#Finally, I want to match in the logger data from each bat's origin section:
trans.daily.temp.avgs$temp.mean.daily.origin.section <- env.data.daily.avg$temp.mean.daily[match((paste(trans.daily.temp.avgs$site, trans.daily.temp.avgs$origin.section, 
                                                                                                        trans.daily.temp.avgs$date)),paste(env.data.daily.avg$site, 
                                                                                                                                           env.data.daily.avg$section,
                                                                                                                                           env.data.daily.avg$date))]
#Of course, not all bats had logger data available from their origin section, so they may have NA values in that column. 

#Now, we can calculate our temperature deviations of interest:
trans.daily.temp.avgs$d.warmest <- trans.daily.temp.avgs$temp.mean.daily - trans.daily.temp.avgs$temp.mean.daily.warmest.section #The daily deviation from the warmest section's daily temperature
trans.daily.temp.avgs$d.coldest <- trans.daily.temp.avgs$temp.mean.daily - trans.daily.temp.avgs$temp.mean.daily.coldest.section#The daily deviation from the coldest section's daily temperature
trans.daily.temp.avgs$d.origin <- trans.daily.temp.avgs$temp.mean.daily - trans.daily.temp.avgs$temp.mean.daily.origin.section #The daily deviation from the origin section's daily temperature 

#And now we will average these data for every individual:
trans.daily.devs.summ <- trans.daily.temp.avgs %>% 
  group_by(trans_id) %>% 
  summarise(mean.temp.dev.warmest.sec = mean(d.warmest, na.rm=T), mean.temp.dev.coldest.sec = mean(d.coldest, na.rm=T), mean.temp.dev.origin.sec = mean(d.origin))

#Match into ind.summ:
ind.summ$mean.temp.dev.warmest.sec <- trans.daily.devs.summ$mean.temp.dev.warmest.sec[match(ind.summ$trans_id, trans.daily.devs.summ$trans_id)]
ind.summ$mean.temp.dev.coldest.sec <- trans.daily.devs.summ$mean.temp.dev.coldest.sec[match(ind.summ$trans_id, trans.daily.devs.summ$trans_id)]
ind.summ$mean.temp.dev.origin.sec <- trans.daily.devs.summ$mean.temp.dev.origin.sec[match(ind.summ$trans_id, trans.daily.devs.summ$trans_id)]





#### Matching in disease data ####

dis.dat.master <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/inf_data_5JUN2023.csv")
#This is the disease data for bats in this study available as of the 5th of June, 2023. This dataset does not include all infection
#data that will ultimately be available, as Jeff Foster et al. are still processing/correcting samples. If the mean_gdL value
#is NA, the sample has not yet been processed and we haven't yet received data back. 

trans_meta<- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv")
#Transmitter metadata database

dis.dat.master$date <- as.Date(dis.dat.master$date, format="%m/%d/%y")
#Fixing date column.

dis.early <- filter(dis.dat.master, season=="hiber_earl")
dis.late <- filter(dis.dat.master, season=="hiber_late")
#separating the master disease database by early and late hibernation. 

dis.early$trans_model <- trans_meta$model[match(dis.early$trans_id, trans_meta$id)]
dis.early$trans_weight_g <- trans_meta$weight_g[match(dis.early$trans_id, trans_meta$id)]
#Bringing in the transmitter model and weight metadata


bad.bands<-select(dis.early, site, trans_id, date, band)
colnames(bad.bands) <- c("site", "trans_id", "date_early", "band_early")
bad.bands$date_late <- dis.late$date[match(bad.bands$trans_id, dis.late$trans_id)]
bad.bands$band_late <- dis.late$band[match(bad.bands$trans_id, dis.late$trans_id)]
bad.bands<-filter(bad.bands, band_early != band_late)
#These are inconsistencies between early and late hibernation in the band # recorded for a bat with the same trans_id. 
#The two bats in South Lake had lost their bands (but not their transmitters) over the course of winter, so their late hibernation band was the one we replaced it with. 

dis.df <- select(dis.early, trans_id, trans_model, trans_weight_g, site, section, state, band, swab_id, sex, mass, uv_orange, uv_right, gd_wall, gd_wall2, wing_score, wing_score2, mean_gdL)
#This is going to be the dataframe with which I combine the ind.summ dataframe. 
#I use the mean_gdL column instead of the mean_gdL column to preserve some information. In the mean_gdL column, if the sample was run but
#no DNA was detected, it is given a value of 0 rather than NA. If the sample was not run yet (we don't have the data back), this value is
#instead an NA. In the true_mean_gdL column, the value is NA in both of those cases, meaning I cannot distinguish between samples without
#detectable DNA and samples that haven't been run yet. This is the only difference in the value of the two columns. Thus, I choose to 
#preserve that information by using mean_gdL. 

colnames(dis.df) <- c("trans_id", "trans.model", "trans.weight.g", "site", "section.early", "state", "band", "swab.id.early", "sex", "mass.early","uv.orange.early", 
                      "uv.right.early", "gd.wall.early","gd.wall2.early", "wing.score.early","wing.score2.early", "gdL.early")
#Re-naming columns so that their early disease information is identifiable.

dis.df$mass.late <- dis.late$mass[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$swab.id.late <- dis.late$swab_id[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$uv.orange.late <- dis.late$uv_orange[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$uv.right.late <- dis.late$uv_right[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$gd.wall.late <- dis.late$gd_wall[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$gd.wall2.late <- dis.late$gd_wall2[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$wing.score.late <- dis.late$wing_score[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$wing.score2.late <- dis.late$wing_score2[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$gdL.late <- dis.late$mean_gdL[match(dis.df$trans_id, dis.late$trans_id)]
#Matching in the disease data matching on transmitter ID. 

dis.df$sex.late <- dis.late$sex[match(dis.df$trans_id, dis.late$trans_id)]
mal.sex <- filter(dis.df, sex != sex.late)
#The sex of this bat miraculously changed over the course of hibernation. 

#Constructing disease severity metrics of interest: mean early/late gd score, mean early/late wing score, mean early/late
#UV score, changes in those metrics, and difference in mass. 

dis.df$uv.orange.early<-as.numeric(dis.df$uv.orange.early)
dis.df$uv.orange.late<-as.numeric(dis.df$uv.orange.late)
dis.df$uv.right.early<-as.numeric(dis.df$uv.right.early)
dis.df$uv.right.late<-as.numeric(dis.df$uv.right.late)
dis.df$gd.wall.early <- as.numeric(dis.df$gd.wall.early)
dis.df$gd.wall.late <- as.numeric(dis.df$gd.wall.late)
dis.df$gd.wall2.early <- as.numeric(dis.df$gd.wall2.early)
dis.df$gd.wall2.late <- as.numeric(dis.df$gd.wall2.late)
dis.df$wing.score.early <- as.numeric(dis.df$wing.score.early)
dis.df$wing.score.late <- as.numeric(dis.df$wing.score.late)
dis.df$wing.score2.early <- as.numeric(dis.df$wing.score2.early)
dis.df$wing.score2.late <- as.numeric(dis.df$wing.score2.late)
#Converting all categorical score data into numerics for averaging. 

dis.df$uv.score.mean.early <- (dis.df$uv.orange.early+dis.df$uv.right.early)/2
#Mean early hibernation uv score
dis.df$uv.score.mean.late <- (dis.df$uv.orange.late+dis.df$uv.right.late)/2
#Mean late hibernation uv score
dis.df$d.uv.score <- dis.df$uv.score.mean.late - dis.df$uv.score.mean.early
#Change between early and late hibernation in mean uv score. 

dis.df$gd.score.mean.early <- (dis.df$gd.wall.early + dis.df$gd.wall2.early)/2
#Mean early hibernation gd score
dis.df$gd.score.mean.late <- (dis.df$gd.wall.late + dis.df$gd.wall2.late)/2
#Mean late hibernation gd score
dis.df$d.gd.score <- dis.df$gd.score.mean.late - dis.df$gd.score.mean.early
#Change between early and late hibernation in mean gd score. 

dis.df$wing.score.mean.early <- (dis.df$wing.score.early + dis.df$wing.score2.early)/2
#Mean early hibernation wing score
dis.df$wing.score.mean.late <- (dis.df$wing.score.late + dis.df$wing.score2.late)/2
#Mean late hibernation wing score
dis.df$d.wing.score <- dis.df$wing.score.mean.late - dis.df$wing.score.mean.early
#Change between early and late hibernation mean wings score. 

dis.df$d.mass <- dis.df$mass.late - dis.df$mass.early
#Change in mass over course of hibernation.

dis.df$sampling.duration.days <- ind.summ$sampling.duration.days[match(dis.df$trans_id, ind.summ$trans_id)]
dis.df$d.mass.daily <- dis.df$d.mass/dis.df$sampling.duration.days
#This is a measure of weight change that corrects for differences in the sampling interval, since bats from all sites were not sampled at the same time.
#it is a "daily weight loss" metric. 

dis.df$d.gdL <- dis.df$gdL.late - dis.df$gdL.early
#Over-winter pathogen growth.

#Need to transform all of the gdL values of 0 (no detectable fungus) to non-zero values by adding an extremely 
#small constant that corresponds to the 40CT qPCR cutoff value:
dis.df$gdL.early <- dis.df$gdL.early + 10^((40-22.04942)/-3.34789)
dis.df$gdL.late <- dis.df$gdL.late + 10^((40-22.04942)/-3.34789)
#Constant added to both early and late hibernation load values. 

dis.df$lgdL.early <- log10(dis.df$gdL.early)
dis.df$lgdL.late <- log10(dis.df$gdL.late)
#log-10 transformed early and late hibernation fungal loads 

dis.df$d.lgdL <- log10(dis.df$d.gdL + 1)
#This is the log10 change in gdL values from early to late hibernation. IMPORTANT: because many bats had a reduction in 
#infection between early and late hibernation, resulting in negative d.gdL values, a constant of 1.0 had to be added
#so that all values could be log-transformed. Therefore, in this d.lgdL value, values of 0.0 correspond to no change,
#whereas values below/above 0.0 correspond to pathogne reduction/growth, respectively. 

## Can now match in all the transmitter summary data for each individual: 
dis.df <- left_join(dis.df, ind.summ[,2:14], by="trans_id")
#Merged. 




#### Variation across sites in behavior metrics ####

dis.df$site <- factor(dis.df$site, levels = c("SOUTH LAKE MINE","GRAPHITE MINE","CP TUNNEL", "BLACKBALL", "ZIMMERMAN", "MEAD MINE","ELROY SPARTA"))
#Make sure site is a factor. I have ordered the levels of this factor from coldest mean torpor bout temp to warmest. 


## Arousal frequency

ar.freq.m <- betareg(arousal.freq.days ~ site, data=dis.df, link="logit");summary(ar.freq.m)
#Beta regression of arousal frequencies. 

ar.freq.summ <- dis.df %>%
  group_by(site) %>%
  summarise(ar.freq.mean = mean(arousal.freq.days, na.rm=T), ar.freq.sd = sd(arousal.freq.days, na.rm=T)) %>%
  mutate(hi.sd = ar.freq.mean + ar.freq.sd, lo.sd = ar.freq.mean - ar.freq.sd)
#Summary table of arousal frequency data across sites. Mean and SD range. 

ar.freq.site.p <- ggMarginal((ggplot(aes(x=site, y=ar.freq.mean), data=ar.freq.summ) +
  geom_jitter(aes(x=site, y=arousal.freq.days, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Arousal Frequency (days)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top")), type="histogram", fill="darkgray", bins=15); ar.freq.site.p
#Plotted



## Mean torpor length

mean.torpor.length.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.torpor.length.mean = mean(mean.torpor.length.days, na.rm=T), mean.torpor.length.sd = sd(mean.torpor.length.days, na.rm=T)) %>%
  mutate(hi.sd = mean.torpor.length.mean + mean.torpor.length.sd, lo.sd = mean.torpor.length.mean - mean.torpor.length.sd)
#Summary table of mean torpor bout length data across sites. Mean and SD range. 

mean.torpor.length.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.torpor.length.mean), data=mean.torpor.length.summ) +
  geom_jitter(aes(x=site, y=mean.torpor.length.days, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Torpor Bout Length (days)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); mean.torpor.length.site.p
#Plotted




## Mean torpor temp

mean.torpor.temp.m <- lm(mean.torpor.temp ~ site, data=dis.df);summary(mean.torpor.temp.m)
#Model of site differences in mean torpor temperatures. 

mean.torpor.temp.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.torpor.temp.mean = mean(mean.torpor.temp, na.rm=T), mean.torpor.temp.sd = sd(mean.torpor.temp, na.rm=T)) %>%
  mutate(hi.sd = mean.torpor.temp.mean + mean.torpor.temp.sd, lo.sd = mean.torpor.temp.mean - mean.torpor.temp.sd)
#Summary table of mean torpor bout temperature data across sites. Mean and SD range. 

mean.torpor.temp.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.torpor.temp.mean), data=mean.torpor.temp.summ) +
  geom_jitter(aes(x=site, y=mean.torpor.temp, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Torpor Bout Temperature (Celsius)") +
    scale_y_continuous(limits=c(0,9), breaks=seq(0,12,1), labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15);mean.torpor.temp.site.p
#Plotted




## Mean change in mean torpor bout temperature

mean.d.torpor.temp.weighted.m <- lm(mean.d.torpor.temp.weighted ~ site, data=dis.df);summary(mean.d.torpor.temp.weighted.m)
#Model of change in torpor bout temperatures across sites (weighted)

#Weighted first:
mean.d.torpor.temp.w.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.d.torpor.temp.w.mean = mean(mean.d.torpor.temp.weighted, na.rm=T), mean.d.torpor.temp.w.sd = sd(mean.d.torpor.temp.weighted, na.rm=T)) %>%
  mutate(hi.sd = mean.d.torpor.temp.w.mean + mean.d.torpor.temp.w.sd, lo.sd = mean.d.torpor.temp.w.mean - mean.d.torpor.temp.w.sd)
#Summary table of change in torpor bout temperature data across sites. Mean and SD range. 

mean.d.torpor.temp.w.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.d.torpor.temp.w.mean), data=mean.d.torpor.temp.w.summ) +
  geom_jitter(aes(x=site, y=mean.d.torpor.temp.weighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Change in Torpor Bout \n Temperature Following Arousal (Weighted)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); mean.d.torpor.temp.w.site.p
#Plotted

#Now unweighted:

mean.d.torpor.temp.unweighted.m <- lm(mean.d.torpor.temp.unweighted ~ site, data=dis.df);summary(mean.d.torpor.temp.unweighted.m)
#Model

mean.d.torpor.temp.uw.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.d.torpor.temp.uw.mean = mean(mean.d.torpor.temp.unweighted, na.rm=T), mean.d.torpor.temp.uw.sd = sd(mean.d.torpor.temp.unweighted, na.rm=T)) %>%
  mutate(hi.sd = mean.d.torpor.temp.uw.mean + mean.d.torpor.temp.uw.sd, lo.sd = mean.d.torpor.temp.uw.mean - mean.d.torpor.temp.uw.sd)
#Summary table of change in torpor bout temperature data across sites. Mean and SD range. 

mean.d.torpor.temp.uw.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.d.torpor.temp.uw.mean), data=mean.d.torpor.temp.uw.summ) +
  geom_jitter(aes(x=site, y=mean.d.torpor.temp.unweighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Change in Torpor Bout \n Temperature Following Arousal (Unweighted)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); mean.d.torpor.temp.uw.site.p
#Plotted






## Mean deviation in mean torpor bout temperature from first torpor bout

mean.temp.dev.w.m <- lm(mean.temp.dev.first.torpor.weighted ~ site, data=dis.df);summary(mean.temp.dev.w.m)
#Model of deviation in torpor bout temperatures from first torpor bout across sites. 

#Weighted first:
mean.temp.dev.w.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.w.mean = mean(mean.temp.dev.first.torpor.weighted, na.rm=T), mean.temp.dev.w.sd = sd(mean.temp.dev.first.torpor.weighted, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.w.mean + mean.temp.dev.w.sd, lo.sd = mean.temp.dev.w.mean - mean.temp.dev.w.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.w.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.temp.dev.w.mean), data=mean.temp.dev.w.summ) +
  geom_jitter(aes(x=site, y=mean.temp.dev.first.torpor.weighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Temp Deviation from \n First Arousal Bout (Weighted)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); mean.temp.dev.w.site.p
#Plotted

#Now unweighted:

mean.temp.dev.uw.m <- lm(mean.temp.dev.first.torpor.unweighted ~ site, data=dis.df);summary(mean.temp.dev.uw.m)
#Model

mean.temp.dev.uw.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.uw.mean = mean(mean.temp.dev.first.torpor.unweighted, na.rm=T), mean.temp.dev.uw.sd = sd(mean.temp.dev.first.torpor.unweighted, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.uw.mean + mean.temp.dev.uw.sd, lo.sd = mean.temp.dev.uw.mean - mean.temp.dev.uw.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.uw.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.temp.dev.uw.mean), data=mean.temp.dev.uw.summ) +
  geom_jitter(aes(x=site, y=mean.temp.dev.first.torpor.unweighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Mean Temp Deviation from \n First Arousal Bout (Unweighted)") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); mean.temp.dev.uw.site.p
#Plotted





#Now looking at average daily temperature deviation from the average daily temperature of the warmest section:

mean.temp.dev.warmest.sec.m <- lm(mean.temp.dev.warmest.sec ~ site, data=dis.df);summary(mean.temp.dev.warmest.sec.m)

mean.temp.dev.warmest.sec.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.warmest.sec.mean = mean(mean.temp.dev.warmest.sec, na.rm=T), mean.temp.dev.warmest.sec.sd = sd(mean.temp.dev.warmest.sec, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.warmest.sec.mean + mean.temp.dev.warmest.sec.sd, lo.sd = mean.temp.dev.warmest.sec.mean - mean.temp.dev.warmest.sec.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.warmest.sec.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.temp.dev.warmest.sec.mean), data=mean.temp.dev.warmest.sec.summ) +
                                         geom_jitter(aes(x=site, y=mean.temp.dev.warmest.sec, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
                                         geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
                                         geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
                                         scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
                                         labs(x=NULL, y="Mean Daily Temp Deviation \n from Warmest Section") +
                                         theme(
                                           axis.text.y = element_text(size=13),
                                           axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                                           axis.title = element_text(size=15),
                                           plot.margin=margin(10,10,0,30),
                                           legend.title = element_text(size=13),
                                           legend.text = element_text(size=13),
                                           legend.position = "top"
                                         )), type="histogram", fill="darkgray", bins=15); mean.temp.dev.warmest.sec.site.p
#Plotted






#Now looking at average daily temperature deviation from the average daily temperature of the coldest section:

mean.temp.dev.coldest.sec.m <- lm(mean.temp.dev.coldest.sec ~ site, data=dis.df);summary(mean.temp.dev.coldest.sec.m)

mean.temp.dev.coldest.sec.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.coldest.sec.mean = mean(mean.temp.dev.coldest.sec, na.rm=T), mean.temp.dev.coldest.sec.sd = sd(mean.temp.dev.coldest.sec, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.coldest.sec.mean + mean.temp.dev.coldest.sec.sd, lo.sd = mean.temp.dev.coldest.sec.mean - mean.temp.dev.coldest.sec.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.coldest.sec.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.temp.dev.coldest.sec.mean), data=mean.temp.dev.coldest.sec.summ) +
                                                  geom_jitter(aes(x=site, y=mean.temp.dev.coldest.sec, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
                                                  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
                                                  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
                                                  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
                                                  labs(x=NULL, y="Mean Daily Temp Deviation \n from Coldest Section") +
                                                  theme(
                                                    axis.text.y = element_text(size=13),
                                                    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                                                    axis.title = element_text(size=15),
                                                    plot.margin=margin(10,10,0,30),
                                                    legend.title = element_text(size=13),
                                                    legend.text = element_text(size=13),
                                                    legend.position = "top"
                                                  )), type="histogram", fill="darkgray", bins=15); mean.temp.dev.coldest.sec.site.p
#Plotted







#Now looking at average daily temperature deviation from the average daily temperature of the origin section:

mean.temp.dev.origin.sec.m <- lm(mean.temp.dev.origin.sec ~ site, data=dis.df);summary(mean.temp.dev.origin.sec.m)

mean.temp.dev.origin.sec.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.origin.sec.mean = mean(mean.temp.dev.origin.sec, na.rm=T), mean.temp.dev.origin.sec.sd = sd(mean.temp.dev.origin.sec, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.origin.sec.mean + mean.temp.dev.origin.sec.sd, lo.sd = mean.temp.dev.origin.sec.mean - mean.temp.dev.origin.sec.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.origin.sec.site.p <- ggMarginal((ggplot(aes(x=site, y=mean.temp.dev.origin.sec.mean), data=mean.temp.dev.origin.sec.summ) +
                                                  geom_jitter(aes(x=site, y=mean.temp.dev.origin.sec, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
                                                  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
                                                  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
                                                  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
                                                  labs(x=NULL, y="Mean Daily Temp Deviation \n from Origin Section") +
                                                  theme(
                                                    axis.text.y = element_text(size=13),
                                                    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                                                    axis.title = element_text(size=15),
                                                    plot.margin=margin(10,10,0,30),
                                                    legend.title = element_text(size=13),
                                                    legend.text = element_text(size=13),
                                                    legend.position = "top"
                                                  )), type="histogram", fill="darkgray", bins=15); mean.temp.dev.origin.sec.site.p
#Plotted


#### Variation across sites in disease metrics ####

## Pathogen loads

#Summary of pathogen load data availability:

num.samples.taken <- dis.dat.master %>%
  mutate(c=1) %>%
  group_by(site, season) %>%
  summarise(samples.taken = sum(c)) %>%
  mutate(uniq = paste(site, season))
#Table of number of swabs taken per site/season

load.data.rec <- dis.dat.master %>%
  filter(!is.na(mean_gdL)) %>%
  mutate(c=1) %>%
  group_by(site, season) %>%
  summarise(data_received = sum(c)) %>%
  mutate(uniq=paste(site, season))
#Table of number of swabs for which data was received per site/season
  
load.data.summ <- num.samples.taken
load.data.summ$load.data.received <- load.data.rec$data_received[match(load.data.summ$uniq, load.data.rec$uniq)]
load.data.summ$load.data.received[is.na(load.data.summ$load.data.received)]<-0
load.data.summ<-select(load.data.summ, site, season, samples.taken, load.data.received)
load.data.summ$missing.samples <- load.data.summ$samples.taken - load.data.summ$load.data.received
#Combined summary dataframe with number of samples missing and proportion received column. 



#Now some density plots of pathogen load data by site/season

dis.dat.master$lgdL<-log10(dis.dat.master$true_mean_gdL)
#This is the log-transformed fungal load data that DOES NOT INCLUDE un-ran samples or samples without detectable fungus. 

p.pl <- list()
site.uniq <- unique(dis.dat.master$site)
for(i in 1:length(site.uniq)){
  p.pl[[i]] <- list()
  dat <- subset(dis.dat.master, site==site.uniq[i])
  p.pl[[i]][[1]] <- ggMarginal((ggplot(aes(x=lgdL, y=season, color=season), data=dat) +
                               geom_jitter(height=0.05)+
                               scale_x_continuous(limits=c(-8,0))+
                               ggtitle(dat$site)), type="density", groupColour = TRUE, groupFill = TRUE)
}
#List of 7 plots (one for each site) with both a jitterplot and marginal density plot of early and late fungal loads (if available)

#Combined plot: 

dis.dat.master$season<-as.factor(dis.dat.master$season)
p.pl2 <- ggplot(aes(x=site, y=lgdL, color=season, group=season), data=dis.dat.master) +
  geom_point(position=position_dodge(width=0.7));p.pl2

## Change in pathogen load: 

non.infected <- filter(dis.df, gdL.early == 1.000000e-08 & gdL.late==1.000000e-08)
#These 7 bats never showed any infection by qPCR. 

d.lgdL.p <- ggplot(aes(x=d.lgdL, y=site, color=site), data=dis.df[dis.df$site!="BLACKBALL" & dis.df$site!="CP TUNNEL" &
                                                                    !(dis.df$trans_id %in% non.infected$trans_id),]) +
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=1) + 
  geom_jitter(height=0.05);d.lgdL.p







## Change in UV score:

d.uv.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.uv.mean = mean(d.uv.score, na.rm=T), d.uv.sd = sd(d.uv.score, na.rm=T)) %>%
  mutate(hi = d.uv.mean + d.uv.sd, lo = d.uv.mean - d.uv.sd)
#Dataframe with the mean +/- standard deviation of the change in UV score from early to late hibernation

d.uv.site.p <- ggMarginal((ggplot(aes(x=site, y=d.uv.mean), data=d.uv.summ) +
  geom_jitter(aes(x=site, y=d.uv.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Change in Mean UV Score")+
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=8); d.uv.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.



## Change in pd score

d.gd.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.gd.mean = mean(d.gd.score, na.rm=T), d.gd.sd = sd(d.gd.score, na.rm=T)) %>%
  mutate(hi = d.gd.mean + d.gd.sd, lo = d.gd.mean - d.gd.sd)
#Dataframe with the mean +/- standard deviation of the change in gd score from early to late hibernation

d.gd.site.p <- ggMarginal((ggplot(aes(x=site, y=d.gd.mean), data=d.gd.summ) +
  geom_jitter(aes(x=site, y=d.gd.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Change in Mean Pd Score") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=8); d.gd.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.






## Change in wing score

d.wing.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.wing.score.mean = mean(d.wing.score, na.rm=T), d.wing.score.sd = sd(d.wing.score, na.rm=T)) %>%
  mutate(hi = d.wing.score.mean + d.wing.score.sd, lo = d.wing.score.mean - d.wing.score.sd)
#Dataframe with the mean +/- standard deviation of the change in wing score from early to late hibernation

d.wing.score.site.p <- ggMarginal((ggplot(aes(x=site, y=d.wing.score.mean), data=d.wing.summ) +
  geom_jitter(aes(x=site, y=d.wing.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Change in Mean Wing Score") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=8); d.wing.score.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.




## Change in mass

d.mass.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.mass.mean = mean(d.mass, na.rm=T), d.mass.sd = sd(d.mass, na.rm=T)) %>%
  mutate(hi = d.mass.mean + d.mass.sd, lo = d.mass.mean - d.mass.sd)
#Dataframe with the mean +/- standard deviation of the change in mass from early to late hibernation

d.mass.site.p <- ggMarginal((ggplot(aes(x=site, y=d.mass.mean), data=d.mass.summ) +
  geom_jitter(aes(x=site, y=d.mass, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
  scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  labs(x=NULL, y="Change in Mass (Grams") +
  theme(
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
    axis.title = element_text(size=15),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13),
    legend.position = "top"
  )), type="histogram", fill="darkgray", bins=15); d.mass.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.



## Daily change in mass 

d.mass.daily.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.mass.daily.mean = mean(d.mass.daily, na.rm=T), d.mass.daily.sd = sd(d.mass.daily, na.rm=T)) %>%
  mutate(hi = d.mass.daily.mean + d.mass.daily.sd, lo = d.mass.daily.mean - d.mass.daily.sd)
#Dataframe with the mean +/- standard deviation of the daily change in mass from early to late hibernation

d.mass.daily.site.p <- ggMarginal((ggplot(aes(x=site, y=d.mass.daily.mean), data=d.mass.daily.summ) +
                               geom_jitter(aes(x=site, y=d.mass.daily, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
                               geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
                               geom_point(size=4, color="Black", fill="White", stroke=1, shape=22)+
                               scale_color_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
                               labs(x=NULL, y="Daily Change in Mass (Grams)") +
                               theme(
                                 axis.text.y = element_text(size=13),
                                 axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                                 axis.title = element_text(size=15),
                                 legend.title = element_text(size=13),
                                 legend.text = element_text(size=13),
                                 legend.position = "top"
                               )), type="histogram", fill="darkgray", bins=15); d.mass.daily.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.



#### Other data summaries ####

#What was the breakdown of sites by sex? 
sex.summ <- dis.df %>%
  mutate(c=1) %>%
  filter(!is.na(arousal.freq.days)) %>%
  group_by(site, sex) %>%
  summarise(num.ind = sum(c))
#Always more males. Blackball and Zimmerman only had 1 female with a working transmitter each. 





# How much spatial variation exists in each site, according to point temperatures in late and early hibernation? 

dis.dat.master <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/midwest_master.csv")
dis.dat.master$date <- as.Date(dis.dat.master$date, format="%m/%d/%y")
#Fixing date column.

dis.dat.master <- dis.dat.master %>%
  filter((swab_type=="BAT" | swab_type=="FAR") & date>"2021-10-01" & date<"2022-05-01") %>%
  filter(site=="SOUTH LAKE MINE"|site=='MEAD MINE'|site=='BLACKBALL'|site=='ZIMMERMAN'|site=='GRAPHITE MINE'|site=='CP TUNNEL'|site=='ELROY SPARTA')
#Filtering the database to just contain the bats and far samples in this study. 

dis.dat.master$site <- factor(dis.dat.master$site, levels = c("SOUTH LAKE MINE","GRAPHITE MINE","CP TUNNEL", "BLACKBALL", "ZIMMERMAN", "MEAD MINE","ELROY SPARTA"))
#re-ordering factor 

early.point.temps <- dis.dat.master %>% filter(season=="hiber_earl") %>%
  select(site, season, swab_type, temp)
late.point.temps <- dis.dat.master %>% filter(season=="hiber_late") %>%
  select(site, season, swab_type, temp)
#separating the master disease database by early and late hibernation

point.temps <- rbind(early.point.temps, late.point.temps)
#Combined dataframe


## Mean point temperature (bats and fars) +/- 1 SD
mean.point.temps <- point.temps %>%
  group_by(site, season) %>%
  summarise(mean.point.temp = mean(temp, na.rm=T), mean.point.temp.sd = sd(temp, na.rm=T)) %>%
  mutate(hi.sd = mean.point.temp + mean.point.temp.sd, lo.sd = mean.point.temp - mean.point.temp.sd)
#Summary table of mean torpor bout temperature data across sites. Mean and SD range. 


mean.point.temps.p <- (ggplot(aes(x=site, y=mean.point.temp, group=season), data=mean.point.temps) +
                                         geom_point(aes(x=site, y=temp, fill=temp, shape=season), position=position_dodge(0.8),data=point.temps, size=2) +
                                         geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), position=position_dodge(0.8),width=0.5, size=0.7) +
                                         geom_point(aes(shape=season), position=position_dodge(0.8), size=4, color="Black", fill="White", stroke=1)+
                                         scale_fill_gradient(low="Blue", high="Red", name="Temperature (Celsius)")+
                                         labs(x=NULL, y="Temperature (Celsius)") +
                                    scale_y_continuous(breaks=seq(0,12,1), labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+
                                         scale_shape_manual(values=c(21,24))+
                                         theme(
                                           axis.text.y = element_text(size=13),
                                           axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                                           axis.title = element_text(size=15),
                                           legend.title = element_text(size=13),
                                           legend.text = element_text(size=13),
                                           legend.position = "top"
                                         ));mean.point.temps.p
#Plotted

#I also want to plot the average torpor bout temperature data on this plot to see how it compares to the variation that was available:
mean.point.temps.p.df <- mean.point.temps
#Making a new combined plotting dataframe. Combines mean point temps from early and late hibernation as well as mean torpor bout temps. 
mean.torpor.temp.summ.p.df <- mean.torpor.temp.summ
mean.torpor.temp.summ.p.df$season <- "Mean Torpor Bout Temperature"
colnames(mean.torpor.temp.summ.p.df) <- c("site","mean.point.temp","mean.point.temp.sd","hi.sd",'lo.sd',"season")
mean.torpor.temp.summ.p.df <- rbind(mean.torpor.temp.summ.p.df, mean.point.temps.p.df)
#This dataframe contains all the mean data. 

#Now need to combine the point temperature dataframe with the mean torpor bout temperature dataframe:
dis.df.nan.temps <- dis.df %>%
  filter(!is.na(mean.torpor.temp)) %>%
  select(site, mean.torpor.temp) %>%
  mutate(season="Mean Torpor Bout Temperature", swab_type=NA)
colnames(dis.df.nan.temps) <- c("site","temp","season","swab_type")
#This dataframe just contains each bat's mean torpor bout temperature. The season and swab_type columns are simply for merging with
#the point temperature data

point.temps.p.df <- rbind(point.temps, dis.df.nan.temps)
#Combined point temp and torpor temp dataframe

mean.torpor.temp.summ.p.df$season <- factor(mean.torpor.temp.summ.p.df$season, levels = c("hiber_earl","Mean Torpor Bout Temperature","hiber_late"))
point.temps.p.df$season <- factor(point.temps.p.df$season, levels = c("hiber_earl","Mean Torpor Bout Temperature","hiber_late"))
#re-ordering factor 

mean.point.temps.p2 <- (ggplot(aes(x=site, y=mean.point.temp, group=season), data=mean.torpor.temp.summ.p.df) +
                         geom_point(aes(x=site, y=temp, fill=temp, shape=season), position=position_dodge(0.5),data=point.temps.p.df, size=4) +
                         geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), position=position_dodge(0.5),width=0.5, size=0.7) +
                         geom_point(aes(shape=season), position=position_dodge(0.5), size=7, color="Black", fill="White", stroke=1)+
                         scale_fill_gradient(low="Blue", high="Red", name="Temperature (Celsius)")+
                         labs(x=NULL, y="Temperature (Celsius)") +
                         scale_y_continuous(breaks=seq(0,12,1), labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+
                         scale_shape_manual(values=c(21,22,24))+
                         theme(
                           axis.text.y = element_text(size=13),
                           axis.text.x = element_text(size=13, angle=70, vjust=1.05, hjust=1.05),
                           axis.title = element_text(size=15),
                           legend.title = element_text(size=13),
                           legend.text = element_text(size=13),
                           legend.position = "top"
                         ));mean.point.temps.p2




#### Behavior ~ Temperature ####

## Arousal frequency: 

m.arousal.freq.temp <- glmer(arousal.freq.days ~ mean.torpor.temp + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df); summary(m.arousal.freq.temp)
#This needs to be modeled as an exponential distribution because it is time between events. Dispersion=1 should make it exponential distribution. 
plot(dis.df$arousal.freq.days ~ dis.df$mean.torpor.temp)
#Doesn't seem like there's any relationship. 


## Mean torpor bout length: 

m.torpor.length.temp <- glmer(mean.torpor.length.days ~ mean.torpor.temp + (1|site), family=Gamma(link="log"), data=dis.df); summary(m.torpor.length.temp, dispersion=1)
#This needs to be modeled as an exponential distribution because it is time between events. Dispersion=1 should make it exponential distribution. 
plot(dis.df$mean.torpor.length.days ~ dis.df$mean.torpor.temp)
#Doesn't seem like there's any relationship. 




## Mean change in temperature bout temperature following arousal 

m.d.torpor.temp.m <- glmmTMB(mean.torpor.temp ~ mean.d.torpor.temp.weighted + (1|site), family=Gamma(link="log"), data=dis.df); summary(m.d.torpor.temp.m)
#Singularity issues when adding sex as random effect. 
plot(allEffects(m.d.torpor.temp.m))
#Can change between weighted and unweighted. 
plot( dis.df$mean.torpor.temp ~ dis.df$mean.d.torpor.temp.weighted)


m.fst <- as.data.frame(expand.grid(mean.d.torpor.temp.weighted=seq(-2,0.6,0.01), site=unique(dis.df$site), sex=unique(dis.df$sex)))
m.fst.yhat <- as.data.frame(predict(m.d.torpor.temp.m, m.fst, se.fit=T, re.form=NA, type='response'))
m.fst.yhat.fin <- cbind(m.fst, m.fst.yhat)
m.fst.yhat.fin <- m.fst.yhat.fin %>%
  group_by(mean.d.torpor.temp.weighted) %>%
  summarise(model.fit = mean(fit), model.se = mean(se.fit)) %>%
  mutate(lo.ci = model.fit-(1.96*model.se), hi.ci = model.fit + (1.96*model.se))
#This dataframe contains model predictions and 95% confidence intervals. 


p.dev.fst.temp2 <- ggplot(aes(x=mean.d.torpor.temp.weighted, y=model.fit), data=m.fst.yhat.fin) +
  geom_ribbon(aes(x=mean.d.torpor.temp.weighted, ymin=lo.ci, ymax=hi.ci), fill="gray", color="black")+
  geom_point(aes(x=mean.d.torpor.temp.weighted, y=mean.torpor.temp, fill=site, shape=sex), size=3, data=dis.df)+
  geom_line(color="black", size=1.5)+
  #scale_fill_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  scale_shape_manual(values=c(21,24), name="Sex")+
  labs(x="Weighted Mean in Change in Mean Torpor \n Bout Temperature Following Arousal (Celsius)", y="Mean Torpor Bout Temperature (Celsius)")+
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme(
    axis.text = element_text(size=13),
    axis.title = element_text(size=19),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13)
  );p.dev.fst.temp2
#Plot. 





## Mean deviation in mean torpor bout temperature from mean of first torpor bout temperature 

m.dev.first.torpor <- glmmTMB(mean.torpor.temp ~ mean.temp.dev.first.torpor.weighted+ (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df); summary(m.dev.first.torpor)
#Can change between weighted and unweighted. 

m.fst.2 <- as.data.frame(expand.grid(mean.temp.dev.first.torpor.weighted=seq(-5,0.3,0.01), site=unique(dis.df$site), sex=unique(dis.df$sex)))
m.fst.2.yhat <- as.data.frame(predict(m.dev.first.torpor , m.fst.2, se.fit=T, re.form=NA, type='response'))
m.fst.2.yhat.fin <- cbind(m.fst.2, m.fst.2.yhat)
m.fst.2.yhat.fin <- m.fst.2.yhat.fin %>%
  group_by(mean.temp.dev.first.torpor.weighted) %>%
  summarise(model.fit = mean(fit), model.se = mean(se.fit)) %>%
  mutate(lo.ci = model.fit-(1.96*model.se), hi.ci = model.fit + (1.96*model.se))
#This dataframe contains model predictions and 95% confidence intervals. 


p.dev.fst.temp2 <- ggplot(aes(x=mean.temp.dev.first.torpor.weighted, y=model.fit), data=m.fst.2.yhat.fin) +
  geom_ribbon(aes(x=mean.temp.dev.first.torpor.weighted, ymin=lo.ci, ymax=hi.ci), fill="gray", color="black")+
  geom_point(aes(x=mean.temp.dev.first.torpor.weighted, y=mean.torpor.temp, fill=site, shape=sex), size=3, data=dis.df)+
  geom_line(color="black", size=1.5)+
  #scale_fill_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  scale_shape_manual(values=c(21,24), name="Sex")+
  labs(x="Weighted Mean in Deviation from First Torpor Bout (Celsius)", y="Mean Torpor Bout Temperature (Celsius)")+
  guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme(
    axis.text = element_text(size=13),
    axis.title = element_text(size=19),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13)
  );p.dev.fst.temp2
#Plot. 


#### UV score ~ Behavior ####

dis.df$d.uv.score <- as.factor(dis.df$d.uv.score)
#Making the change in UV score a factor because treating it as numeric seems inappropriate. 

uv.m1 <- lmer(mean.torpor.temp ~ d.uv.score + (1|site) + (1|sex), data=dis.df);summary(uv.m1)
plot(allEffects(uv.m1))
plot(dis.df$mean.torpor.temp ~ dis.df$d.uv.score)
#This doesn't really seem like anything. 

#### Pd score ~ behavior ####

dis.df$d.gd.score <- as.factor(dis.df$d.gd.score)
#Making the change in gd score a factor because treating it as numeric seems inappropriate. 

gd.m1 <- lmer(mean.torpor.temp ~ d.gd.score + (1|site) + (1|sex), data=dis.df);summary(gd.m1)
plot(allEffects(gd.m1))
plot(dis.df$mean.torpor.temp ~ dis.df$d.gd.score)
#This doesn't really seem like anything. 


#### Wing score ~ behavior ####
#### Weight loss ~ behavior ####

#For modeling purposes, it's going to be easiest if the mass change metric was in positive space, so transforming it here:
dis.df$d.mass.daily.p <- -1*dis.df$d.mass.daily
#Positive values for loss in mass. 



## Arousal frequency

m.dmass1 <- glmmTMB(d.mass.daily.p ~ arousal.freq.days + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df);summary(m.dmass1)
#crossed random effects. 
plot(allEffects(m.dmass1))
plot(dis.df$d.mass.daily.p ~ dis.df$arousal.freq.days)
#Pretty clear negative negative relationship between arousal frequency and the amount of mass lost. 
#In other words, if you woke up more frequently, you lost more mass. 

m.dmass1.df <- as.data.frame(expand.grid(arousal.freq.days=seq(5,25,0.01), site=unique(dis.df$site), sex=unique(dis.df$sex)))
m.dmass1.yhat <- as.data.frame(predict(m.dmass1, m.dmass1.df, se.fit=T, re.form=NA, type='response'))
m.dmass1.yhat.fin <- cbind(m.dmass1.df, m.dmass1.yhat)
m.dmass1.yhat.fin <- m.dmass1.yhat.fin %>%
  group_by(arousal.freq.days) %>%
  summarise(model.fit = mean(fit), model.se = mean(se.fit)) %>%
  mutate(lo.ci = model.fit-(1.96*model.se), hi.ci = model.fit + (1.96*model.se))
#This dataframe contains model predictions and 95% confidence intervals. 

p.dmass1 <- ggplot(aes(x=arousal.freq.days, y=model.fit), data=m.dmass1.yhat.fin) +
  geom_ribbon(aes(x=arousal.freq.days, ymin=lo.ci, ymax=hi.ci), fill="gray", color="black")+
  geom_point(aes(x=arousal.freq.days, y=d.mass.daily.p, fill=mean.torpor.temp, shape=sex), size=3, color="black", data=dis.df)+
  geom_line(color="black", size=1.5)+
  scale_fill_gradient(low="Blue", high="Red", name="Mean Torpor Bout Temperature")+
  scale_shape_manual(values=c(21,24), name="Sex")+
  labs(x="Arousal Frequency (days)", y="Daily Mass Loss (grams)")+
  theme(
    axis.text = element_text(size=13),
    axis.title = element_text(size=19),
    plot.margin=margin(10,10,0,30),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13)
  );p.dmass1
#Plot. 


m.dmass2 <- glmer(d.mass.daily.p ~ log10(arousal.freq.days)*mean.torpor.temp + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df);summary(m.dmass2)
#Convergence issues when using non-transformed arousal frequency value. I don't know if logging that metric
#is appropriate, but it allows model to run. Suggests no association. 





## Mean torpor temperature

m.dmass3 <- glmmTMB(d.mass.daily.p ~ mean.torpor.temp + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df);summary(m.dmass3)
plot(dis.df$d.mass.daily.p ~ dis.df$mean.torpor.temp)
#Appears to be no association. 





## Mean change in mean torpor bout temperature following arousal (weighted)

m.dmass4 <- glmmTMB(d.mass.daily.p ~ mean.d.torpor.temp.weighted + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df);summary(m.dmass4)
#Nothing by itself
plot(dis.df$d.mass.daily.p ~ dis.df$mean.d.torpor.temp.weighted)

#What about interacting with temperature? 
m.dmass5 <- glmmTMB(d.mass.daily.p ~ mean.d.torpor.temp.weighted*mean.torpor.temp + (1|site) + (1|sex), family=Gamma(link="log"), data=dis.df);summary(m.dmass5)
#No interaction




#### Pathogen growth ~ behavior ####
