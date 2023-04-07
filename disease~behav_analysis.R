library(tidyverse)

dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/arousals_torpors_working.csv")
#This dataframe contains a list and summary of every individual's arousal and torpor events. 

dat$start.datetime <- as.POSIXct(strptime(dat$start.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
dat$end.datetime <- as.POSIXct(strptime(dat$end.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
#Formatting date correctly

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

ind.summ$arousal.freq.days <- as.numeric(ind.summ$sampling.duration.days/ind.summ$num.arousals)
#This is the arousal frequency or, on average, how many days separated arousal bouts. 


#Moving on to average torpor bout length.

#Again, if a torpor bout was the first, last, or both the first and last events for a bat, these torpor bouts cannot be 
#used to calculate average lengths because they could have been artificially changed by handling. 
tt.average.torpor.lengths <- dat %>%
  filter(trans_id %in% tt & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num >1 & event_num < max(event_num)) %>%
  summarise(mean.torpor.length.days = mean(event.length.days))
#This dataframe contains the mean torpor bout length for the tt bats. 

ta.average.torpor.lengths <- dat %>%
  filter(trans_id %in% ta & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num >1) %>%
  summarise(mean.torpor.length.days = mean(event.length.days))
#This dataframe contains the mean torpor bout length for the ta bats. 

at.average.torpor.lengths <- dat %>%
  filter(trans_id %in% at & behavior=="Torpor") %>%
  group_by(trans_id) %>%
  filter(event_num < max(event_num)) %>%
  summarise(mean.torpor.length.days = mean(event.length.days))
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

for(i in 2:nrow(dat)) {if(dat[i,3] == "Arousal") {dat[i,15] <- NA}
  else{if(dat[i,2] != dat[i-1,2]) {dat[i,15] <- NA}
    else{if(dat[i,2] != dat[i-2,2]) {dat[i,15] <- NA}
      else{dat[i,15] <- dat[i,7] - dat[i-2,7]}
      }
    }
  }
#This for-loop calculates the change in mean torpor bout temperature following arousal events. Of course, no change in temp
#value is given to the first torpor bout. 

mean.d.torpor.temps.weighted <- dat %>%
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
  summarise(mean.temp.dev.first.torpor.weighted = weighted.mean(x=temp.dev.first.torpor, w=event.length.days, na.rm=T))
#This dataframe contains the weighted average (again weighted by torpor bout length) of the deviations of each torpor
#bout's mean temperature from that of the first torpor bout. 
mean.temp.dev.first.torpor.unweighted <- dat %>%
  group_by(trans_id) %>%
  summarise(mean.temp.dev.first.torpor.unweighted = mean(temp.dev.first.torpor, na.rm=T))
#This is the unweighted analog of above

ind.summ$mean.temp.dev.first.torpor.weighted <- NA
ind.summ$mean.temp.dev.first.torpor.weighted <- mean.temp.dev.first.torpor.weighted$mean.temp.dev.first.torpor.weighted[match(ind.summ$trans_id, mean.temp.dev.first.torpor.weighted$trans_id)]
ind.summ$mean.temp.dev.first.torpor.weighted <- NA
ind.summ$mean.temp.dev.first.torpor.unweighted <- mean.temp.dev.first.torpor.unweighted$mean.temp.dev.first.torpor.unweighted[match(ind.summ$trans_id, mean.temp.dev.first.torpor.unweighted$trans_id)]
#Bringing the weighted and unweighted averages into the individual summary dataframe. 





#### Matching in disease data ####

dis.dat.master <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/midwest_master.csv")
#this dataframe doesn't have any infection data. That will have to be incorporated later once it's available. 

dis.dat.master$date <- as.Date(dis.dat.master$date, format="%m/%d/%y")
#Fixing date column.

dis.dat.master <- dis.dat.master %>%
  filter(!is.na(trans_id) & date>"2021-10-01" & date<"2022-05-01") %>%
  filter(site=="SOUTH LAKE MINE"|site=='MEAD MINE'|site=='BLACKBALL'|site=='ZIMMERMAN'|site=='GRAPHITE MINE'|site=='CP TUNNEL'|site=='ELROY SPARTA')
#Filtering the database to just contain the bats in this study. 

dis.early <- filter(dis.dat.master, season=="hiber_earl")
dis.late <- filter(dis.dat.master, season=="hiber_late")
#separating the master disease database by early and late hibernation. 

poo<-select(dis.early, site, trans_id, date, band)
colnames(poo) <- c("site", "trans_id", "date_early", "band_early")
poo$date_late <- dis.late$date[match(poo$trans_id, dis.late$trans_id)]
poo$band_late <- dis.late$band[match(poo$trans_id, dis.late$trans_id)]
poo<-filter(poo, band_early != band_late)
#These are inconsistencies between early and late hibernation in the band # recorded for a bat with the same trans_id. 

dis.df <- select(dis.early, trans_id, site, section, state, band, swab_id, sex, uv_orange, uv_right, gd_wall, gd_wall2, wing_score, wing_score2)
#This is going to be the dataframe with which I combine the ind.summ dataframe. 

colnames(dis.df) <- c("trans_id", "site", "section_early", "state", "band", )

