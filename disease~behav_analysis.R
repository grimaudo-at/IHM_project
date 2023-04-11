library(tidyverse)
library(lme4)
library(effects)
library(emmeans)

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
ind.summ$mean.temp.dev.first.torpor.unweighted <- NA
ind.summ$mean.temp.dev.first.torpor.unweighted <- mean.temp.dev.first.torpor.unweighted$mean.temp.dev.first.torpor.unweighted[match(ind.summ$trans_id, mean.temp.dev.first.torpor.unweighted$trans_id)]
#Bringing the weighted and unweighted averages into the individual summary dataframe. 





#### Matching in disease data ####

dis.dat.master <- read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/midwest_master.csv")
trans_meta<- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv")
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

dis.df <- select(dis.early, trans_id, trans_model, trans_weight_g, site, section, state, band, swab_id, sex, mass, uv_orange, uv_right, gd_wall, gd_wall2, wing_score, wing_score2)
#This is going to be the dataframe with which I combine the ind.summ dataframe. 

colnames(dis.df) <- c("trans_id", "trans.model", "trans.weight.g", "site", "section.early", "state", "band", "swab.id.early", "sex", "mass.early","uv.orange.early", 
                      "uv.right.early", "gd.wall.early","gd.wall2.early", "wing.score.early","wing.score2.early")
#Re-naming columns so that their early disease information is identifiable.

dis.df$mass.late <- dis.late$mass[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$swab.id.late <- dis.late$swab_id[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$uv.orange.late <- dis.late$uv_orange[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$uv.right.late <- dis.late$uv_right[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$gd.wall.late <- dis.late$gd_wall[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$gd.wall2.late <- dis.late$gd_wall2[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$wing.score.late <- dis.late$wing_score[match(dis.df$trans_id, dis.late$trans_id)]
dis.df$wing.score2.late <- dis.late$wing_score2[match(dis.df$trans_id, dis.late$trans_id)]
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

## Can now match in all the transmitter summary data for each individual: 
dis.df <- left_join(dis.df, ind.summ[,2:11], by="trans_id")
#Merged. 

#### Variation across sites in behavior metrics ####

dis.df$site <- factor(dis.df$site, levels = c("SOUTH LAKE MINE","GRAPHITE MINE","CP TUNNEL", "BLACKBALL", "ZIMMERMAN", "MEAD MINE","ELROY SPARTA"))
#Make sure site is a factor. I have ordered the levels of this factor from coldest mean torpor bout temp to warmest. 


## Arousal frequency

ar.freq.summ <- dis.df %>%
  group_by(site) %>%
  summarise(ar.freq.mean = mean(arousal.freq.days, na.rm=T), ar.freq.sd = sd(arousal.freq.days, na.rm=T)) %>%
  mutate(hi.sd = ar.freq.mean + ar.freq.sd, lo.sd = ar.freq.mean - ar.freq.sd)
#Summary table of arousal frequency data across sites. Mean and SD range. 

ar.freq.site.p <- ggplot(aes(x=site, y=ar.freq.mean), data=ar.freq.summ) +
  geom_jitter(aes(x=site, y=arousal.freq.days, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); ar.freq.site.p
#Plotted



## Mean torpor length

mean.torpor.length.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.torpor.length.mean = mean(mean.torpor.length.days, na.rm=T), mean.torpor.length.sd = sd(mean.torpor.length.days, na.rm=T)) %>%
  mutate(hi.sd = mean.torpor.length.mean + mean.torpor.length.sd, lo.sd = mean.torpor.length.mean - mean.torpor.length.sd)
#Summary table of mean torpor bout length data across sites. Mean and SD range. 

mean.torpor.length.site.p <- ggplot(aes(x=site, y=mean.torpor.length.mean), data=mean.torpor.length.summ) +
  geom_jitter(aes(x=site, y=mean.torpor.length.days, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.torpor.length.site.p
#Plotted




## Mean torpor temp

mean.torpor.temp.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.torpor.temp.mean = mean(mean.torpor.temp, na.rm=T), mean.torpor.temp.sd = sd(mean.torpor.temp, na.rm=T)) %>%
  mutate(hi.sd = mean.torpor.temp.mean + mean.torpor.temp.sd, lo.sd = mean.torpor.temp.mean - mean.torpor.temp.sd)
#Summary table of mean torpor bout temperature data across sites. Mean and SD range. 

mean.torpor.temp.site.p <- ggplot(aes(x=site, y=mean.torpor.temp.mean), data=mean.torpor.temp.summ) +
  geom_jitter(aes(x=site, y=mean.torpor.temp, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.torpor.temp.site.p
#Plotted




## Mean change in mean torpor bout temperature

#Weighted first:
mean.d.torpor.temp.w.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.d.torpor.temp.w.mean = mean(mean.d.torpor.temp.weighted, na.rm=T), mean.d.torpor.temp.w.sd = sd(mean.d.torpor.temp.weighted, na.rm=T)) %>%
  mutate(hi.sd = mean.d.torpor.temp.w.mean + mean.d.torpor.temp.w.sd, lo.sd = mean.d.torpor.temp.w.mean - mean.d.torpor.temp.w.sd)
#Summary table of change in torpor bout temperature data across sites. Mean and SD range. 

mean.d.torpor.temp.w.site.p <- ggplot(aes(x=site, y=mean.d.torpor.temp.w.mean), data=mean.d.torpor.temp.w.summ) +
  geom_jitter(aes(x=site, y=mean.d.torpor.temp.weighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.d.torpor.temp.w.site.p
#Plotted

#Now unweighted:
mean.d.torpor.temp.uw.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.d.torpor.temp.uw.mean = mean(mean.d.torpor.temp.unweighted, na.rm=T), mean.d.torpor.temp.uw.sd = sd(mean.d.torpor.temp.unweighted, na.rm=T)) %>%
  mutate(hi.sd = mean.d.torpor.temp.uw.mean + mean.d.torpor.temp.uw.sd, lo.sd = mean.d.torpor.temp.uw.mean - mean.d.torpor.temp.uw.sd)
#Summary table of change in torpor bout temperature data across sites. Mean and SD range. 

mean.d.torpor.temp.uw.site.p <- ggplot(aes(x=site, y=mean.d.torpor.temp.uw.mean), data=mean.d.torpor.temp.uw.summ) +
  geom_jitter(aes(x=site, y=mean.d.torpor.temp.unweighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.d.torpor.temp.uw.site.p
#Plotted






## Mean deviation in mean torpor bout temperature from first torpor bout

#Weighted first:
mean.temp.dev.w.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.w.mean = mean(mean.temp.dev.first.torpor.weighted, na.rm=T), mean.temp.dev.w.sd = sd(mean.temp.dev.first.torpor.weighted, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.w.mean + mean.temp.dev.w.sd, lo.sd = mean.temp.dev.w.mean - mean.temp.dev.w.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.w.site.p <- ggplot(aes(x=site, y=mean.temp.dev.w.mean), data=mean.temp.dev.w.summ) +
  geom_jitter(aes(x=site, y=mean.temp.dev.first.torpor.weighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.temp.dev.w.site.p
#Plotted

#Now unweighted:
mean.temp.dev.uw.summ <- dis.df %>%
  group_by(site) %>%
  summarise(mean.temp.dev.uw.mean = mean(mean.temp.dev.first.torpor.unweighted, na.rm=T), mean.temp.dev.uw.sd = sd(mean.temp.dev.first.torpor.unweighted, na.rm=T)) %>%
  mutate(hi.sd = mean.temp.dev.uw.mean + mean.temp.dev.uw.sd, lo.sd = mean.temp.dev.uw.mean - mean.temp.dev.uw.sd)
#Summary table of deviation in torpor bout temperature from first torpor bout data across sites. Mean and SD range. 

mean.temp.dev.uw.site.p <- ggplot(aes(x=site, y=mean.temp.dev.uw.mean), data=mean.temp.dev.uw.summ) +
  geom_jitter(aes(x=site, y=mean.temp.dev.first.torpor.unweighted, color=mean.torpor.temp), data=dis.df, size=2, height=0, width=0.2) +
  geom_errorbar(aes(ymin=lo.sd, ymax=hi.sd), width=0.2, size=0.7) +
  geom_point(color='Black', size=4) +
  scale_color_gradient(low="Blue", high="Red"); mean.temp.dev.uw.site.p
#Plotted



#### Variation across sites in disease metrics ####


## Change in UV score:

d.uv.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.uv.mean = mean(d.uv.score, na.rm=T), d.uv.sd = sd(d.uv.score, na.rm=T)) %>%
  mutate(hi = d.uv.mean + d.uv.sd, lo = d.uv.mean - d.uv.sd)
#Dataframe with the mean +/- standard deviation of the change in UV score from early to late hibernation

d.uv.site.p <- ggplot(aes(x=site, y=d.uv.mean), data=d.uv.summ) +
  geom_jitter(aes(x=site, y=d.uv.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Blue")+
  scale_color_gradient(low="Blue", high="Red"); d.uv.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.



## Change in pd score

d.gd.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.gd.mean = mean(d.gd.score, na.rm=T), d.gd.sd = sd(d.gd.score, na.rm=T)) %>%
  mutate(hi = d.gd.mean + d.gd.sd, lo = d.gd.mean - d.gd.sd)
#Dataframe with the mean +/- standard deviation of the change in gd score from early to late hibernation

d.gd.site.p <- ggplot(aes(x=site, y=d.gd.mean), data=d.gd.summ) +
  geom_jitter(aes(x=site, y=d.gd.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Blue")+
  scale_color_gradient(low="Blue", high="Red"); d.gd.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.




## Change in wing score

d.wing.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.wing.score.mean = mean(d.wing.score, na.rm=T), d.wing.score.sd = sd(d.wing.score, na.rm=T)) %>%
  mutate(hi = d.wing.score.mean + d.wing.score.sd, lo = d.wing.score.mean - d.wing.score.sd)
#Dataframe with the mean +/- standard deviation of the change in wing score from early to late hibernation

d.wing.score.site.p <- ggplot(aes(x=site, y=d.wing.score.mean), data=d.wing.summ) +
  geom_jitter(aes(x=site, y=d.wing.score, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Blue")+
  scale_color_gradient(low="Blue", high="Red"); d.wing.score.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.




## Change in mass

d.mass.summ <- dis.df %>%
  group_by(site) %>%
  summarise(d.mass.mean = mean(d.mass, na.rm=T), d.mass.sd = sd(d.mass, na.rm=T)) %>%
  mutate(hi = d.mass.mean + d.mass.sd, lo = d.mass.mean - d.mass.sd)
#Dataframe with the mean +/- standard deviation of the change in mass from early to late hibernation

d.mass.site.p <- ggplot(aes(x=site, y=d.mass.mean), data=d.mass.summ) +
  geom_jitter(aes(x=site, y=d.mass, color=mean.torpor.temp), size=2, width=0.2, height=0.01, data=dis.df) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2, size=0.7) +
  geom_point(size=4, color="Blue")+
  scale_color_gradient(low="Blue", high="Red"); d.mass.site.p
#Plotting variation across sites. There are gray points here because these are bats that were re-captured but their transmitter wasn't working.


#### Disease ~ arousal frequency ####


#### Disease ~ mean torpor bout temp ####
#### Disease ~ mean change in torpor bout temp ####
#### Disease ~ mean deviation from temp of first torpor bout ####
