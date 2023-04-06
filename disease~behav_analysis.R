library(tidyverse)

dat <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/arousals_torpors_working.csv")
#This dataframe contains a list and summary of every individual's arousal and torpor events. 

dat$start.datetime <- as.POSIXct(strptime(dat$start.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
dat$end.datetime <- as.POSIXct(strptime(dat$end.datetime, "%Y-%m-%d %H:%M:%S",tz='EST'))
#Formatting date correctly

#### Building independent variables of interest #####

#First, going to start simply with the number of arousals, arousal frequency, average torpor bout length. 

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
  filter(trans_id %in% ta) %>%
  group_by(trans_id) %>%
  filter(event_num == 1 | event_num == max(event_num)-1) %>%
  mutate(sampling.duration = max(end.datetime) - min(start.datetime)) %>%
  summarise(sampling.duration=mean(sampling.duration))
#This dataframe contains the sampling durations for the at bats. 



