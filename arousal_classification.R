library(tidyverse)
master <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv")
#Calibrated transmitter temperature master dataset 

master$site<-as.factor(master$site)
master$trans_id<-as.factor(master$trans_id)
master$date <- as.Date(master$date, format="%Y-%m-%d")
master$time <- format(strptime(master$time, "%H:%M:%S"), "%H:%M:%S")
master$serial_num <- as.factor(master$serial_num)
master$logger_model <- as.factor(master$logger_model)
master$datetime <- as.POSIXct(strptime(paste(master$date, master$time), "%Y-%m-%d %H:%M:%S",tz='EST'))
#Formatting data correctly

master<-filter(master, temp < 50)
#This removes several instances where the loggers hay-wired and read extremely high values. This only happened for two loggers, #104 and 108, both at
#CP tunnel. These high values result in mis-classifications of arousals. 

#### Identifying arousals and torpor data ####
master <- master %>%
  group_by(trans_id) %>%
  mutate(max.temp = max(temp))
#Creating a column for each bat's maximum skin temperature
master$behavior <- NA
master$behavior[master$temp >= (master$max.temp-10)] <- "Arousal"
master$behavior[is.na(master$behavior)] <- "Torpor"
#This identifies each arousal bout using the same methodology as Reeder et al. 2012. However, while it accurately identifies arousals, it does
#mis-classify much of the 'ramping up' and 'ramping down' temperatures on either end of an arousal as torpor bouts. To address this, the below 
#for-loop looks at the temperature values immediately before the arousal event, at position i. If i is greater than or equal to i-1 + 0.5, it is 
#classified as part of the arousal event. Similarly, if the temperature value immediately following the arousal event, i, is greater than or equal
#to (i+1)+0.5, it additionally is classified as part of the arousal event. The below code runs the loop twice, meaning that the two datapoints
#immediately before and after the arousal can be classified as part of that arousal if they meet the above conditions. 

x<-1
start.time<-Sys.time()
repeat{
  print(x)
  x = x+1
  for(i in 2:nrow(master)) {if(master[i,1] == master[i-1,1]) {if(master[i,10]=="Torpor" & master[i+1,10]=="Arousal" & master[i,5] >= (master[i-1,5]+1)) {master[i,10]<-"Arousal"}
    else{if(master[i,10]=="Torpor" & master[i-1,10]=="Arousal" & master[i,5] >= (master[i+1,5]+1)) {master[i,10]<-"Arousal"}
      else{master[i,10]<-master[i,10]}}}
    else{master[i,10] <- "Torpor"}}
  if(x==6){break}
}
end.time<-Sys.time()
time.taken<-end.time-start.time;time.taken
#This chunk takes about 1.2 hours to run. Modify the x== argument to run the loop more times. The number of times it will run will be equal to x-1

master$event_num <- NA
master[1,11] <- 1
for(i in 2:nrow(master)) {if(master[i,2] == master[i-1,2]) {if(master[i,10] == master[i-1,10]) {master[i,11] <- master[i-1,11]}
  else{master[i,11]<-(master[i-1,11]+1)}}
  else{master[i,11] <- 1}}
#This above for-loop counts events. Takes a long time to run. Each torpor or arousal is an "event", and its "event number" is its position in the sequence of torpors/arousals. 

#write.csv(master, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv", row.names = F)

p <- list()
logg <- unique(master$trans_id)
for(i in 1:length(logg)){
p[[i]] <- list()
dat <- subset(master, trans_id==logg[i])
p[[i]][[1]] <- ggplot(dat, aes(datetime,temp, color=behavior)) + 
geom_line(color="black")+
geom_point()+
scale_y_continuous(limits=c(-2,30))+
labs(x="Date", y=expression("Temperature " (degree*C)~" "))+
scale_color_manual(values=c("red","blue"))+
ggtitle(paste(dat$trans_id))+
theme(panel.background = element_blank(),
     legend.position = "top",
    axis.text = element_text(size=20, color='grey16'), 
   axis.title = element_text(size=25, color="grey16"),
  axis.ticks = element_line(size=0.8, color="grey57"),
 panel.border = element_rect(colour = "grey57", fill=NA, size=0.8))
}
#list of ggplots to plot for quality control. 



#### Construcing arousal/torpor summary dataframe ####

tr.sum <- master %>%
  group_by(site, trans_id, behavior, event_num) %>%
  summarise(start.datetime = min(datetime), end.datetime=max(datetime), mean.temp = mean(temp), min.temp = min(temp), max.temp = max(temp)) %>%
  mutate(event.length = as.numeric(difftime(end.datetime, start.datetime, units="days"))) %>%
  filter(trans_id != "104" & trans_id != "108")

tr.sum$arousal <-NA
tr.sum$arousal[tr.sum$behavior=="Arousal"]<-1

tr.sum.freq <- tr.sum %>%
  group_by(trans_id) %>%
  summarise(start.time = min(start.datetime), end.time = max(end.datetime)) %>%
  mutate(sampling.duration = end.time - start.time)
tr.sum.freq$sampling.duration <- as.numeric(tr.sum.freq$sampling.duration)

poo <- tr.sum %>%
  group_by(site, trans_id) %>%
  summarise(n.ar = sum(arousal, na.rm=T)) %>%
  filter(trans_id != "104" & trans_id != "108")
poo$sampling.duration <- tr.sum.freq$sampling.duration[match(poo$trans_id, tr.sum.freq$trans_id)]
poo$ar.freq <- poo$sampling.duration/poo$n.ar

tr.sum2<- tr.sum %>%
  group_by(trans_id) %>%
  filter(event_num != min(event_num) & event_num != max(event_num))

p<-ggplot(aes(x=site, y=event.length), data=tr.sum2[tr.sum2$behavior=="Torpor",]) +
  geom_boxplot()+
  geom_jitter(aes(color=mean.temp), width=0.1)+
  labs(x="Site",y="Torpor bout length (days)")+
  scale_color_gradient(low="blue", high="red");p

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/torpor_length~site.PNG",p,scale=3,width=8,height=6 ,units="cm",dpi=600)

p2<-ggplot(aes(x=site, y=n.ar), data=poo) +
  geom_boxplot()+
  geom_jitter(width=0.1)+
  labs(x="Site", y="Number of arousals");p2

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/n_arousals~site.PNG",p2,scale=3,width=8,height=6 ,units="cm",dpi=600)

p3 <- ggplot(aes(x=site, y=ar.freq), data=poo) +
  geom_boxplot()+
  geom_jitter(width=0.1)+
  labs(x="Site", y="Arousal frequency (days between arousals)");p3

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/arousal_freq~site.PNG",p3,scale=3,width=8,height=6 ,units="cm",dpi=600)

p4<-ggplot(aes(x=start.datetime, y=event.length, color=mean.temp), data=tr.sum2[tr.sum2$behavior=="Torpor",]) +
  geom_point()+
  facet_wrap(~site)+
  labs(x="Date", y="Torpor bout length (days)")+
  scale_color_gradient(low="blue", high="red");p4

#ggsave(file="/Users/alexg8/Dropbox/Grimaudo_WNS_Project/figs/IHM Project/Exploratory/torpor_length~date*site.PNG",p4,scale=3,width=8,height=6 ,units="cm",dpi=600)
