library(tidyverse)
master <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_raw_working.csv")
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
#This chunk takes about 4 hours to run. Modify the x== argument to run the loop more times. The number of times it will run will be equal to x-1

master$event_num <- NA
master[1,11] <- 1
for(i in 2:nrow(master)) {if(master[i,2] == master[i-1,2]) {if(master[i,10] == master[i-1,10]) {master[i,11] <- master[i-1,11]}
  else{master[i,11]<-(master[i-1,11]+1)}}
  else{master[i,11] <- 1}}
#This above for-loop counts events. Takes a long time to run. Each torpor or arousal is an "event", and its "event number" is its position in the sequence of torpors/arousals. 

#write.csv(master, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_working.csv", row.names = F)


#### Summarizing event data ####

master <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_raw_working.csv")

master<-filter(master, temp < 50)
#This removes several instances where the loggers hay-wired and read extremely high values. This only happened for two loggers, #104 and 108, both at
#CP tunnel. These high values result in mis-classifications of arousals. 

master$site<-as.factor(master$site)
master$trans_id<-as.factor(master$trans_id)
master$date <- as.Date(master$date, format="%Y-%m-%d")
master$time <- format(strptime(master$time, "%H:%M:%S"), "%H:%M:%S")
master$serial_num <- as.factor(master$serial_num)
master$logger_model <- as.factor(master$logger_model)
master$datetime <- as.POSIXct(strptime(paste(master$date, master$time), "%Y-%m-%d %H:%M:%S",tz='EST'))

mean.torpor.temps <- master %>%
  group_by(site, trans_id, behavior) %>%
  summarise(mean.torpor.temp = mean(temp)) %>%
  filter(behavior != "Arousal")
#This dataframe contains each bat's average torpor temperature. 
  
tr.sum <- master %>%
  group_by(site, trans_id, behavior, event_num) %>%
  summarise(start.datetime = min(datetime), end.datetime=max(datetime), mean.temp = mean(temp), median.temp = median(temp), min.temp = min(temp), max.temp = max(temp), sd.temp = sd(temp)) %>%
  mutate(event.length = as.numeric(difftime(end.datetime, start.datetime, units="hours")), temp_range = max.temp - min.temp) %>%
  group_by(trans_id) %>%
  arrange(event_num, .by_group = T)
#This table is a summary of each of each individual's arousal and torpor bouts. 

tr.sum$mean.torpor.temp <- mean.torpor.temps$mean.torpor.temp[match(tr.sum$trans_id, mean.torpor.temps$trans_id)]
#Matching in the mean torpor temperature data. 


#write.csv(tr.sum, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/arousals_torpors_working.csv", row.names=F)


#### Attempting to classify movements and locations #####
## THE BELOW CHUNK OF CODE ASSIGNS 'LOCATIONS' TO TORPOR EVENTS, ATTEMPTING TO CLASSIFY MOVEMENTS TO DIFFERENT SECTIONS. THIS IS AN IMPERFECT PROCESS AND CERTAINLY
## DOES NOT CORRECTLY IDENTIFY ALL TRUE MOVEMENTS AND MIS-CLASSIFIES OTHERS. 
tr.sum$location <- NA
tr.sum[1,14] <- 1
tr.sum[1,14] <- 1
tr.sum$sd.temp[is.na(tr.sum$sd.temp)]<-0
#This column is going to keep track of how many different locations an individual has used and when. First value needs to be added to run loop

for(i in 3:nrow(tr.sum)) {if(tr.sum[i,2] != tr.sum[i-1,2]) {if(tr.sum[i,3] == "Torpor") {tr.sum[i,14] <- 1} else{tr.sum[i,14] <- NA}}
  else{if(tr.sum[i,3] != "Torpor") {tr.sum[i,14] <- NA} 
    else{if(tr.sum[i,2] != tr.sum[i-2,2]) {tr.sum[i,14] <- 1}
      else{if(tr.sum[i,8] > (tr.sum[i-2,8] + tr.sum[i-2,11])) {tr.sum[i,14] <- (tr.sum[i-2,14] + 1)}
        else{if(tr.sum[i,8] < (tr.sum[i-2,8] - tr.sum[i-2,11])) {tr.sum[i,14] <- (tr.sum[i-2,14]+1)}
          else{tr.sum[i,14] <- tr.sum[i-2,14]}}}}}}
#This loop identifies the locations in which all bouts of torpor occurred, if they were new from the previous torpor bout. If they weren't new 
#locations, then they were assigned the same location value as the previous torpor bout. New locations were identified if the average temperature
#during the associated torpor bout was outside of 1 standard deviation of the average temperature of the previous torpor bout. 

tr.sum$uniq <- paste(tr.sum$trans_id, tr.sum$event_num)
master$uniq <- paste(master$trans_id, master$event_num)
#Unique columns for matching in the location information
master$location <- tr.sum$location[match(master$uniq, tr.sum$uniq)]
#Bringin in location data to master df
master$location2[master$location%%2 == 1]<-"a"
master$location2[master$location%%2 != 1]<-"b"
master$location2[is.na(master$location2)]<-"c"
#The above is just meant for illustration purposes. It is an odd/even identifier used to color by to illustrate movements. 


p <- list()
logg <- unique(master$trans_id)
for(i in 1:length(logg)){
  p[[i]] <- list()
  dat <- subset(master, trans_id==logg[i])
  p[[i]][[1]] <- ggplot(dat, aes(datetime,temp, color=location2)) + 
    geom_line(color="black")+
    geom_point(alpha=0.9)+
    scale_y_continuous(limits=c(-2,30))+
    labs(x="Date", y=expression("Temperature " (degree*C)~" "))+
    scale_color_manual(values=c("a"="red","b"="blue","c"="black"))+
    ggtitle(paste(dat$trans_id))+
    theme(panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(size=20, color='grey16'), 
          axis.title = element_text(size=25, color="grey16"),
          axis.ticks = element_line(size=0.8, color="grey57"),
          panel.border = element_rect(colour = "grey57", fill=NA, size=0.8))
}
#list of ggplots to plot for quality control. 


