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

#### Identifying arousals and torpor data ####
master <- master %>%
  group_by(trans_id) %>%
  mutate(max.temp = max(temp))
#Creating a column for each bat's maximum skin temperature
master$arousal <- NA
master$arousal[master$temp >= (master$max.temp-10)] <- 1
master$arousal[is.na(master$arousal)] <- 0
#This identifies each arousal bout using the same methodology as Reeder et al. 2012. However, while it accurately identifies arousals, it does
#mis-classify much of the 'ramping up' and 'ramping down' temperatures on either end of an arousal as torpor bouts. To address this, the below 
#for-loop looks at the temperature values immediately before the arousal event, at position i. If i is greater than or equal to i-1 + 0.5, it is 
#classified as part of the arousal event. Similarly, if the temperature value immediately following the arousal event, i, is greater than or equal
#to (i+1)+0.5, it additionally is classified as part of the arousal event. The below code runs the loop twice, meaning that the two datapoints
#immediately before and after the arousal can be classified as part of that arousal if they meet the above conditions. 

x<-1
repeat{
  print(x)
  x = x+1
  for(i in 2:nrow(master)) {if(master[i,1] == master[i-1,1]) {if(master[i,10]==0 & master[i+1,10]==1 & master[i,5] >= (master[i-1,5]+1)) {master[i,10]<-1}
    else{if(master[i,10]==0 & master[i-1,10]==1 & master[i,5] >= (master[i+1,5]+1)) {master[i,10]<-1}
      else{master[i,10]<-master[i,10]}}}
    else{master[i,10] <- 0}}
  if(x==3){break}
}
#Modify the x== argument to run the loop more times. The number of times it will run will be equal to x-1

#### Construcing arousal/torpor summary dataframe ####