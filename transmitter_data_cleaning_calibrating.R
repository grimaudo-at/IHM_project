library(tidyverse)
library(berryFunctions)
library(chron)
library(ggplot2)


#### Cleaning and binding ####
dfs <- list.files(path = c("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Transmitter data"), pattern = "*.csv", full.names = T) %>%
  lapply(function(x) {read.csv(x, header=F)})
#List of all transmitter dataframes

l <- lapply(dfs, function(x) {if(str_detect(x[1,1], 'DS2422')=="TRUE") {"DS1922L"}
  else{if(str_detect(x[1,1], 'DS1922L')=="TRUE") {"DS1922L"}
    else{"DS1921G"}}})
#This is a list of each logger's model. Relative position in the list should match that of dfs. 

l2 <- lapply(dfs, function(x) {str_sub(x[2,1],-8,-1)})
#This is a list of the last 8 digits of each logger's serial number. Again, relative position in list should match that of dfs. 

dfs2 <- lapply(dfs, function(x) {if(str_detect(x[1,1], 'DS2422')=="TRUE") {x<-x[-c(1:21),]}
  else{if(str_detect(x[1,1], 'DS1922L')=="TRUE") {x<-x[-c(1:21),]}
    else{x<-x[-c(1:16),]}}})
#This removes metadata rows from each dataframe. Because DS1922L and DS1921G models record different numbers of rows of metadata, differing numbers
#of rows need to be filtered away. 

dfs2<-lapply(dfs2, function(x) {as.data.frame(x)})
#Makes all elements in list into dataframes again. 

for(i in seq_along(dfs2)) {
  dfs2[[i]]$thing <- as.factor(rep(c("datetime","unit","value")))
} 
#Adds a column to each dataframe that will then become column headers. 

for(i in seq_along(dfs2)) {
  dfs2[[i]]$g <- rep(c(seq(1:(nrow(dfs2[[i]])/3))), each=3)
}
#Making repeating values to identify unique observations
pivot <- function(t) {pivot_wider(t, names_from=thing, values_from = x, id_cols=g)}
dfs2<-lapply(dfs2, pivot)
#Coercing dataframes to wider format
for(i in seq_along(dfs2)) {
  dfs2[[i]] <- select(dfs2[[i]], datetime, value)
}
#Removing unnecessary columns
dfs2<-lapply(dfs2, function(x) {as.data.frame(x)})
#Makes all elements in list into dataframes again. 

for(i in seq_along(dfs2)) {
  dfs2[[i]]$serial_num <- l2[[i]]
}
#Bringing in serial # data
for(i in seq_along(dfs2)) {
  dfs2[[i]]$logger_model <- l[[i]]
}
#Bringing in logger model data

logger_meta <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv")
#Reading in transmitter metadata database

logger_meta <- filter(logger_meta, functional=="Y")
#Filtering away all un-recovered or non-functional transmitters. THIS IS IMPORTANT TO DO BEFORE MATCHING IN METADATA BECAUSE THERE MAY BE INNACURATE
#SERIAL NUMBER VALUES FOR UN-RETRIEVED OR NON-FUNCTIONAL TRANSMITTERS. 

for(i in seq_along(dfs2)) {
  dfs2[[i]]$site <- logger_meta$site[match(dfs2[[i]]$serial_num, logger_meta$serial)]
}
#Matching in site data

for(i in seq_along(dfs2)) {
  dfs2[[i]]$trans_id <- logger_meta$id[match(dfs2[[i]]$serial_num, logger_meta$serial)]
}
#Matching in transmitter ID data

for(i in seq_along(dfs2)) {
  dfs2[[i]]$cal.group <- logger_meta$calibration_group[match(dfs2[[i]]$serial_num, logger_meta$serial)]
}
#Matching in calibration group data

master <- do.call("rbind", dfs2)
#Dataframe containing all transmitter data. 

f.tim <- filter(master, serial_num == "7BB15441")
#This logger's time format is different than the rest, so I am going to remove it and then loop it back in again after fixing time data for rest. 
master<- filter(master, serial_num != "7BB15441")
master$time2<-NA
for(i in 1:nrow(master)) {if(str_detect(master[i,1], 'AM')=="TRUE") {master[i,8] <- "AM"}
  else{master[i,8] <- "PM"}}
master<-separate(master, datetime, into=c("date","time"),sep=" ")
master$time <- paste(master$time, master$time2, sep=" ")
master<-subset(master, select=-c(time2))
master$time <- format(strptime(master$time, "%I:%M:%S %p"), "%H:%M:%S")
master$date <- as.Date(master$date, format="%m/%d/%y")
#Above chunk of code is correcting format of dates and times. 

f.tim <- separate(f.tim, datetime, into=c("date","time"), sep=" ")
f.tim$time <- format(strptime(f.tim$time, "%H:%M"), "%H:%M:%S")
f.tim$date <- as.Date(f.tim$date, format="%m/%d/%Y")
#Above chunk of code is correcting format of dates and times for previously removed logger, 7BB15441

master<-rbind(master, f.tim)
#re-binding master and 7BB15441

master$value <- as.numeric(master$value)
#Converting temperature data to numeric

#### Calibrating ####

t.coup <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_calibration_thermocouples/all_T_values.csv")
#Data from thermocouples
t.coup$start_date <- as.Date(t.coup$start_date, format="%m/%d/%y")
t.coup$end_date <- as.Date(t.coup$end_date, format="%m/%d/%y")
t.coup$start_time <- format(strptime(t.coup$start_time, "%H:%M"), "%H:%M:%S")
t.coup$end_time <- format(strptime(t.coup$end_time, "%H:%M"), "%H:%M:%S")

g1 <- filter(master, cal.group=="1")
g2 <- filter(master, cal.group=="2")
g3 <- filter(master, cal.group=="3")
#These three dataframes contain the transmitter data from each of the three calibration groups. 

g1$T <- NA
g2$T <- NA
g3$T <- NA
#Columns in which I will identify data recorded during each of the 3 calibration trials. 

g1$T[g1$date==t.coup[1,3] & g1$time >= t.coup[1,5] & g1$time <= t.coup[1,6]] <- "T1"
g1$T[g1$date==t.coup[2,3] & g1$time >= t.coup[2,5] & g1$time <= t.coup[2,6]] <- "T2"
g1$T[g1$date==t.coup[3,3] & g1$time >= t.coup[3,5] & g1$time <= t.coup[3,6]] <- "T3"
g2$T[g2$date==t.coup[4,3] & g2$time >= t.coup[4,5] & g2$time <= t.coup[4,6]] <- "T1"
g2$T[g2$date==t.coup[5,4] & g2$time >= "00:00:00" & g2$time <= t.coup[5,6]] <- "T2"
g2$T[g2$date==t.coup[6,3] & g2$time >= t.coup[6,5] & g2$time <= t.coup[6,6]] <- "T3"
g3$T[g3$date==t.coup[7,3] & g3$time >= t.coup[7,5] & g3$time <= t.coup[7,6]] <- "T1"
#Identifying transmitter data recorded by each logger in each calibration trial. 

g1.a <- g1 %>%
  filter(!is.na(T)) %>%
  group_by(serial_num, cal.group, T) %>%
  summarise(Dc = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
g2.a <- g2 %>%
  filter(!is.na(T)) %>%
  group_by(serial_num, cal.group, T) %>%
  summarise(Dc = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
g3.a <- g3 %>%
  filter(!is.na(T)) %>%
  group_by(serial_num, cal.group, T) %>%
  summarise(Dc = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
#Summarizing each transmitter's D values for each calibration period. Uniq is a column for matching
t.coup$uniq <- paste(t.coup$group, t.coup$T, sep="_")
#Identical column for matching in the thermocouple dataframe

g1.a$Tc <- t.coup$value[match(g1.a$uniq, t.coup$uniq)]
g2.a$Tc <- t.coup$value[match(g2.a$uniq, t.coup$uniq)]
g3.a$Tc <- t.coup$value[match(g3.a$uniq, t.coup$uniq)]
#Bringing in T data from thermocouples. 

#In order to run these data through the calibration equations, it would be useful to have this data in wide format with each logger being a row. 
#I have enough calibration data to succesfully calibrate transmitters in groups 1 and 2, but not 3. Group 3's calibration curve will have to be 
#estimated (see next section).

g12.a<-rbind(g1.a, g2.a)
#Binding together those loggers that can be calibrated

g12.a$D<-NA
g12.a$D[g12.a$T=="T1"] <- "D1"
g12.a$D[g12.a$T=="T2"] <- "D2"
g12.a$D[g12.a$T=="T3"] <- "D3"
#Adding in data on D timepoints

g12.a.twide <- pivot_wider(g12.a, names_from=c(T), values_from=c(Tc), id_cols = serial_num)
g12.a.dwide <- pivot_wider(g12.a, names_from=c(D), values_from=c(Dc), id_cols = serial_num)
g12.a.wide <- merge(g12.a.twide, g12.a.dwide, by="serial_num")
g12.a.wide$cal.group <- g12.a$cal.group[match(g12.a.wide$serial_num, g12.a$serial_num)]
#Now it's in wide format.

#### Group 3 calibration ####

#First, I will calibrate loggers in groups 1 and 2 (T and D values stored in dataframe g12.a.wide) because they have data for all 3 calibration phases.
#Then, below, I will use information from these two groups' calibration curve to estimate the calibration for those in group 3, for which only the 
#first calibration phase was completed. 

#Calibrating group 1 and 2:

g12.a.wide$a <- (g12.a.wide$D1 - g12.a.wide$D3) + ((g12.a.wide$D2 - g12.a.wide$D1)*(g12.a.wide$T3 - g12.a.wide$T1))/(g12.a.wide$T2 - g12.a.wide$T1)
g12.a.wide$b <- g12.a.wide$T3^2 - g12.a.wide$T2*g12.a.wide$T3 + g12.a.wide$T2*g12.a.wide$T1 - g12.a.wide$T3*g12.a.wide$T1
g12.a.wide$c <- -(g12.a.wide$a/g12.a.wide$b)
g12.a.wide$B <- (g12.a.wide$D2 - g12.a.wide$D1 - g12.a.wide$c*((g12.a.wide$T2^2)-(g12.a.wide$T1^2)))/(g12.a.wide$T2-g12.a.wide$T1)
g12.a.wide$A <- g12.a.wide$D1 - g12.a.wide$B*g12.a.wide$T1 - g12.a.wide$c*(g12.a.wide$T1^2)
#Using the above values, we should be able to correct all temp data from groups 1 and 2

#Now, I need to estimate D2 and D3 values for all group 3 loggers based on the difference between their T1 and D1 values. I will use the below
#regression on the deviation values for D1, D2, and D3 in groups 1 and 2. Basically, the idea here is to use information on groups 1 and 2 and 
#how their transmitters deviated from the temps recorded by the thermocouple to estimate the same deviation for group 3. I have this deviation data for
#group 3 between their T1 and D1 values, so I will use group 1 and 2's deviation data to see if we can predict deviation at T2 and T3 based on their 
#deviation at T1. 

#First, need to calculate deviation data:
g12.a.wide$d.D1 <- g12.a.wide$T1 - g12.a.wide$D1
g12.a.wide$d.D2 <- g12.a.wide$T2 - g12.a.wide$D2
g12.a.wide$d.D3 <- g12.a.wide$T3 - g12.a.wide$D3

dev.plot1 <- ggplot(aes(x=d.D1, y=d.D2), data=g12.a.wide)+
  geom_point()+
  geom_smooth(method="lm", color="blue")+
  geom_abline(a=0, b=1);dev.plot1 
#blue line is a regression line, black is a 1:1 line. 
dev.plot2 <- ggplot(aes(x=d.D1, y=d.D3), data=g12.a.wide)+
  geom_point()+
  geom_smooth(method="lm", color="blue")+
  geom_abline(a=0, b=1);dev.plot2
#Same plot as above, but now comparing d.D1 to d.D2.
#Major take-away: deviation between T1 and D1 has some predictive power on future deviations. However, it is not a 1:1 line, but rather slightly negative.
#The residuals of this regression can actually be quite large (sometimes an entire degree), meaning that just deviation data alone is now powerfully
#predictive. However, for our purposes, it will have to suffice.

d12.m <- lm (d.D2 ~ d.D1, data=g12.a.wide);summary(d12.m) #slope of line estimated to be 0.66
d13.m <- lm (d.D3 ~ d.D1, data=g12.a.wide);summary(d13.m) #slope of line estimated to be 0.73
#Regressions I'll use to estimate group 3's D2 and D3 values

g3.a$D<-NA
g3.a$D[g3.a$T=="T1"] <- "D1"
g3.a$D[g3.a$T=="T2"] <- "D2"
g3.a$D[g3.a$T=="T3"] <- "D3"
#Need this data to force g3 into wide

g3.a.twide <- pivot_wider(g3.a, names_from=c(T), values_from=c(Tc), id_cols = serial_num)
g3.a.dwide <- pivot_wider(g3.a, names_from=c(D), values_from=c(Dc), id_cols = serial_num)
g3.a.wide <- merge(g3.a.twide, g3.a.dwide, by="serial_num")
g3.a.wide$cal.group <- g3.a$cal.group[match(g3.a.wide$serial_num, g3.a$serial_num)]
#Now it's in wide format.

g3.a.wide$T2 <- 23 - (2-g3.a.wide$T1)
g3.a.wide$T3 <- 37 - (2-g3.a.wide$T1)
#Here, I am estimating what T2 and T3 were basically assuming that the difference between T1 and the incubator temperature setting (the temperatures
#described for quadratic curve in Reeder et al. 2012) was constant across the 2, 23, and 37 degrees C trials. 

g3.a.wide$d.D1 <- g3.a.wide$T1 - g3.a.wide$D1
#Deviation between D1 and T1

g3.a.wide$D2 <- g3.a.wide$T2 - (predict.lm(d12.m, newdata = g3.a.wide))
g3.a.wide$D3 <- g3.a.wide$T3 - (predict.lm(d13.m, newdata = g3.a.wide))
#Predicted D2 and D3 values using the regressions of deviations. 

g3.a.wide$a <- (g3.a.wide$D1 - g3.a.wide$D3) + ((g3.a.wide$D2 - g3.a.wide$D1)*(g3.a.wide$T3 - g3.a.wide$T1))/(g3.a.wide$T2 - g3.a.wide$T1)
g3.a.wide$b <- g3.a.wide$T3^2 - g3.a.wide$T2*g3.a.wide$T3 + g3.a.wide$T2*g3.a.wide$T1 - g3.a.wide$T3*g3.a.wide$T1
g3.a.wide$c <- -(g3.a.wide$a/g3.a.wide$b)
g3.a.wide$B <- (g3.a.wide$D2 - g3.a.wide$D1 - g3.a.wide$c*((g3.a.wide$T2^2)-(g3.a.wide$T1^2)))/(g3.a.wide$T2-g3.a.wide$T1)
g3.a.wide$A <- g3.a.wide$D1 - g3.a.wide$B*g3.a.wide$T1 - g3.a.wide$c*(g3.a.wide$T1^2)
#Using the above values, we should be able to correct all temp data from group 3

g3.a.wide <- subset(g3.a.wide, select=-c(d.D1))
g12.a.wide <- subset(g12.a.wide, select=-c(d.D1, d.D2, d.D3))
cal.df <- rbind(g3.a.wide, g12.a.wide)
#This dataframe now contains all necessary parameters for calibrating raw data. 

master$A <- cal.df$A[match(master$serial_num, cal.df$serial_num)]
master$B <- cal.df$B[match(master$serial_num, cal.df$serial_num)]
master$c <- cal.df$c[match(master$serial_num, cal.df$serial_num)]
#Necessary parameters now looped into master df

master$value_cal <- (master$value - master$A - master$B*master$value - master$c*(master$value^2)) + master$value
#This column is the calibrated temperature data. 


#### Dataset trimming and writing ####

obs_per <- logger_meta %>%
  mutate(date_deployed = as.Date(date_deployed, format="%m/%d/%y")) %>%
  mutate(date_retrieved = as.Date(date_retrieved, format="%m/%d/%y")) %>%
  group_by(site) %>%
  summarise(date_deployed = mean(date_deployed), date_retrieved = mean(date_retrieved)) %>%
  mutate(sampling_length = date_retrieved - date_deployed)
#This dataframe contains information on when transmitters were deployed/retrieved

master$date_deployed <- obs_per$date_deployed[match(master$site, obs_per$site)]
master$date_retrieved <- obs_per$date_retrieved[match(master$site, obs_per$site)]
master <- master %>%
  filter(!(date<=date_deployed) & !(date>=date_retrieved))
#This dataframe now does not contain any data outside of the period that the transmitter was on a bat. 

master$datetime <- paste(master$date, master$time, sep=" ")
master$datetime <- as.POSIXct(strptime(master$datetime, "%Y-%m-%d %H:%M:%S",tz='EST')) 
#This column is necessary so that R can plot data over time on hourly scales. 

#p <- list()
#logg <- unique(master$uniq)
#for(i in 1:length(logg)){
  #p[[i]] <- list()
  #dat <- subset(master, uniq==logg[i])
  #p[[i]][[1]] <- ggplot(dat, aes(datetime,value_cal)) + 
   # geom_line(color="blue", size=0.5) + 
    #geom_line(aes(x=datetime, y=value), color="black", size=0.5, data=dat)+
    #scale_y_continuous(limits=c(-2,30))+
    #labs(x="Date", y=expression("Temperature " (degree*C)~" "), color="Legend")+
    #scale_color_manual(values=colors)+
    #ggtitle(paste(dat$uniq))+
    #theme(panel.background = element_blank(),
     #     legend.position = "top",
      #    axis.text = element_text(size=20, color='grey16'), 
       #   axis.title = element_text(size=25, color="grey16"),
        #  axis.ticks = element_line(size=0.8, color="grey57"),
         # panel.border = element_rect(colour = "grey57", fill=NA, size=0.8))
#}
#This stores a separate ggplot for each transmitter in a list called 'p'. Used for building .pdf of all data

mastercsv <- select(master, site, trans_id, date, time, value_cal, serial_num, logger_model, )
colnames(mastercsv) <- c("site","trans_id", "date","time","temp", "serial_num","logger_model")
#write.csv(mastercsv, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_raw_working.csv",row.names = F)
#DO NOT WRITE THIS CSV UNLESS ABSOLUTELY CERTAIN YOU WANT TO DO SO. THE SCRIPT "arousal_classification.R" adds two important columns, "behavior" and
#"event number" to this dataframe and then writes it as a .csv under the same name and location. Writing over this .csv in this script will force you
#to re-run the arousal classification script to retrieve those values, which takes several hours to run the for loops. 
