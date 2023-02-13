library(tidyverse)
library(berryFunctions)
library(chron)

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

master$time2<-NA
for(i in 1:nrow(master)) {if(str_detect(master[i,1], 'AM')=="TRUE") {master[i,8] <- "AM"}
  else{master[i,8] <- "PM"}}
master<-separate(master, datetime, into=c("date","time"),sep=" ")
master$time <- paste(master$time, master$time2, sep=" ")
master<-subset(master, select=-c(time2))
master$time <- format(strptime(master$time, "%I:%M:%S %p"), "%H:%M:%S")
master$date <- as.Date(master$date, format="%m/%d/%y")
#Above chunk of code is correcting format of dates and times. 

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
  summarise(D = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
g2.a <- g2 %>%
  filter(!is.na(T)) %>%
  group_by(serial_num, cal.group, T) %>%
  summarise(D = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
g3.a <- g3 %>%
  filter(!is.na(T)) %>%
  group_by(serial_num, cal.group, T) %>%
  summarise(D = mean(value)) %>%
  mutate(uniq = paste(cal.group, T, sep="_"))
#Summarizing each transmitter's D values for each calibration period. Uniq is a column for matching
t.coup$uniq <- paste(t.coup$group, t.coup$T, sep="_")
#Identical column for matching in the thermocouple dataframe

g1.a$Tc <- t.coup$value[match(g1.a$uniq, t.coup$uniq)]
g2.a$Tc <- t.coup$value[match(g2.a$uniq, t.coup$uniq)]
g3.a$Tc <- t.coup$value[match(g3.a$uniq, t.coup$uniq)]
#Bringing in T data from thermocouples. 




