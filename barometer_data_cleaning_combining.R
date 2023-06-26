library(tidyverse)
library(stringr)
library(weathermetrics)

#First, creating a list of all the barometer dataframes:
dfs <- list.files(path = c("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Barometer data"), pattern = "*.csv", full.names = T) %>% #extracting file pathways
  lapply(function(x) {read.csv(x, header=F)}) %>% #reading in the .csv files as dataframes
  lapply(function(x) {x = x[-1:-2,-1]}) %>% #doing a little trimming away of metadata columns/rows.
  lapply(setNames, c("datetime", "pressure_hg", "temp", "rh")) #re-naming all the columns. 

#Currently, none of the dataframes have site metadata. However, that info is stored in their file name, so the below code uses string recognition from the package
#'stringr' to extract that site name from the list of file pathways: 
file.names <- list.files(path = c("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Barometer data"), pattern = "*.csv", full.names = T) #list of file pathways
site.pattern <- "(?<=/)[^/]+(?=_BAROMETER)" #This is 'regular expression pattern' that is used to extract the site name. 
file.sites <- vector(mode='list', length=6) #This is an empty list in which we will store the site names. 
for(i in 1:length(file.names)) {
  file.sites[i] = regmatches(file.names[i], regexpr(site.pattern, file.names[i], perl = TRUE))
  file.sites[i] = gsub("_"," ",file.sites[i])
} #This for-loop extracts the site names from the file pathways and stores them in their correct position in file.sites. It also replaces "_" in site names with " ". 


#In the below chunk of code, I am looping site metadata into each object in the dfs list, as well as adding a 'date' column and formatting the 
#'datetime' column. However, for some reason, the .csv file for South Lake is not reading the datetime column the same way the others are being 
#read and needs to be formatted differently. Below, I start with formatting the South Lake dataframe, and then move onto the others, which 
#are in positions 1:4 and 6 in the dfs list:
dfs[[5]]$site = file.sites[[5]] #Creates a new column, 'site,' with the site metadata extracted from its file pathway. 
dfs[[5]]$date = as.Date(str_extract(dfs[[5]]$datetime, "\\S+"), format="%m/%d/%y") #Separates the 'date' data from the 'datetime' column and makes it its own, formatted column
dfs[[5]]$datetime = as.POSIXct(strptime(dfs[[5]]$datetime, "%m/%d/%y %H:%M",tz='EST')) #Converts datetime into a POSIXct time object. Necessary for plotting. 
#^^^SOUTH LAKE dataframe
for(i in 1:4) {
  dfs[[i]]$site = file.sites[[i]] #Creates a new column, 'site,' in each dataframe with the site metadata extracted from its file pathway. 
  dfs[[i]]$date = as.Date(str_extract(dfs[[i]]$datetime, "\\d{2}/\\d{2}/\\d{2}"), format="%m/%d/%y") #Separates the 'date' data from the 'datetime' column and makes it its own, formatted column
  dfs[[i]]$datetime = as.POSIXct(strptime(dfs[[i]]$datetime, "%m/%d/%y %I:%M:%S %p",tz='EST')) #Converts datetime into a POSIXct time object. Necessary for plotting. 
}
#^^^Dataframes 1-4

dfs[[6]]$site = file.sites[[6]] #Creates a new column, 'site,' in each dataframe with the site metadata extracted from its file pathway. 
dfs[[6]]$date = as.Date(str_extract(dfs[[6]]$datetime, "\\d{2}/\\d{2}/\\d{2}"), format="%m/%d/%y") #Separates the 'date' data from the 'datetime' column and makes it its own, formatted column
dfs[[6]]$datetime = as.POSIXct(strptime(dfs[[6]]$datetime, "%m/%d/%y %I:%M:%S %p",tz='EST')) #Converts datetime into a POSIXct time object. Necessary for plotting. 

#^^^Dataframe 6. 

#The temperature data in Elroy Sparta is in Fahrenheit instead of Celsius, like the other dataframes. Need to correct that:
dfs[[3]]$temp <- fahrenheit.to.celsius(as.numeric(dfs[[3]]$temp), round=3)


#Now we can trim the barometer dataframes to not include any observations outside of the study period. 
start.end.dates <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/transmitter_metadata.csv") %>%
  mutate(date_deployed = as.Date(date_deployed, format="%m/%d/%y")) %>%
  mutate(date_retrieved = as.Date(date_retrieved, format="%m/%d/%y")) %>%
  group_by(site) %>%
  summarise(date_deployed = mean(date_deployed), date_retrieved = mean(date_retrieved)) %>%
  filter(site != "GRAPHITE")
start.dates<-as.list(start.end.dates$date_deployed)
end.dates<-as.list(start.end.dates$date_retrieved)
#This dataframe contains the start and end dates of the study period for each site. It is constructed from the transmitter metadata file.
#I then create two lists of start and end dates for each site, which are ordered in the same order as the dataframes in dfs

for(i in 1:length(dfs)) {
  dfs[[i]] = filter(dfs[[i]], date>start.dates[[i]] & date<end.dates[[i]])
}
#This for-loop trims each dataframe to not contain data outside the study period. Importantly, it trims away the day of barometer deployment/recovery. 

master <- do.call("rbind", dfs)
#Combined barometer database

master$section <- NA
#I need to add in the section metadata.
master$section[master$site=="BLACKBALL"]<-"R5"
master$section[master$site=="CP TUNNEL"]<-"1000"
master$section[master$site=="ELROY SPARTA"]<-"B"
master$section[master$site=="MEAD MINE"]<-"A"
master$section[master$site=="SOUTH LAKE"]<-"C"
master$section[master$site=="ZIMMERMAN"]<-3

#write.csv(master, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/barometer_master.csv", row.names = F)
