data <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/camera_data_IHM.csv")

#View(data)

unique(data$site)
unique(data$section)
unique(data$photo_time)
# we need to fix naming issues in site

data$site[data$site=="BLACKBALL "] <- "BLACKBALL"
data$site[data$site=="MEADMINE "] <- "MEADMINE"
data$site[data$site=="ZIMMERMAN "] <- "ZIMMERMAN"
# We fixed the naming issue

str(data)
#Need to change the date to be read as a date. 
#Also want to change site to factor. 

data$site <- as.factor(data$site)
data$section <- as.factor(data$section)
data$date <- as.Date(data$date, format="%m/%d/%y")
#Fixed our data type issue

head(data)
tail(data)
#looks good

dim(data)
#cool

names(data)
#looks good

nrow(data)
ncol(data)
#looks

cp.tunnel <- data[data$site=="CP TUNNEL ",]
#This dataframe is just data from CP tunnel

cp.summary <- aggregate(cluster_size~date, FUN=sum, data=cp.tunnel)
View(cp.summary)
#This is a summary of just cp tunnel

summ <- aggregate(cluster_size~site+section+date, FUN=sum, data=data)
View(summ)
#Summarized all of the data

plot(summ$cluster_size ~ summ$date)
#Base R plot

summ$uniq <- NA
summ$uniq <- paste(summ$site, summ$section)

library(ggplot2)



plot1 <- ggplot(aes(x=date, y=cluster_size, color=site, linetype=section), data=summ[summ$site=='MEADMINE',]) +
  geom_line()
plot1

geom_line()
geom_bar()
geom_histogram()



