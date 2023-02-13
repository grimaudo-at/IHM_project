library(tidyverse)
library(berryFunctions)

ex <- as.data.frame(read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Transmitter data/GRAPHITE_41.csv", header=F))
ex2 <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Transmitter data/GRAPHITE_30.csv", header=F)
ex3 <- read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/Transmitter data/MEAD_MINE_192.csv", header=F)

dfs <- list(ex, ex2, ex3)

l <- lapply(dfs, function(x) {if(str_detect(x[1,1], 'DS2422')=="TRUE") {"DS1922L"}
  else{"DS1921G"}})

l2 <- lapply(dfs, function(x) {str_sub(x[2,1],-8,-1)})

dfs2 <- lapply(dfs, function(x) {if(str_detect(x[1,1], 'DS2422')=="TRUE") {x<-x[-c(1:21),]}
  else{x<-x[-c(1:16),]}})

dfs2<-lapply(dfs2, function(x) {as.data.frame(x)})

for(i in seq_along(dfs2)) {
  dfs2[[i]]$thing <- rep(c("datetime","unit","value"))
} 

pivot <- function(t) {pivot_wider(t, names_from=thing, values_from = x)}
wide.dat <- lapply(dfs2, pivot)

poo<-as.data.frame(wide.dat[2]);View(poo)
