dat<-read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv")

library(tidyverse)
g.21<-filter(dat, model=="DS1921G")
l.22<-filter(dat, model=="DS1922L")

section<-c("1","2","3")

set.seed(201)

g.21$section[sample(1:nrow(g.21), nrow(g.21), replace=F)] = rep(section, c(12,12,11))
l.22$section[sample(1:nrow(l.22), nrow(l.22), replace=F)] = rep(section, c(12,12,11))
#randomly assigning sections to loggers

dat.gt<-rbind(g.21, l.22)

dat.gt.1<-filter(dat.gt, section=="1")
dat.gt.2<-filter(dat.gt, section=="2")
dat.gt.3<-filter(dat.gt, section=="3")

dat.gt.1$id<-sample(1:24, replace=F)
dat.gt.2$id<-sample(25:48, replace=F)
dat.gt.3$id<-sample(49:70, replace=F)

dat.gt<-rbind(dat.gt.1, dat.gt.2, dat.gt.3)
#write.csv(dat.gt, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv", row.names = F)

dat.el<-filter(dat, is.na(weight_g))

el.21<-filter(dat.el, model=="DS1921G")
el.22<-filter(dat.el, model=="DS1922L")

section<-c("1","2")
set.seed(201)

el.21$section[sample(1:nrow(el.21), nrow(el.21), replace=F)] = rep(section, c(2,3))
el.22$section[sample(1:nrow(el.22), nrow(el.22), replace=F)] = rep(section, c(8,7))

el<-rbind(el.21,el.22)

el.1<-filter(el, section=="1")
el.2<-filter(el, section=="2")

el.1$id<-sample(71:80, replace=F)
el.2$id<-sample(81:90, replace=F)

dat.el<-rbind(el.1,el.2)
dat.gt<-filter(dat, !is.na(weight_g))

dat<-rbind(dat.gt, dat.el)
dat$section<-NA
dat$site[!is.na(dat$weight_g)]<-"GRAPHITE"
dat$site[is.na(dat$weight_g)]<-"ELROY SPARTA"
#write.csv(dat, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv", row.names = F)


dat<-read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv")
dat.gt<-filter(dat, site=="GRAPHITE"|site=="ELROY SPARTA")
dat.ngt<-filter(dat, site!='GRAPHITE'&site!="ELROY SPARTA")
library(tidyverse)
g.21<-filter(dat.ngt, model=="DS1921G")
l.22<-filter(dat.ngt, model=="DS1922L")

site<-c("CP TUNNEL", "SOUTH LAKE", "MEAD MINE", "BLACKBALL", "ZIMMERMAN")
g.21$site[sample(1:nrow(g.21), nrow(g.21), replace=F)] = rep(site, c(5,25,20,20,24))
l.22$site[sample(1:nrow(l.22), nrow(l.22), replace=F)] = rep(site, c(15,25,20,20,25))

asn<-rbind(g.21, l.22)

asn$id[asn$site=="CP TUNNEL"]<-sample(91:110, replace=F)
asn$id[asn$site=="SOUTH LAKE"]<-sample(111:160, replace=F)
asn$id[asn$site=="MEAD MINE"]<-sample(161:200, replace=F)
asn$id[asn$site=="BLACKBALL"]<-sample(201:240, replace=F)
asn$id[asn$site=="ZIMMERMAN"]<-sample(241:289, replace=F)

dat.new<-rbind(asn, dat.gt)
#write.csv(dat.new, "/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv", row.names = F)
