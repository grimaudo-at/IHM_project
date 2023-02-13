library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lme4)
library(effects)

count<-read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/kate_cluster_master.csv")
mw<-read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/midwest_working.csv")

count$species<-toupper(count$species)
count$section<-toupper(count$section)
#capitalizing species column

count$date<-as.Date(count$date, format="%m/%d/%y")
#formatiting date

count<-separate(count, col=date, into=c('y','m','d'), sep="-", remove=F)
count$m<-as.numeric(count$m)
count$y<-as.numeric(count$y)
early<-filter(count, m>5)
late<-filter(count, m<5)
early$wyear<-early$y+1
late$wyear<-late$y
early$season<-"hiber_earl"
late$season<-"hiber_late"
count<-rbind(early, late)
#winter year

count<-filter(count, species=="MYLU")
count<-filter(count, section!="ALL")
bb<-filter(count, site=='BLACKBALL')
zm<-filter(count, site=="ZIMMERMAN")
sl<-filter(count, site=="SOUTH LAKE MINE")
el<-filter(count, site=="ELROY SPARTA")
nd<-filter(count, site=="NEDA MINE")
tm<-filter(count, site=="TAYLOR MINE")

#filtering out just MYLU blackball and zimmerman. Also dropping section names of "All", which are summareis of counts

bb.a<-aggregate(clustersize~site+wyear+season+section, data=bb, FUN=sum)
zm.a<-aggregate(clustersize~site+wyear+season+section, data=zm, FUN=sum)
sl.a<-aggregate(clustersize~site+wyear+season+section, data=sl, FUN=sum)
el.a<-aggregate(clustersize~site+wyear+season+section, data=el, FUN=sum)
nd.a<-aggregate(clustersize~site+wyear+season+section, data=nd, FUN=sum)
tm.a<-aggregate(clustersize~site+wyear+season+section, data=tm, FUN=sum)
#aggreating counts

bb.tot.a<-aggregate(clustersize~site+wyear+season, data=bb.a, FUN=sum)
zm.tot.a<-aggregate(clustersize~site+wyear+season, data=zm.a, FUN=sum)
sl.tot.a<-aggregate(clustersize~site+wyear+season, data=sl.a, FUN=sum)
el.tot.a<-aggregate(clustersize~site+wyear+season, data=el.a, FUN=sum)
nd.tot.a<-aggregate(clustersize~site+wyear+season, data=nd.a, FUN=sum)
tm.tot.a<-aggregate(clustersize~site+wyear+season, data=tm.a, FUN=sum)
#these dataframes contain total population sizes for each season in each site and winter year

bb.tot.a$unique<-paste(bb.tot.a$site, bb.tot.a$wyear, bb.tot.a$season)
zm.tot.a$unique<-paste(zm.tot.a$site, zm.tot.a$wyear, zm.tot.a$season)
sl.tot.a$unique<-paste(sl.tot.a$site, sl.tot.a$wyear, sl.tot.a$season)
el.tot.a$unique<-paste(el.tot.a$site, el.tot.a$wyear, el.tot.a$season)
nd.tot.a$unique<-paste(nd.tot.a$site, nd.tot.a$wyear, nd.tot.a$season)
tm.tot.a$unique<-paste(tm.tot.a$site, tm.tot.a$wyear, tm.tot.a$season)
bb.a$unique<-paste(bb.a$site, bb.a$wyear, bb.a$season)
zm.a$unique<-paste(zm.a$site, zm.a$wyear, zm.a$season)
sl.a$unique<-paste(sl.a$site, sl.a$wyear, sl.a$season)
el.a$unique<-paste(el.a$site, el.a$wyear, el.a$season)
nd.a$unique<-paste(nd.a$site, nd.a$wyear,nd.a$season)
tm.a$unique<-paste(tm.a$site, tm.a$wyear,tm.a$season)
#unique columns for matching

bb.a$total.mylu<-bb.tot.a$clustersize[match(bb.a$unique, bb.tot.a$unique)]
zm.a$total.mylu<-zm.tot.a$clustersize[match(zm.a$unique, zm.tot.a$unique)]
sl.a$total.mylu<-sl.tot.a$clustersize[match(sl.a$unique, sl.tot.a$unique)]
el.a$total.mylu<-el.tot.a$clustersize[match(el.a$unique, el.tot.a$unique)]
nd.a$total.mylu<-nd.tot.a$clustersize[match(nd.a$unique, nd.tot.a$unique)]
tm.a$total.mylu<-tm.tot.a$clustersize[match(tm.a$unique, tm.tot.a$unique)]
#bringing total count data into section aggregate dataframe

nd.a$prop.tot<-nd.a$clustersize/nd.a$total.mylu
#proportion value for Neda


### ZIMMERMAN

zm.a<- filter(zm.a, section=="2" | section=="3"|section=="3B"|section=="3R"|section=="4R"
              |section=="4"|section=="5"|section=="6"|section=="7"|section=="8"|section=="8A"|
                section=="9"|section=="11"|section=="9_6"|section=="9_10"|section=="12"|section=="10"
              |section=="13"|section=="14"|section=="15"|section=="22")


zm.earl<-filter(zm.a, season=="hiber_earl")
zm.late<-filter(zm.a, season=="hiber_late")
zm.a$unique<-paste(zm.a$wyear,zm.a$section)
zm.earl$unique<-paste(zm.earl$wyear,zm.earl$section)
zm.late$unique<-paste(zm.late$wyear,zm.late$section)
zm.a$earl_count<-zm.earl$clustersize[match(zm.a$unique, zm.earl$unique)]
zm.a$late_count<-zm.late$clustersize[match(zm.a$unique, zm.late$unique)]
zm.a$d.clust<-zm.a$late_count-zm.a$earl_count
zm.a$d.dir[zm.a$d.clust>0]<-"Increase"
zm.a$d.dir[zm.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize

zm.a<-filter(zm.a, wyear>=2019)

zm<-ggplot(aes(x=section, y=clustersize, fill=d.dir), data=subset(zm.a))+
  geom_bar(stat="identity")+
  ggtitle("Zimmerman")+
  ylab("Number of Bats")+
  xlab("Section")+
  facet_grid(wyear~season, scales = "free")+
  scale_fill_manual(values=c("red","darkgreen"))+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13)
  );zm





#### BLACKBALL

#lots of crazy section names; only want to use those that are on our maps:

bb.2020<-filter(bb.a, wyear==2020)
bb.sec<-c(unique(bb.2020$section))
#this vector contains all correct section names

bb.a<-subset(bb.a, section %in% bb.sec)

bb.earl<-filter(bb.a, season=="hiber_earl")
bb.late<-filter(bb.a, season=="hiber_late")
bb.a$unique<-paste(bb.a$wyear,bb.a$section)
bb.earl$unique<-paste(bb.earl$wyear,bb.earl$section)
bb.late$unique<-paste(bb.late$wyear,bb.late$section)
bb.a$earl_count<-bb.earl$clustersize[match(bb.a$unique, bb.earl$unique)]
bb.a$late_count<-bb.late$clustersize[match(bb.a$unique, bb.late$unique)]
bb.a$d.clust<-bb.a$late_count-bb.a$earl_count
bb.a$d.dir[bb.a$d.clust>0]<-"Increase"
bb.a$d.dir[bb.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize

bb.a<-filter(bb.a, wyear>=2019)

bb<-ggplot(aes(x=section, y=clustersize, fill=d.dir), data=subset(bb.a))+
  geom_bar(stat="identity")+
  ggtitle("Blackball")+
  ylab("Number of Bats")+
  xlab("Section")+
  facet_grid(wyear~season, scales = "free")+
  scale_fill_manual(values=c("red","darkgreen"))+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13)
  );bb




#### SOUTH LAKE
sl<-ggplot(aes(x=section, y=clustersize), data=subset(sl.a))+
  geom_bar(stat="identity")+
  ggtitle("South Lake")+
  ylab("Number of Bats")+
  xlab("Section")+
  facet_grid(wyear~season, scales = "free")+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13)
  );sl


#### ELROY SPARTA
el.earl<-filter(el.a, season=="hiber_earl")
el.late<-filter(el.a, season=="hiber_late")
el.a$unique<-paste(el.a$wyear,el.a$section)
el.earl$unique<-paste(el.earl$wyear,el.earl$section)
el.late$unique<-paste(el.late$wyear,el.late$section)
el.a$earl_count<-el.earl$clustersize[match(el.a$unique, el.earl$unique)]
el.a$late_count<-el.late$clustersize[match(el.a$unique, el.late$unique)]
el.a$d.clust<-el.a$late_count-el.a$earl_count
el.a$d.dir[el.a$d.clust>0]<-"Increase"
el.a$d.dir[el.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize


el<-ggplot(aes(x=section, y=clustersize, fill=d.dir), data=subset(el.a))+
  geom_bar(stat="identity")+
  ggtitle("ELROY SPARTA")+
  ylab("Number of Bats")+
  xlab("Section")+
  facet_grid(wyear~season, scales = "free")+
  scale_fill_manual(values=c("red","darkgreen"))+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13)
  );el


#### NEDA MINE

nd.earl<-filter(nd.a, season=="hiber_earl")
nd.late<-filter(nd.a, season=="hiber_late")
nd.a$unique<-paste(nd.a$wyear,nd.a$section)
nd.earl$unique<-paste(nd.earl$wyear,nd.earl$section)
nd.late$unique<-paste(nd.late$wyear,nd.late$section)
nd.a$earl_count<-nd.earl$clustersize[match(nd.a$unique, nd.earl$unique)]
nd.a$late_count<-nd.late$clustersize[match(nd.a$unique, nd.late$unique)]
nd.a$d.clust<-nd.a$late_count-nd.a$earl_count
nd.a$d.dir[nd.a$d.clust>0]<-"Increase"
nd.a$d.dir[nd.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize

nd.a.early<-filter(nd.a, season=="hiber_earl" & wyear=="2021")
nd.a.early<-select(nd.a.early, section, prop.tot)
nd.a.early <- nd.a.early[order(-nd.a.early$prop.tot),]
nd.lvl<-nd.a.early$section
nd.a$section <- factor(nd.a$section, levels=nd.lvl)
#re-ordering level of section factor for neda

nd.a$season[nd.a$season=="hiber_earl"]<-"Early Hibernation"
nd.a$season[nd.a$season=="hiber_late"]<-"Late Hibernation"
#Changing for labeling. 

nd<-ggplot(aes(x=section, y=prop.tot, fill=d.dir), data=subset(nd.a, wyear!="2022" & !is.na(section)))+
  geom_bar(stat="identity")+
  ylab("Proportion of total MYLU in site")+
  xlab("Sections")+
  facet_grid(wyear~season)+
  scale_fill_manual(values=c("red","darkgreen"))+
  labs(fill="Over-winter change in # of bats")+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13),
    axis.text.x = element_blank(),
  );nd


#### TAYLOR MINE

tm.earl<-filter(tm.a, season=="hiber_earl")
tm.late<-filter(tm.a, season=="hiber_late")
tm.a$unique<-paste(tm.a$wyear,tm.a$section)
tm.earl$unique<-paste(tm.earl$wyear,tm.earl$section)
tm.late$unique<-paste(tm.late$wyear,tm.late$section)
tm.a$earl_count<-tm.earl$clustersize[match(tm.a$unique, tm.earl$unique)]
tm.a$late_count<-tm.late$clustersize[match(tm.a$unique, tm.late$unique)]
tm.a$d.clust<-tm.a$late_count-tm.a$earl_count
tm.a$d.dir[tm.a$d.clust>0]<-"Increase"
tm.a$d.dir[tm.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize


tm<-ggplot(aes(x=section, y=clustersize, fill=d.dir), data=subset(tm.a))+
  geom_bar(stat="identity")+
  geom_text(aes(label = d.clust), vjust = -0.2)+
  ggtitle("TAYLOR MINE")+
  ylab("Number of Bats")+
  xlab("Section")+
  facet_grid(wyear~season, scales="free")+
  scale_fill_manual(values=c("red","darkgreen"))+
  theme(
    plot.title = element_text(size=20),
    axis.title = element_text(size=18),
    strip.text = element_text(size=15),
    axis.text = element_text(size=13)
  );tm



##### SECTION CHANGE IN TEMP STUFF

count<-filter(count, section!="" | section !="")
count.a<-aggregate(clustersize~site+wyear+season+section, data=count, FUN=sum)
#aggregating counts by site, wyear, season, and section for all sites

count.tot.a<-aggregate(clustersize~site+wyear+season, data=count.a, FUN=sum)
#this dataframe contains total population sizes for each season in each site and winter year

count.tot.a$unique<-paste(count.tot.a$site, count.tot.a$wyear, count.tot.a$season)
count.a$unique<-paste(count.a$site, count.a$wyear,count.a$season)
#unique columns for matching

count.a$total.mylu<-count.tot.a$clustersize[match(count.a$unique, count.tot.a$unique)]
#bringing total count data into section aggregate dataframe

count.earl<-filter(count.a, season=="hiber_earl")
count.late<-filter(count.a, season=="hiber_late")
count.a$unique<-paste(count.a$site, count.a$wyear,count.a$section)
count.earl$unique<-paste(count.earl$site, count.earl$wyear,count.earl$section)
count.late$unique<-paste(count.late$site, count.late$wyear,count.late$section)
count.a$earl_count<-count.earl$clustersize[match(count.a$unique, count.earl$unique)]
count.a$late_count<-count.late$clustersize[match(count.a$unique, count.late$unique)]
count.a$d.clust<-count.a$late_count-count.a$earl_count
count.a$d.dir[count.a$d.clust>0]<-"Increase"
count.a$d.dir[count.a$d.clust<0]<-"Decrease"
#characterizing sections based on if they increase or decrease in clustersize


count.a<-subset(count.a, select=-c(season, clustersize))
#removing unecessary columns

#there are several sections with counts of NA for a season. 
#To be conservative, going to assume these sections were not counted instead of being counts of 0. 
count.a<-filter(count.a, !is.na(earl_count) & !is.na(late_count))
#Removing NA counts. This leaves me with 1,184 site-sections-years with count data for both the start and end of hib

#Want to calculate a within-winter lambda, but don't want to divide by 0. Adding a constant of 1:
count.a$earl_count<-count.a$earl_count+1
count.a$late_count<-count.a$late_count+1

count.a$lambda<-count.a$late_count/count.a$earl_count
#Within-winter lambda

hist(count.a$lambda, breaks=40)
min(count.a$lambda)
max(count.a$lambda)
#serious variation in lambdas 

#Now going to calculate a change in temperature variable using MYLU roosting temperature. 
#This will be section-specific. Basically the change in avg roosting temperatures from early to late hib per section, not 
#necessarily on the same bats. 
 
#Need site, section, wyear, season-level averages

mw<-select(mw, site, section, date, temp, swab_type, species, gdL)
#Relevant columns

mw<-filter(mw, species=="MYLU")
#just mylu roosting temps

mw$date<-as.Date(mw$date, format="%m/%d/%y")
#formatiting date

mw<-separate(mw, col=date, into=c('y','m','d'), sep="-", remove=F)
mw$m<-as.numeric(mw$m)
mw$y<-as.numeric(mw$y)
mw.early<-filter(mw, m>=8)
mw.late<-filter(mw, m<5)
#only data between august and may, no summer
mw.early$wyear<-mw.early$y+1
mw.late$wyear<-mw.late$y
mw.early$season<-"hiber_earl"
mw.late$season<-"hiber_late"
mw<-rbind(mw.early, mw.late)
#winter year


mw.temp.a<-aggregate(temp~site+section+wyear+season, data=mw, FUN=mean)
mw.gdL.a<-aggregate(gdL~site+section+wyear+season, data=mw, FUN=mean)
#calculating avg roostig temps and loads for site-section-wyear-season
mw.temp.a.earl<-filter(mw.temp.a, season=="hiber_earl")
mw.temp.a.late<-filter(mw.temp.a, season=="hiber_late")
mw.temp.a.earl$unique<-paste(mw.temp.a.earl$site, mw.temp.a.earl$wyear, mw.temp.a.earl$section)
mw.temp.a.late$unique<-paste(mw.temp.a.late$site, mw.temp.a.late$wyear, mw.temp.a.late$section)
mw.gdL.a.earl<-filter(mw.gdL.a, season=="hiber_earl")
mw.gdL.a.late<-filter(mw.gdL.a, season=="hiber_late")
mw.gdL.a.earl$unique<-paste(mw.gdL.a.earl$site, mw.gdL.a.earl$wyear, mw.gdL.a.earl$section)
mw.gdL.a.late$unique<-paste(mw.gdL.a.late$site, mw.gdL.a.late$wyear, mw.gdL.a.late$section)
#unique ID's for merging

count.a$earl_temp<-mw.temp.a.earl$temp[match(count.a$unique, mw.temp.a.earl$unique)]
count.a$late_temp<-mw.temp.a.late$temp[match(count.a$unique, mw.temp.a.earl$unique)]
count.a$earl_load<-mw.gdL.a.earl$gdL[match(count.a$unique, mw.gdL.a.earl$unique)]
count.a$late_load<-mw.gdL.a.late$gdL[match(count.a$unique, mw.gdL.a.late$unique)]
#Bringing temp data into count df

count.a$d.temp<-count.a$late_temp-count.a$earl_temp
#change in roosting temp column

count.a$l_earl_load<-log10(count.a$earl_load)
count.a$l_late_load<-log10(count.a$late_load)
count.a$d.lgdL<-count.a$l_late_load-count.a$l_earl_load
#log10 change in loads column

count.a$log_lambda<-log10(count.a$lambda)
#log10 lambda

plot(count.a$log_lambda~count.a$d.temp)
#nothing LOL

#Want YOA data

yoa<-read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/1_yoa_midwest.csv")
count.a$yoa<-yoa$yoa[match(count.a$site, yoa$site)]
#yoa column

count.a$epi_year<-count.a$wyear-count.a$yoa
#epi year

g1<-ggplot(aes(x=d.temp, y=log_lambda, size=late_count), data=count.a)+
  geom_point()+
  facet_wrap(~epi_year);g1

bb.2<-filter(count.a, site=="BLACKBALL")

bb.p<-ggplot(aes(x=d.temp, y=log_lambda, size=late_count), data=bb.2)+
  geom_point()+
  facet_wrap(~epi_year);bb.p

#what about loads?
plot(count.a$d.lgdL~count.a$d.temp)


g2<-ggplot(aes(x=d.dir, y=d.lgdL, size=late_count), data=count.a)+
  geom_point();g2










#### INDIVIDUAL-LEVEL DATA
mw<-read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/DATA/midwest_working.csv")

mw<-select(mw, site, section, date, temp, swab_type, species, gdL, band)
#Relevant columns

mw<-filter(mw, species=="MYLU")
#just mylu roosting temps

mw$date<-as.Date(mw$date, format="%m/%d/%y")
#formatiting date

mw<-separate(mw, col=date, into=c('y','m','d'), sep="-", remove=F)
mw$m<-as.numeric(mw$m)
mw$y<-as.numeric(mw$y)
mw.early<-filter(mw, m>=8)
mw.late<-filter(mw, m<5)
#only data between august and may, no summer
mw.early$wyear<-mw.early$y+1
mw.late$wyear<-mw.late$y
mw.early$season<-"hiber_earl"
mw.late$season<-"hiber_late"
#winter year

mw.early<-filter(mw.early, !is.na(band))
mw.late<-filter(mw.late, !is.na(band))
mw.early$unique<-paste(mw.early$site, mw.early$wyear, mw.early$band)
mw.late$unique<-paste(mw.late$site, mw.late$wyear, mw.late$band)
#unique ID's for matching

mw.ind<-mw.early
mw.ind$gdL.late<-mw.late$gdL[match(mw.ind$unique, mw.late$unique)]
mw.ind$temp.late<-mw.late$temp[match(mw.ind$unique, mw.late$unique)]
#matching in loads and temp

mw.ind$d.temp<-mw.ind$temp.late-mw.ind$temp
#Individual's change in roosting temp
mw.ind$d.gdL<-mw.ind$gdL.late-mw.ind$gdL
#Individual's change in loads

mw.ind<-filter(mw.ind, !is.na(d.temp) & !is.na(d.gdL))
#we have 217 observations on individual bats for both a change in roosting temp and loads

mw.ind$lgdL<-log10(mw.ind$gdL)
mw.ind$lgdL.late<-log10(mw.ind$gdL.late)
mw.ind$d.lgdL<-mw.ind$lgdL.late-mw.ind$lgdL
#change in log10 loads

plot(mw.ind$d.lgdL~mw.ind$d.temp)
plot(mw.ind$lgdL.late~mw.ind$d.temp)
plot(mw.ind$d.lgdL~mw.ind$temp.late)
plot(mw.ind$d.gdL~mw.ind$d.temp)


m1<-lmer(d.lgdL~d.temp*temp+(1|site), data=mw.ind)
summary(m1)
plot(allEffects(m1))
#Nothing?

Anova(m1)

mw.ind$yoa<-yoa$yoa[match(mw.ind$site, yoa$site)]
mw.ind$epi_year<-mw.ind$wyear-mw.ind$yoa
#epidemic year


#any effect of epidemic year?
g3<-ggplot(aes(x=d.temp, y=d.lgdL), data=mw.ind)+
  geom_point()+
  facet_wrap(~epi_year);g3

m2<-lmer(d.lgdL~d.temp*epi_year+(1|site), data=mw.ind)
summary(m2)
plot(allEffects(m2))
#Nothing

## What if we give a gdL value equivalent to 40ct for those bats that had band info read but had NA gdL values?
mw.early.ct<-mw.early
mw.late.ct<-mw.late
mw.early.ct$gdL[is.na(mw.early.ct$gdL)]<-10^((40-22.04942)/-3.34789)
mw.late.ct$gdL[is.na(mw.late.ct$gdL)]<-10^((40-22.04942)/-3.34789)

mw.ind.ct<-mw.early.ct
mw.ind.ct$gdL.late<-mw.late.ct$gdL[match(mw.ind.ct$unique, mw.late.ct$unique)]
mw.ind.ct$temp.late<-mw.late.ct$temp[match(mw.ind.ct$unique, mw.late.ct$unique)]
#matching in loads and temp

mw.ind.ct$d.temp<-mw.ind.ct$temp.late-mw.ind.ct$temp
#Individual's change in roosting temp
mw.ind.ct$d.gdL<-mw.ind.ct$gdL.late-mw.ind.ct$gdL
#Individual's change in loads

mw.ind.ct<-filter(mw.ind.ct, !is.na(d.temp) & !is.na(d.gdL))
#we have 217 observations on individual bats for both a change in roosting temp and loads

mw.ind.ct$lgdL<-log10(mw.ind.ct$gdL)
mw.ind.ct$lgdL.late<-log10(mw.ind.ct$gdL.late)
mw.ind.ct$d.lgdL<-mw.ind.ct$lgdL.late-mw.ind.ct$lgdL
#change in log10 loads

plot(mw.ind.ct$d.lgdL~mw.ind.ct$d.temp)
plot(mw.ind.ct$lgdL.late~mw.ind.ct$d.temp)
plot(mw.ind.ct$d.lgdL~mw.ind.ct$temp.late)
plot(mw.ind.ct$d.gdL~mw.ind.ct$d.temp)

m3<-lmer(lgdL.late~d.temp+(1|site), data=mw.ind.ct)
summary(m3)
plot(allEffects(m3))
#nothing?


