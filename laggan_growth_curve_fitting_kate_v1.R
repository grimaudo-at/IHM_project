## finding R's for isolates 

#10/10/21
#UPDATED: 1/10/22

# next steps: go through and adjust target R for temps:3,7,11,15,21
  #annotate with notes for warnings and all NA's 

# all NA NO: TY2012,




#clear environment  
rm(list=ls())

#loading in packages
library(ggplot2)
library(MASS)
library(reshape2)
library(Matrix)
#devtools::install_github("RobinHankin/Brobdingnag")
library(tidyverse)
library(lme4)
library(effects)
library(emmeans)
library(dplyr)
library(tidyr)
library(viridis)
library(viridisLite)
library(cowplot)

#reading in data 
#cult = read.csv("DATA/nl_merge_culture_temp.csv") #kate

cult = read.csv("/Users/alexg8/Dropbox/MIDWEST_WNS/nichole_model_fitting/nl_merge_culture_temp_1_6_22.csv")
View(cult)



head(cult)
unique(cult$plate_id)
unique(cult$sampling_pt)
unique(cult$temp)
unique(cult$isolate_id)


# setting lgdL0 
cult$lgdL0=cult$lgdL
cult$lgdL0[cult$lgdL0=="NA"]=NA
cult$lgdL0[is.na(cult$lgdL0)==T] = -5.75
unique(cult$lgdL0)


# subset to only dilutions of 100
cult$dilution=as.factor(cult$dilution)
cult=subset(cult, dilution=="100"); cult$dilution = droplevels(cult$dilution); unique(cult$dilution)

#removing old assay dates
cult=subset(cult, assay_date.1!="2021-06-17"); cult$assay_date.1 = droplevels(cult$assay_date.1); unique(cult$assay_date.1)
cult=subset(cult, assay_date.1!="2021-06-18"); cult$assay_date.1 = droplevels(cult$assay_date.1); unique(cult$assay_date.1)
cult = subset(cult, assay_id.1!="P1A_1to100_dilution_Laggan_Blehert_07Aug2021_KLP")

##look at the data
ggplot(data=cult, aes(x=day, y=lgdL0, color=isolate_id))+
  geom_point()+
  stat_summary()+
  facet_wrap(~temp)

#setting starting plate across all temps:
cult.start = cult %>%
  filter(temp=="X")
new.cult = cult %>%
  filter(temp!="X")

start.temp1 = cult.start
start.temp1$temp=1
start.temp3 = cult.start
start.temp3$temp=3
start.temp5 = cult.start
start.temp5$temp=5
start.temp7 = cult.start
start.temp7$temp=7
start.temp9 = cult.start
start.temp9$temp=9
start.temp11 = cult.start
start.temp11$temp=11
start.temp13 = cult.start
start.temp13$temp=13
start.temp15 = cult.start
start.temp15$temp=15
start.temp18 = cult.start
start.temp18$temp=18
start.temp21 = cult.start
start.temp21$temp=21

all.start = rbind(start.temp1, start.temp3, start.temp5,start.temp7, start.temp9, start.temp11, 
                  start.temp13, start.temp15, start.temp18, start.temp21)

dim(all.start)
dim(cult.start)

new.cult = rbind(all.start, new.cult)
dim(new.cult)

#manually adjusting dilutions for start plate
new.cult$gdL2 = new.cult$gdL
new.cult$gdL2[new.cult$day==0]= new.cult$gdL2[new.cult$day==0]/100

#new lgdL with proper dilutions 
new.cult$lgdL = log10(new.cult$gdL2)
new.cult$log.num.conidia = 5.34 + (new.cult$lgdL)*1.05
new.cult$num.conidia = 10^(new.cult$log.num.conidia)

#making NA's in gdL to be 0 (gdL0)
new.cult$gdL0 = new.cult$gdL
new.cult$gdL0[is.na(new.cult$gdL0)==T]=0

#removing blanks from dataset
new.cult = subset(new.cult,isolate_id!="BLNK")

####Looking at a few summary models ####
#setting day and temp as factors
new.cult$day.f = as.factor(new.cult$day)
new.cult$temp.f = as.factor(new.cult$temp)

new.cult$day.f = as.factor(new.cult$day)
new.cult$temp.f = as.factor(new.cult$temp)

library(lme4)
mod1 = lmer(lgdL~temp.f*isolate_id +(1|day.f), data=new.cult)
summary(mod1)

library(effects)
plot(allEffects(mod1))
emmeans(mod1, ~ temp.f|isolate_id)


mod2 = lmer(lgdL~isolate_id +(1|day.f) +(1|temp.f), data=new.cult)
summary(mod2)

library(effects)
plot(allEffects(mod2))
emip = emmeans(mod2,"isolate_id" ) 
pairs(emip)

### labelling pre and post isolates ###
new.cult$isolate_id = as.character(new.cult$isolate_id)
new.cult$isolate_type = new.cult$isolate_id
new.cult$collection.year = parse_number(new.cult$isolate_id)
new.cult$collection.year[new.cult$collection.year==1]=NA
new.cult$collection.year[new.cult$collection.year==2]=NA
new.cult$collection.year[new.cult$isolate_id=="ATCC"]=2008
new.cult$site = substr(new.cult$isolate_id, 1, 2)
new.cult$site[new.cult$site=="AT"] = "WH"

new.cult$stage [new.cult$isolate_id=="WH2008"]="pre"
new.cult$stage [new.cult$isolate_id=="ATCC"]="pre"
new.cult$stage [new.cult$isolate_id=="HL2007"]=NA
new.cult$stage [new.cult$isolate_id=="SC2007"]="pre"
new.cult$stage [new.cult$isolate_id=="BM2009"]="pre"
new.cult$stage [new.cult$isolate_id=="FS2009"]="pre"

new.cult$stage [new.cult$isolate_id=="WH2020"]="post"
new.cult$stage [new.cult$isolate_id=="SC2020"]="post"
new.cult$stage [new.cult$isolate_id=="BM2019"]="post"
new.cult$stage [new.cult$isolate_id=="FS2020"]="post"

### labelling internatinal isolates ###
new.cult$isolate_type="international"
new.cult$isolate_type[new.cult$collection.year>1000]="us"

#### summary models of collection year and stage ####
mod3 = lmer(lgdL~collection.year +(1|day.f) +(1|temp.f), data=new.cult)
summary(mod3)
plot(allEffects(mod3))

mod4 = lmer(lgdL~site*collection.year +(1|day.f) +(1|temp.f), data=new.cult)
summary(mod4)


mod5 = lmer(lgdL~site*stage +(1|day.f) +(1|temp.f), data=new.cult)
summary(mod5)

library(effects)
plot(allEffects(mod5))
emip = emmeans(mod5,~stage|site ) 
pairs(emip)

mod6 = lmer(lgdL~isolate_type +(1|day.f) +(1|temp.f), data=new.cult)
summary(mod6)

plot(allEffects(mod5))
emip = emmeans(mod5,~stage|site ) 
pairs(emip)


##### FITTING LOGISTIC MODELS ####

new.cult$gdL
summary(new.cult$gdL)

tidy.cult = new.cult %>%
  group_by(temp,isolate_id)%>%
  summarise(gdL = mean(gdL, na.rm=T), max.gdL = max(gdL, na.rm=T), min.gdL = min(gdL, na.rm=T))
tidy.cult

tidy.cult2 = new.cult %>%
  group_by(isolate_id)%>%
  summarise(max.gdL = max(gdL, na.rm=T))
tidy.cult2
View(tidy.cult2)

#Set K first. We tried 0.5 to 1.5 in 0.05 inteverals and 0.85 minimized RMSE (see below)
K=283

#######KATE's FOR LOOP - NLS fits ##############
new.cult2 = subset(new.cult, isolate_id!="BLNK")
new.cult2 = subset(new.cult, isolate_id!="BLNK ")
unique.temps = unique(new.cult2$temp)
unique.isolates = unique(new.cult2$isolate_id)
total.groups = expand.grid(unique.temps = unique.temps, unique.isolates=unique.isolates)
total.groups$comb.group = paste(total.groups$unique.isolates, total.groups$unique.temps, sep ="_")
new.cult2$group = paste(new.cult2$isolate_id, new.cult2$temp, sep = "_")
N0=NA
r=NA
r.se = NA
dat.tog = data.frame(r, r.se, comb.group = NA)

for (i in 1:nrow(total.groups) ){
  Data = new.cult2 %>%
    filter(day>0 & new.cult2$group==total.groups$comb.group[i]) %>%
    arrange(day)
  #[new.cult2$Day>0 & new.cult2$group = total.groups$comb.group[i]]
  start = cult.start[cult.start$isolate_id==total.groups$unique.isolates[i],]
  N0 = start$gdL
  N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
  try (
  Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
                data=Data,
                start=list(r=0.6), trace=TRUE),
  #error = function(e){
  #  message(paste("An error occurred for item", i, stuff[[i]],":\n"), e)
  #}
  )
  summary(Testset1)
  # Recording summary coefs #
  r = summary(Testset1)$coef[1] 
  r.se = summary(Testset1)$coef[2]
  dat = data.frame(r, r.se, comb.group = total.groups$comb.group[i])
  dat.tog = bind_rows(dat.tog, dat)
}

#dat.tog2 = dat.tog %>%
#  tidyr::separate_wider_delim(comb.group, "_", names = c("isolate_type", "temp"))

x = colsplit(dat.tog$comb.group, "_", c("isolate_id", "temp"))
dat.tog$isolate_id = x$isolate_id
dat.tog$temp = x$temp

dat.tog.k283 = dat.tog


write.csv(dat.tog.k283, "nichole_model_fitting/dat.tog.k283.csv")

##can start here##
dat.tog = read.csv("Dropbox/MIDWEST_WNS/nichole_model_fitting/dat.tog.k283.csv")
head(dat.tog)
unique(dat.tog$isolate_id)

##make some graphs of the NLS fits
dat.tog = dat.tog %>%
  mutate(time.period = str_sub(isolate_id, 3, 6))%>%
  mutate(isolate.place = str_sub(isolate_id, 1, 2))

dat.tog$time.period.n = as.numeric(dat.tog$time.period) 

dat.tog2 = dat.tog %>%
  filter(time.period.n>3)%>%
  drop_na(time.period.n) %>%
  mutate(phase = ifelse(time.period.n>2013, "contemporary", "historic"))

dat.tog2$time.period.nf=as.factor(dat.tog2$time.period.n)

library(ggthemes)
g1 = ggplot(data=dat.tog2, aes(x=temp, y=r, color=phase))+
  geom_point(size=2)+
  geom_smooth(method = "gam",se=F, size=1.5)+
  #geom_smooth(method="loess", span=1, se=F)+
  ylim(0, 0.75)+
  #scale_colour_manual(name="Treatment", values=c("#06B1AB","#055387"))+
  scale_colour_manual(name="Treatment", values=c("#076DB0","#E92024"))+
  #xlim(0,18)+
  #coord_cartesian(ylim = c(0, 0.75))+
  xlab(expression("Temperature "( degree*C)))+
  ylab ("Growth rate")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text.x=element_text(hjust=1,size=12),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        axis.line=element_line(),
        legend.text = element_text(size=12), 
        legend.position = "top")
  
g1


  

#######KATE's FOR LOOP - BRMS fits ##############
library(brms)
new.cult2 = subset(new.cult, isolate_id!="BLNK")
new.cult2 = subset(new.cult, isolate_id!="BLNK ")
unique.temps = unique(new.cult2$temp)
unique.isolates = unique(new.cult2$isolate_id)
total.groups = expand.grid(unique.temps = unique.temps, unique.isolates=unique.isolates)
total.groups$comb.group = paste(total.groups$unique.isolates, total.groups$unique.temps, sep ="_")
new.cult2$group = paste(new.cult2$isolate_id, new.cult2$temp, sep = "_")
N0=NA
r=NA
r.se = NA
k=NA
k.se=NA
dat.tog = data.frame(r, r.se, comb.group = NA, k, k.se)
dim(new.cult2)
new.cult3 = new.cult2 %>%
  group_by(swab_id, temp, day,group,isolate_type, collection.year, site, stage)%>%
  summarise(gdL = mean(gdL, na.rm=T))
dim(new.cult3)

for (i in 1:nrow(total.groups) ){
#for (i in 112:nrow(total.groups) ){
  Data = new.cult2 %>%
    filter(new.cult2$day>0 & new.cult2$group==total.groups$comb.group[i]) %>%
    arrange(day)
  Data = Data %>%
    group_by(swab_id,temp, day,group,isolate_type, collection.year, site, stage)%>%
    summarise(gdL = mean(gdL, na.rm=T))
  #[new.cult2$Day>0 & new.cult2$group = total.groups$comb.group[i]]
  start = cult.start[cult.start$isolate_id==total.groups$unique.isolates[i],]
  N0 = start$gdL
  N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
  Data$N0 =sample(N0, nrow(Data), replace=T)
 # try (
  hillprior <- c(
    set_prior("normal(0, 2)", nlpar = "r", lb=0, ub=2),
   # set_prior("normal(283, 2)", nlpar = "K"),
    set_prior("uniform(282, 284)", nlpar = "K", lb=282, ub=284),
    set_prior("normal(0, 10)", class="sigma"))
  
  log.fit.formula = bf(
        gdL ~ K/
          (1+((K/N0)-1)*exp(-r*day)), 
        #nonlinear variables
        r + K ~ 1, 
        #nonlinear fit
        nl = TRUE)
    
Testset.bayes <- brm(log.fit.formula, 
                     family = gaussian(),
                     data = Data, 
                     prior = hillprior,
                     control = list(adapt_delta = 0.95, max_treedepth = 13), 
                     cores=2)
summary(Testset.bayes)
  # Recording summary coefs #
  r = fixef(Testset.bayes)[1,1]
  r.se = fixef(Testset.bayes)[1,2]
  k = fixef(Testset.bayes)[2,1]
  k.se = fixef(Testset.bayes)[2,2]
  dat = data.frame(r, r.se, comb.group = total.groups$comb.group[i], k, k.se)
  dat.tog = bind_rows(dat.tog, dat)
}

unique(dat.tog$comb.group)
#write.csv(dat.tog, "nichole_model_fitting/dat.tog.bayes_08MAR2023.csv")

x = colsplit(dat.tog$comb.group, "_", c("isolate_id", "temp"))
dat.tog$isolate_id = x$isolate_id
dat.tog$temp = x$temp

library(stringr)
library(tidyverse)
dat.tog = dat.tog %>%
  mutate(time.period = str_sub(comb.group, 3, 6))%>%
  mutate(isolate.place = str_sub(comb.group, 1, 2))


dat.tog$time.period.n = as.numeric(dat.tog$time.period) 

dat.tog2 = dat.tog %>%
  filter(time.period.n>3)%>%
  drop_na(time.period.n) %>%
  mutate(phase = ifelse(time.period.n>2013, "contemporary", "historic"))

dat.tog2$time.period.nf=as.factor(dat.tog2$time.period.n)

library(ggthemes)
g1 = ggplot(data=dat.tog2, aes(x=temp, y=r, color=phase))+
  geom_point(size=2)+
  geom_smooth(method = "gam",se=F, size=1.5)+
  #geom_smooth(method="loess", span=1, se=F)+
  ylim(0, 0.75)+
  #scale_colour_manual(name="Treatment", values=c("#06B1AB","#055387"))+
  scale_colour_manual(name="Treatment", values=c("#076DB0","#E92024"))+
  #xlim(0,18)+
  #coord_cartesian(ylim = c(0, 0.75))+
  xlab(expression("Temperature "( degree*C)))+
  ylab ("Growth rate")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text.x=element_text(hjust=1,size=12),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        axis.line=element_line(),
        legend.text = element_text(size=12), 
        legend.position = "top")

g1


##plot another way
g1 = ggplot(data=dat.tog2, aes(x=as.factor(temp), y=r, fill=phase))+
  geom_boxplot(position = position_dodge2(.5))+
  geom_point(size=2, position = position_dodge2(1), aes(color=isolate.place))+
  
#  ylim(0, 0.75)+
  #scale_colour_manual(name="Treatment", values=c("#06B1AB","#055387"))+
  #scale_colour_manual(name="Treatment", values=c("#076DB0","#E92024"))+
  #xlim(0,18)+
  #coord_cartesian(ylim = c(0, 0.75))+
  xlab(expression("Temperature "( degree*C)))+
  ylab ("Growth rate")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text.x=element_text(hjust=1,size=12),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        axis.line=element_line(),
        legend.text = element_text(size=12), 
        legend.position = "top")

g1


##plot dat tog and group by isolate
g3 = ggplot(data=dat.tog, aes(x=temp, y=r, color=isolate_id))+
  geom_point(size=2)+
 # geom_smooth(method = "gam",se=F, size=1.5)+
  geom_smooth(method="loess", span=1, se=F)+
  ylim(0, 0.75)+
  #scale_colour_manual(name="Treatment", values=c("#06B1AB","#055387"))+
  #scale_colour_manual(name="Treatment", values=c("#076DB0","#E92024"))+
  #xlim(0,18)+
  #coord_cartesian(ylim = c(0, 0.75))+
  xlab(expression("Temperature "( degree*C)))+
  ylab ("Growth rate")+
  theme_bw()+
  theme(axis.title=element_text(size=15),
        axis.text.x=element_text(hjust=1,size=12),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        axis.line=element_line(),
        legend.text = element_text(size=12), 
        legend.position = "top")

g3


#quick analysis
lm1 = lm(r~as.factor(temp)*phase, data =dat.tog2 )
summary(lm1)
anova(lm1)
library(emmeans)
emip = emmeans(lm1, ~as.factor(temp)*phase)
emip
pairs(emip)

##next steps - curve fitting
View(dat.tog.k283)


#######NICHOLE'S MILLION LINES OF CODE!##################
#### HL2007 temps 1, 5, 9, 13, 18 ####
# HL2007 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t1 = summary(Testset1)$coef[1] 
HL2007.se.t1 = summary(Testset1)$coef[2] 

# HL2007 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t3 = summary(Testset1)$coef[1] 
HL2007.se.t3 = summary(Testset1)$coef[2] 


# HL2007 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t5 = summary(Testset1)$coef[1] 
HL2007.se.t5 = summary(Testset1)$coef[2] 

# HL2007 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.51), trace=TRUE) #QARNINGS
summary(Testset1)
# Recording summary coefs #
HL2007.r.t7 = summary(Testset1)$coef[1] 
HL2007.se.t7 = summary(Testset1)$coef[2] 

# HL2007 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)   # WARNININGS
summary(Testset1)
# Recording summary coefs #
HL2007.r.t9 = summary(Testset1)$coef[1] 
HL2007.se.t9 = summary(Testset1)$coef[2] 

# HL2007 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)   
summary(Testset1)
# Recording summary coefs #
HL2007.r.t11 = summary(Testset1)$coef[1] 
HL2007.se.t11 = summary(Testset1)$coef[2] 

# HL2007 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t13 = summary(Testset1)$coef[1] 
HL2007.se.t13 = summary(Testset1)$coef[2] 

# HL2007 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t15 = summary(Testset1)$coef[1] 
HL2007.se.t15 = summary(Testset1)$coef[2] 


# HL2007 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# Recording summary coefs #
HL2007.r.t18 = summary(Testset1)$coef[1] 
HL2007.se.t18 = summary(Testset1)$coef[2]

# HL2007 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="HL2007",]
start = cult.start[cult.start$isolate_id=="HL2007",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.1), trace=TRUE) #WARNINGS
summary(Testset1)
# Recording summary coefs #
HL2007.r.t21 = summary(Testset1)$coef[1] 
HL2007.se.t21 = summary(Testset1)$coef[2] 


# logging HL2007
#setting temp in Haile's dataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(HL2007.r.t1, HL2007.r.t3,HL2007.r.t5, HL2007.r.t7, HL2007.r.t9,HL2007.r.t11, HL2007.r.t13,HL2007.r.t15, HL2007.r.t18, HL2007.r.t21)
se = c(HL2007.se.t1, HL2007.se.t3,HL2007.se.t5, HL2007.se.t7, HL2007.se.t9,HL2007.se.t11, HL2007.se.t13,HL2007.se.t15, HL2007.se.t18, HL2007.se.t21 )
HL2007.r.se.tall = data.frame(temp, r, se)
HL2007.r.se.tall$isolate_id = "HL2007"



#### SC2007 temps 1, 5, 9, 13, 18 ####
#SC2007 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t1 = summary(Testset1)$coef[1] 
SC2007.se.t1  = summary(Testset1)$coef[2] 

#SC2007 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t3 = summary(Testset1)$coef[1] 
SC2007.se.t3  = summary(Testset1)$coef[2] 

#SC2007 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t5 = summary(Testset1)$coef[1] 
SC2007.se.t5  = summary(Testset1)$coef[2] 

#SC2007 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t7 = summary(Testset1)$coef[1] 
SC2007.se.t7  = summary(Testset1)$coef[2] 

#SC2007 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t9 = summary(Testset1)$coef[1] 
SC2007.se.t9  = summary(Testset1)$coef[2] 

#SC2007 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t11 = summary(Testset1)$coef[1] 
SC2007.se.t11  = summary(Testset1)$coef[2] 

#SC2007 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t13 = summary(Testset1)$coef[1] 
SC2007.se.t13  = summary(Testset1)$coef[2] 

#SC2007 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t15 = summary(Testset1)$coef[1] 
SC2007.se.t15  = summary(Testset1)$coef[2] 

#SC2007 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2007.r.t18 = summary(Testset1)$coef[1] 
SC2007.se.t18  = summary(Testset1)$coef[2]

#SC2007 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="SC2007",]
start = cult.start[cult.start$isolate_id=="SC2007",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.2), trace=TRUE)  #WARNINGS
summary(Testset1)
# recording summary coefs #
SC2007.r.t21 = summary(Testset1)$coef[1] 
SC2007.se.t21  = summary(Testset1)$coef[2]

# Logging SC2007 isolate #
#setting temp indataframe 
#temp = c(1,3,5,7,9,11,13,18,21); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(SC2007.r.t1, SC2007.r.t3,SC2007.r.t5, SC2007.r.t7, SC2007.r.t9, SC2007.r.t11, SC2007.r.t13, SC2007.r.t15, SC2007.r.t18, SC2007.r.t21)
se = c(SC2007.se.t1, SC2007.se.t3,SC2007.se.t5, SC2007.se.t7, SC2007.se.t9,SC2007.se.t11, SC2007.se.t13,SC2007.se.t15, SC2007.se.t18, SC2007.se.t21 )
SC2007.r.se.tall = data.frame(temp, r, se)
SC2007.r.se.tall$isolate_id = "SC2007"


#### WH2008 temps 1,3, 5,7, 9,11, 13, 15, 18,21 ####
# WH2008 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t1 = summary(Testset1)$coef[1] 
WH2008.se.t1  = summary(Testset1)$coef[2]

# WH2008 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t3 = summary(Testset1)$coef[1] 
WH2008.se.t3  = summary(Testset1)$coef[2]

# WH2008 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t5 = summary(Testset1)$coef[1] 
WH2008.se.t5 = summary(Testset1)$coef[2]

# WH2008 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t7 = summary(Testset1)$coef[1] 
WH2008.se.t7 = summary(Testset1)$coef[2]

# WH2008 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t9 = summary(Testset1)$coef[1] 
WH2008.se.t9 = summary(Testset1)$coef[2]

# WH2008 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t11 = summary(Testset1)$coef[1] 
WH2008.se.t11 = summary(Testset1)$coef[2]

# WH2008 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t13 = summary(Testset1)$coef[1] 
WH2008.se.t13 = summary(Testset1)$coef[2]

# WH2008 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2008.r.t15 = summary(Testset1)$coef[1] 
WH2008.se.t15 = summary(Testset1)$coef[2]

# WH2008 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.225), trace=TRUE)  # WARNINGS
summary(Testset1)
# recording summary coefs #
WH2008.r.t18 = summary(Testset1)$coef[1] 
WH2008.se.t18 = summary(Testset1)$coef[2]

# WH2008 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="WH2008",]
start = cult.start[cult.start$isolate_id=="WH2008",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
WH2008.r.t21 = summary(Testset1)$coef[1] 
WH2008.se.t21 = summary(Testset1)$coef[2]

# Logging WH2008 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(WH2008.r.t1, WH2008.r.t3,WH2008.r.t5, WH2008.r.t7, WH2008.r.t9,WH2008.r.t11, WH2008.r.t13,WH2008.r.t15, WH2008.r.t18, WH2008.r.t21)
se = c(WH2008.se.t1, WH2008.se.t3,WH2008.se.t5, WH2008.se.t7, WH2008.se.t9,WH2008.se.t11, WH2008.se.t13,WH2008.se.t15, WH2008.se.t18, WH2008.se.t21 )
WH2008.r.se.tall = data.frame(temp, r, se)
WH2008.r.se.tall$isolate_id = "WH2008"

#### ATCC temps 1, 5, 7, 9, 11, 13, 15, 18, 21 ####
# ATCC T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t1 = summary(Testset1)$coef[1] 
ATCC.se.t1 = summary(Testset1)$coef[2]

# ATCC T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t3 = summary(Testset1)$coef[1] 
ATCC.se.t3 = summary(Testset1)$coef[2]

# ATCC T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t5 = summary(Testset1)$coef[1] 
ATCC.se.t5 = summary(Testset1)$coef[2]

# ATCC T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t7 = summary(Testset1)$coef[1] 
ATCC.se.t7 = summary(Testset1)$coef[2]

# ATCC T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t9 = summary(Testset1)$coef[1] 
ATCC.se.t9 = summary(Testset1)$coef[2]

# ATCC T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t11 = summary(Testset1)$coef[1] 
ATCC.se.t11 = summary(Testset1)$coef[2]

# ATCC T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t13 = summary(Testset1)$coef[1] 
ATCC.se.t13 = summary(Testset1)$coef[2]

# ATCC T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t15 = summary(Testset1)$coef[1] 
ATCC.se.t15 = summary(Testset1)$coef[2]

# ATCC T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
ATCC.r.t18 = summary(Testset1)$coef[1] 
ATCC.se.t18 = summary(Testset1)$coef[2]

# ATCC T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="ATCC",]
start = cult.start[cult.start$isolate_id=="ATCC",]
N0 = start$gdL
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.01), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
ATCC.r.t21 = summary(Testset1)$coef[1] 
ATCC.se.t21 = summary(Testset1)$coef[2]

# Logging ATCC isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(ATCC.r.t1, ATCC.r.t3, ATCC.r.t5, ATCC.r.t7, ATCC.r.t9, ATCC.r.t11, ATCC.r.t13, ATCC.r.t15, ATCC.r.t18, ATCC.r.t21)
se = c(ATCC.se.t1, ATCC.se.t3, ATCC.se.t5, ATCC.se.t7, ATCC.se.t9, ATCC.se.t11, ATCC.se.t13, ATCC.se.t15, ATCC.se.t18, ATCC.se.t21 )
ATCC.r.se.tall = data.frame(temp, r, se)
ATCC.r.se.tall$isolate_id = "ATCC"


#### FS2009 temps 1, 3, 5, 7, 9, 11, 13, 15, 18, 21 ####
# FS2009 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t1 = summary(Testset1)$coef[1] 
FS2009.se.t1 = summary(Testset1)$coef[2]

# FS2009 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t3 = summary(Testset1)$coef[1] 
FS2009.se.t3 = summary(Testset1)$coef[2]

# FS2009 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t5 = summary(Testset1)$coef[1] 
FS2009.se.t5 = summary(Testset1)$coef[2]

# FS2009 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t7 = summary(Testset1)$coef[1] 
FS2009.se.t7 = summary(Testset1)$coef[2]

# FS2009 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t9 = summary(Testset1)$coef[1] 
FS2009.se.t9 = summary(Testset1)$coef[2]

# FS2009 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t11 = summary(Testset1)$coef[1] 
FS2009.se.t11 = summary(Testset1)$coef[2]

# FS2009 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t13 = summary(Testset1)$coef[1] 
FS2009.se.t13 = summary(Testset1)$coef[2]

# FS2009 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2009.r.t15 = summary(Testset1)$coef[1] 
FS2009.se.t15 = summary(Testset1)$coef[2]


# FS2009 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.48), trace=TRUE)  # WARNINGS
summary(Testset1)
# recording summary coefs #
FS2009.r.t18 = summary(Testset1)$coef[1] 
FS2009.se.t18 = summary(Testset1)$coef[2]

# FS2009 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="FS2009",]
start = cult.start[cult.start$isolate_id=="FS2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.1), trace=TRUE)  #WARNINGS
summary(Testset1)
# recording summary coefs #
FS2009.r.t21 = summary(Testset1)$coef[1] 
FS2009.se.t21 = summary(Testset1)$coef[2]

# Logging FS2009 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(FS2009.r.t1, FS2009.r.t3,FS2009.r.t5, FS2009.r.t7, FS2009.r.t9,FS2009.r.t11, FS2009.r.t13,FS2009.r.t15, FS2009.r.t18, FS2009.r.t21)
se = c(FS2009.se.t1, FS2009.se.t3,FS2009.se.t5, FS2009.se.t7, FS2009.se.t9,FS2009.se.t11, FS2009.se.t13,FS2009.se.t15, FS2009.se.t18, FS2009.se.t21 )
FS2009.r.se.tall = data.frame(temp, r, se)
FS2009.r.se.tall$isolate_id = "FS2009"


#### BM2009 temps 1, 3, 5, 7, 9, 11, 13, 15, 18 ####
# BM2009 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BM2009.r.t1 = summary(Testset1)$coef[1] 
BM2009.se.t1 = summary(Testset1)$coef[2]

# BM2009 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BM2009.r.t3 = summary(Testset1)$coef[1] 
BM2009.se.t3 = summary(Testset1)$coef[2]

# BM2009 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.45), trace=TRUE)  #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2009.r.t5 = summary(Testset1)$coef[1] 
BM2009.se.t5 = summary(Testset1)$coef[2]

# BM2009 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
BM2009.r.t7 = summary(Testset1)$coef[1] 
BM2009.se.t7 = summary(Testset1)$coef[2]

# BM2009 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
BM2009.r.t9 = summary(Testset1)$coef[1] 
BM2009.se.t9 = summary(Testset1)$coef[2]


# BM2009 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
BM2009.r.t11 = summary(Testset1)$coef[1] 
BM2009.se.t11 = summary(Testset1)$coef[2]

# BM2009 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2009.r.t13 = summary(Testset1)$coef[1] 
BM2009.se.t13 = summary(Testset1)$coef[2]

# BM2009 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2009.r.t15 = summary(Testset1)$coef[1] 
BM2009.se.t15 = summary(Testset1)$coef[2]

# BM2009 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)  #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2009.r.t18 = summary(Testset1)$coef[1] 
BM2009.se.t18 = summary(Testset1)$coef[2]

# BM2009 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="BM2009",]
start = cult.start[cult.start$isolate_id=="BM2009",]
N0 = start$gdL

Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.4), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2009.r.t21 = summary(Testset1)$coef[1] 
BM2009.se.t21 = summary(Testset1)$coef[2]

# Logging BM2009 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(BM2009.r.t1, BM2009.r.t3,BM2009.r.t5, BM2009.r.t7, BM2009.r.t9,BM2009.r.t11, BM2009.r.t13,BM2009.r.t15, BM2009.r.t18, BM2009.r.t21)
se = c(BM2009.se.t1, BM2009.se.t3,BM2009.se.t5, BM2009.se.t7, BM2009.se.t9,BM2009.se.t11, BM2009.se.t13,BM2009.se.t15, BM2009.se.t18, BM2009.se.t21 )
BM2009.r.se.tall = data.frame(temp, r, se)
BM2009.r.se.tall$isolate_id = "BM2009"

#### TY2012 temps 1, 5, 9, 13, 18 ####
# TY2012 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  # ALL NA'S
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t1 = summary(Testset1)$coef[1]   
TY2012.se.t1 = summary(Testset1)$coef[2]

# TY2012 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  # ALL NA'S
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)),  
              data=Data,
              start=list(r=0.3), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t3 = summary(Testset1)$coef[1]   
TY2012.se.t3 = summary(Testset1)$coef[2]

# TY2012 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  # ALL NA'S
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t5 = summary(Testset1)$coef[1]   
TY2012.se.t5 = summary(Testset1)$coef[2]

# TY2012 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)  # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t7 = summary(Testset1)$coef[1]   
TY2012.se.t7 = summary(Testset1)$coef[2]

# TY2012 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  # ALL NA'S
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t9 = summary(Testset1)$coef[1]   
TY2012.se.t9 = summary(Testset1)$coef[2]

# TY2012 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t11 = summary(Testset1)$coef[1]   
TY2012.se.t11 = summary(Testset1)$coef[2]

# TY2012 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  # ALL NA'S
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t13 = summary(Testset1)$coef[1]   
TY2012.se.t13 = summary(Testset1)$coef[2]

# TY2012 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL 
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) # WARNINGS # ERROR
summary(Testset1)
# recording summary coefs #
TY2012.r.t15 = summary(Testset1)$coef[1]   
TY2012.se.t15 = summary(Testset1)$coef[2]


# TY2012 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  #ALL NAS
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
TY2012.r.t18 = summary(Testset1)$coef[1]   
TY2012.se.t18 = summary(Testset1)$coef[2]

# TY2012 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="TY2012",]
start = cult.start[cult.start$isolate_id=="TY2012",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)  # WARNINGS 
summary(Testset1)
# recording summary coefs #
TY2012.r.t21 = summary(Testset1)$coef[1]   
TY2012.se.t21 = summary(Testset1)$coef[2]

# Logging TY2012 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(TY2012.r.t1, TY2012.r.t3,TY2012.r.t5, TY2012.r.t7, TY2012.r.t9, TY2012.r.t11, TY2012.r.t13, TY2012.r.t15, TY2012.r.t18, TY2012.r.t21)
se = c(TY2012.se.t1, TY2012.se.t3, TY2012.se.t5, TY2012.se.t7, TY2012.se.t9, TY2012.se.t11, TY2012.se.t13, TY2012.se.t15, TY2012.se.t18, TY2012.se.t21 )
TY2012.r.se.tall = data.frame(temp, r, se)
TY2012.r.se.tall$isolate_id = "TY2012"

#### PY2019 temps 1, 5, 9, 13, 18 ####
# PY2019 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.47), trace=TRUE) #WARNING
summary(Testset1)
# recording summary coefs #
PY2019.r.t1 = summary(Testset1)$coef[1]   
PY2019.se.t1 = summary(Testset1)$coef[2]

# PY2019 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.35), trace=TRUE) #WARNING
summary(Testset1)
# recording summary coefs #
PY2019.r.t3 = summary(Testset1)$coef[1]   
PY2019.se.t3 = summary(Testset1)$coef[2]


# PY2019 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
PY2019.r.t5 = summary(Testset1)$coef[1]   
PY2019.se.t5 = summary(Testset1)$coef[2]

# PY2019 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
PY2019.r.t7 = summary(Testset1)$coef[1]   
PY2019.se.t7 = summary(Testset1)$coef[2]

# PY2019 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
PY2019.r.t9 = summary(Testset1)$coef[1]   
PY2019.se.t9 = summary(Testset1)$coef[2]

# PY2019 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
PY2019.r.t11 = summary(Testset1)$coef[1]   
PY2019.se.t11 = summary(Testset1)$coef[2]

# PY2019 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
PY2019.r.t13 = summary(Testset1)$coef[1]   
PY2019.se.t13 = summary(Testset1)$coef[2]

# PY2019 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
PY2019.r.t15 = summary(Testset1)$coef[1]   
PY2019.se.t15 = summary(Testset1)$coef[2]

# PY2019 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
PY2019.r.t18 = summary(Testset1)$coef[1]   
PY2019.se.t18 = summary(Testset1)$coef[2]

# PY2019 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="PY2019",]
start = cult.start[cult.start$isolate_id=="PY2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.21), trace=TRUE) #WARNING
summary(Testset1)
# recording summary coefs #
PY2019.r.t21 = summary(Testset1)$coef[1]   
PY2019.se.t21 = summary(Testset1)$coef[2]

# Logging PY2019 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(PY2019.r.t1, PY2019.r.t3, PY2019.r.t5, PY2019.r.t7, PY2019.r.t9, PY2019.r.t11, PY2019.r.t13, PY2019.r.t15, PY2019.r.t18, PY2019.r.t21)
se = c(PY2019.se.t1, PY2019.se.t3, PY2019.se.t5, PY2019.se.t7, PY2019.se.t9, PY2019.se.t11, PY2019.se.t13, PY2019.se.t15, PY2019.se.t18, PY2019.se.t21 )
PY2019.r.se.tall = data.frame(temp, r, se)
PY2019.r.se.tall$isolate_id = "PY2019"

#### NR2019 temps 1, 5, 9, 13, 18 ####
# NR2019 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t1 = summary(Testset1)$coef[1]   
NR2019.se.t1 = summary(Testset1)$coef[2]

# NR2019 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t3 = summary(Testset1)$coef[1]   
NR2019.se.t3 = summary(Testset1)$coef[2]

# NR2019 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t5 = summary(Testset1)$coef[1]   
NR2019.se.t5 = summary(Testset1)$coef[2]

# NR2019 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t7 = summary(Testset1)$coef[1]   
NR2019.se.t7 = summary(Testset1)$coef[2]

# NR2019 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.8), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t9 = summary(Testset1)$coef[1]   
NR2019.se.t9 = summary(Testset1)$coef[2]

# NR2019 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.8), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t11 = summary(Testset1)$coef[1]   
NR2019.se.t11 = summary(Testset1)$coef[2]


# NR2019 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t13 = summary(Testset1)$coef[1]   
NR2019.se.t13 = summary(Testset1)$coef[2]

# NR2019 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t15 = summary(Testset1)$coef[1]   
NR2019.se.t15 = summary(Testset1)$coef[2]

# NR2019 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.7), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
NR2019.r.t18 = summary(Testset1)$coef[1]   
NR2019.se.t18 = summary(Testset1)$coef[2]

# NR2019 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="NR2019",]
start = cult.start[cult.start$isolate_id=="NR2019",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.35), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
NR2019.r.t21 = summary(Testset1)$coef[1]   
NR2019.se.t21 = summary(Testset1)$coef[2]

# Logging NR2019 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(NR2019.r.t1, NR2019.r.t3, NR2019.r.t5, NR2019.r.t7, NR2019.r.t9, NR2019.r.t11, NR2019.r.t13, NR2019.r.t15, NR2019.r.t18, NR2019.r.t21)
se = c(NR2019.se.t1, NR2019.se.t3, NR2019.se.t5, NR2019.se.t7, NR2019.se.t9, NR2019.se.t11, NR2019.se.t13, NR2019.se.t15, NR2019.se.t18, NR2019.se.t21 )
NR2019.r.se.tall = data.frame(temp, r, se)
NR2019.r.se.tall$isolate_id = "NR2019"

#### BM2019 temps 1, 5, 9, 13, 18 ####
# BM2019 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL 
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2019.r.t1 = summary(Testset1)$coef[1]   
BM2019.se.t1 = summary(Testset1)$coef[2]

# BM2019 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2019.r.t3 = summary(Testset1)$coef[1]   
BM2019.se.t3 = summary(Testset1)$coef[2]

# BM2019 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.42), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2019.r.t5 = summary(Testset1)$coef[1]   
BM2019.se.t5 = summary(Testset1)$coef[2]

# BM2019 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2019.r.t7 = summary(Testset1)$coef[1]   
BM2019.se.t7 = summary(Testset1)$coef[2]

# BM2019 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2019.r.t9 = summary(Testset1)$coef[1]   
BM2019.se.t9 = summary(Testset1)$coef[2]

# BM2019 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2019.r.t11 = summary(Testset1)$coef[1]   
BM2019.se.t11 = summary(Testset1)$coef[2]


# BM2019 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2019.r.t13 = summary(Testset1)$coef[1]   
BM2019.se.t13 = summary(Testset1)$coef[2]

# BM2019 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BM2019.r.t15 = summary(Testset1)$coef[1]   
BM2019.se.t15 = summary(Testset1)$coef[2]


# BM2019 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.005), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BM2019.r.t18 = summary(Testset1)$coef[1]   
BM2019.se.t18 = summary(Testset1)$coef[2]

# BM2019 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="BM2019",]
start = cult.start[cult.start$isolate_id=="BM2019",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.1), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BM2019.r.t21 = summary(Testset1)$coef[1]   
BM2019.se.t21 = summary(Testset1)$coef[2]

# Logging BM2019 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(BM2019.r.t1, BM2019.r.t3, BM2019.r.t5, BM2019.r.t7, BM2019.r.t9, BM2019.r.t11, BM2019.r.t13, BM2019.r.t15, BM2019.r.t18, BM2019.r.t21)
se = c(BM2019.se.t1, BM2019.se.t3, BM2019.se.t5, BM2019.se.t7, BM2019.se.t9, BM2019.se.t11, BM2019.se.t13, BM2019.se.t15, BM2019.se.t18, BM2019.se.t21 )
BM2019.r.se.tall = data.frame(temp, r, se)
BM2019.r.se.tall$isolate_id = "BM2019"


#### WL2020 temps 1, 5, 9, 13, 18 ####
# WL2020 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0045), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t1 = summary(Testset1)$coef[1]   
WL2020.se.t1 = summary(Testset1)$coef[2]

# WL2020 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0045), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t3 = summary(Testset1)$coef[1]   
WL2020.se.t3 = summary(Testset1)$coef[2]


# WL2020 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0045), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t5 = summary(Testset1)$coef[1]   
WL2020.se.t5 = summary(Testset1)$coef[2]

# WL2020 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0045), trace=TRUE)  #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t7 = summary(Testset1)$coef[1]   
WL2020.se.t7 = summary(Testset1)$coef[2]


# WL2020 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t9 = summary(Testset1)$coef[1]   
WL2020.se.t9 = summary(Testset1)$coef[2]

# WL2020 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE)  #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t11 = summary(Testset1)$coef[1]   
WL2020.se.t11 = summary(Testset1)$coef[2]

# WL2020 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.005), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t13 = summary(Testset1)$coef[1]   
WL2020.se.t13 = summary(Testset1)$coef[2]


# WL2020 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.005), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t15 = summary(Testset1)$coef[1]   
WL2020.se.t15 = summary(Testset1)$coef[2]


# WL2020 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t18 = summary(Testset1)$coef[1]   
WL2020.se.t18 = summary(Testset1)$coef[2]

# WL2020 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="WL2020",]
start = cult.start[cult.start$isolate_id=="WL2020",]
N0 = start$gdL  # ALL NAs
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE)  #WARNINGS #ERROR
summary(Testset1)
# recording summary coefs #
WL2020.r.t21 = summary(Testset1)$coef[1]   
WL2020.se.t21 = summary(Testset1)$coef[2]

# Logging WL2020 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(WL2020.r.t1, WL2020.r.t3, WL2020.r.t5, WL2020.r.t7, WL2020.r.t9, WL2020.r.t11, WL2020.r.t13, WL2020.r.t15, WL2020.r.t18, WL2020.r.t21)
se = c(WL2020.se.t1, WL2020.se.t3, WL2020.se.t5, WL2020.se.t7, WL2020.se.t9, WL2020.se.t11, WL2020.se.t13, WL2020.se.t15, WL2020.se.t18, WL2020.se.t21 )
WL2020.r.se.tall = data.frame(temp, r, se)
WL2020.r.se.tall$isolate_id = "WL2020"


#### WH2020 temps 1, 5, 9, 13, 18 ####
# WH2020 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2020.r.t1 = summary(Testset1)$coef[1]   
WH2020.se.t1 = summary(Testset1)$coef[2]

# WH2020 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
WH2020.r.t3 = summary(Testset1)$coef[1]   
WH2020.se.t3 = summary(Testset1)$coef[2]

# WH2020 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t5 = summary(Testset1)$coef[1]   
WH2020.se.t5 = summary(Testset1)$coef[2]

# WH2020 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t7 = summary(Testset1)$coef[1]   
WH2020.se.t7 = summary(Testset1)$coef[2]

# WH2020 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t9 = summary(Testset1)$coef[1]   
WH2020.se.t9 = summary(Testset1)$coef[2]

# WH2020 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t11 = summary(Testset1)$coef[1]   
WH2020.se.t11 = summary(Testset1)$coef[2]

# WH2020 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t13 = summary(Testset1)$coef[1]   
WH2020.se.t13 = summary(Testset1)$coef[2]

# WH2020 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t15 = summary(Testset1)$coef[1]   
WH2020.se.t15 = summary(Testset1)$coef[2]


# WH2020 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0001), trace=TRUE) # WARNINGS
summary(Testset1)
# recording summary coefs #
WH2020.r.t18 = summary(Testset1)$coef[1]   
WH2020.se.t18 = summary(Testset1)$coef[2]

# WH2020 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="WH2020",]
start = cult.start[cult.start$isolate_id=="WH2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.0001), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
WH2020.r.t21 = summary(Testset1)$coef[1]   
WH2020.se.t21 = summary(Testset1)$coef[2]

# Logging WH2020 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(WH2020.r.t1, WH2020.r.t3, WH2020.r.t5, WH2020.r.t7, WH2020.r.t9, WH2020.r.t11, WH2020.r.t13, WH2020.r.t15, WH2020.r.t18, WH2020.r.t21)
se = c(WH2020.se.t1, WH2020.se.t3, WH2020.se.t5, WH2020.se.t7, WH2020.se.t9, WH2020.se.t11, WH2020.se.t13, WH2020.se.t15, WH2020.se.t18, WH2020.se.t21 )
WH2020.r.se.tall = data.frame(temp, r, se)
WH2020.r.se.tall$isolate_id = "WH2020"


#### FS2020 temps 1, 5, 9, 13, 18 ####
# FS2020 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2020.r.t1 = summary(Testset1)$coef[1]   
FS2020.se.t1 = summary(Testset1)$coef[2]

# FS2020 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2020.r.t3 = summary(Testset1)$coef[1]   
FS2020.se.t3 = summary(Testset1)$coef[2]

# FS2020 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2020.r.t5 = summary(Testset1)$coef[1]   
FS2020.se.t5 = summary(Testset1)$coef[2]

# FS2020 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
FS2020.r.t7 = summary(Testset1)$coef[1]   
FS2020.se.t7 = summary(Testset1)$coef[2]

# FS2020 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE) # WARNINGS
summary(Testset1)
# recording summary coefs #
FS2020.r.t9 = summary(Testset1)$coef[1]   
FS2020.se.t9 = summary(Testset1)$coef[2]

# FS2020 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
FS2020.r.t11 = summary(Testset1)$coef[1]   
FS2020.se.t11 = summary(Testset1)$coef[2]


# FS2020 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
FS2020.r.t13 = summary(Testset1)$coef[1]   
FS2020.se.t13 = summary(Testset1)$coef[2]

# FS2020 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
FS2020.r.t15 = summary(Testset1)$coef[1]   
FS2020.se.t15 = summary(Testset1)$coef[2]

# FS2020 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) # WARNINGS
summary(Testset1)
# recording summary coefs #
FS2020.r.t18 = summary(Testset1)$coef[1]   
FS2020.se.t18 = summary(Testset1)$coef[2]

# FS2020 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="FS2020",]
start = cult.start[cult.start$isolate_id=="FS2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
FS2020.r.t21 = summary(Testset1)$coef[1]   
FS2020.se.t21 = summary(Testset1)$coef[2]

# Logging FS2020 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(FS2020.r.t1, FS2020.r.t3, FS2020.r.t5, FS2020.r.t7, FS2020.r.t9, FS2020.r.t11, FS2020.r.t13, FS2020.r.t15, FS2020.r.t18, FS2020.r.t21)
se = c(FS2020.se.t1, FS2020.se.t3, FS2020.se.t5, FS2020.se.t7, FS2020.se.t9, FS2020.se.t11, FS2020.se.t13, FS2020.se.t15, FS2020.se.t18, FS2020.se.t21 )
FS2020.r.se.tall = data.frame(temp, r, se)
FS2020.r.se.tall$isolate_id = "FS2020"


#### SC2020 temps 1, 5, 9, 13, 18 ####
# SC2020 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t1 = summary(Testset1)$coef[1]   
SC2020.se.t1 = summary(Testset1)$coef[2]

# SC2020 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t3 = summary(Testset1)$coef[1]   
SC2020.se.t3 = summary(Testset1)$coef[2]


# SC2020 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t5 = summary(Testset1)$coef[1]   
SC2020.se.t5 = summary(Testset1)$coef[2]

# SC2020 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t7 = summary(Testset1)$coef[1]   
SC2020.se.t7 = summary(Testset1)$coef[2]


# SC2020 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t9 = summary(Testset1)$coef[1]   
SC2020.se.t9 = summary(Testset1)$coef[2]

# SC2020 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t11 = summary(Testset1)$coef[1]   
SC2020.se.t11 = summary(Testset1)$coef[2]

# SC2020 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t13 = summary(Testset1)$coef[1]   
SC2020.se.t13 = summary(Testset1)$coef[2]

# SC2020 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
SC2020.r.t15 = summary(Testset1)$coef[1]   
SC2020.se.t15 = summary(Testset1)$coef[2]

# SC2020 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.1), trace=TRUE) # WARNINGS
summary(Testset1)
# recording summary coefs #
SC2020.r.t18 = summary(Testset1)$coef[1]   
SC2020.se.t18 = summary(Testset1)$coef[2]

# SC2020 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="SC2020",]
start = cult.start[cult.start$isolate_id=="SC2020",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.56), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
SC2020.r.t21 = summary(Testset1)$coef[1]   
SC2020.se.t21 = summary(Testset1)$coef[2]

# Logging SC2020 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(SC2020.r.t1, SC2020.r.t3, SC2020.r.t5, SC2020.r.t7, SC2020.r.t9, SC2020.r.t11, SC2020.r.t13, SC2020.r.t15, SC2020.r.t18, SC2020.r.t21)
se = c(SC2020.se.t1, SC2020.se.t3, SC2020.se.t5, SC2020.se.t7, SC2020.se.t9, SC2020.se.t11, SC2020.se.t13, SC2020.se.t15, SC2020.se.t18, SC2020.se.t21 )
SC2020.r.se.tall = data.frame(temp, r, se)
SC2020.r.se.tall$isolate_id = "SC2020"


#### BB2020 temps 1, 5, 9, 13, 18 ####
# BB2020 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t1 = summary(Testset1)$coef[1]   
BB2020.se.t1 = summary(Testset1)$coef[2]

# BB2020 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BB2020.r.t3 = summary(Testset1)$coef[1]   
BB2020.se.t3 = summary(Testset1)$coef[2]

# BB2020 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t5 = summary(Testset1)$coef[1]   
BB2020.se.t5 = summary(Testset1)$coef[2]

# BB2020 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t7 = summary(Testset1)$coef[1]   
BB2020.se.t7 = summary(Testset1)$coef[2]


# BB2020 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t9 = summary(Testset1)$coef[1]   
BB2020.se.t9 = summary(Testset1)$coef[2]

# BB2020 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t11 = summary(Testset1)$coef[1]   
BB2020.se.t11 = summary(Testset1)$coef[2]

# BB2020 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t13 = summary(Testset1)$coef[1]   
BB2020.se.t13 = summary(Testset1)$coef[2]

# BB2020 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BB2020.r.t15 = summary(Testset1)$coef[1]   
BB2020.se.t15 = summary(Testset1)$coef[2]


# BB2020 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE) # WARNING
summary(Testset1)
# recording summary coefs #
BB2020.r.t18 = summary(Testset1)$coef[1]   
BB2020.se.t18 = summary(Testset1)$coef[2]

# BB2020 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="BB2020",]
start = cult.start[cult.start$isolate_id=="BB2020",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.55), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
BB2020.r.t21 = summary(Testset1)$coef[1]   
BB2020.se.t21 = summary(Testset1)$coef[2]

# Logging BB2020 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(BB2020.r.t1, BB2020.r.t3, BB2020.r.t5, BB2020.r.t7, BB2020.r.t9, BB2020.r.t11, BB2020.r.t13, BB2020.r.t15, BB2020.r.t18, BB2020.r.t21)
se = c(BB2020.se.t1, BB2020.se.t3, BB2020.se.t5, BB2020.se.t7, BB2020.se.t9, BB2020.se.t11, BB2020.se.t13, BB2020.se.t15, BB2020.se.t18, BB2020.se.t21 )
BB2020.r.se.tall = data.frame(temp, r, se)
BB2020.r.se.tall$isolate_id = "BB2020"


#### BC2021 temps 1, 5, 9, 13, 18 ####
# BC2021 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t1 = summary(Testset1)$coef[1]   
BC2021.se.t1 = summary(Testset1)$coef[2]

# BC2021 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t3 = summary(Testset1)$coef[1]   
BC2021.se.t3 = summary(Testset1)$coef[2]


# BC2021 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BC2021.r.t5 = summary(Testset1)$coef[1]   
BC2021.se.t5 = summary(Testset1)$coef[2]

# BC2021 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t7 = summary(Testset1)$coef[1]   
BC2021.se.t7 = summary(Testset1)$coef[2]


# BC2021 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
BC2021.r.t9 = summary(Testset1)$coef[1]   
BC2021.se.t9 = summary(Testset1)$coef[2]

# BC2021 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t11 = summary(Testset1)$coef[1]   
BC2021.se.t11 = summary(Testset1)$coef[2]


# BC2021 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
BC2021.r.t13 = summary(Testset1)$coef[1]   
BC2021.se.t13 = summary(Testset1)$coef[2]

# BC2021 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t15 = summary(Testset1)$coef[1]   
BC2021.se.t15 = summary(Testset1)$coef[2]


# BC2021 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE)
summary(Testset1)
# recording summary coefs #
BC2021.r.t18 = summary(Testset1)$coef[1]   
BC2021.se.t18 = summary(Testset1)$coef[2]

# BC2021 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="BC2021",]
start = cult.start[cult.start$isolate_id=="BC2021",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
BC2021.r.t21 = summary(Testset1)$coef[1]   
BC2021.se.t21 = summary(Testset1)$coef[2]


# Logging BC2021 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(BC2021.r.t1, BC2021.r.t3, BC2021.r.t5, BC2021.r.t7, BC2021.r.t9, BC2021.r.t11, BC2021.r.t13, BC2021.r.t15, BC2021.r.t18, BC2021.r.t21)
se = c(BC2021.se.t1, BC2021.se.t3, BC2021.se.t5, BC2021.se.t7, BC2021.se.t9, BC2021.se.t11, BC2021.se.t13, BC2021.se.t15, BC2021.se.t18, BC2021.se.t21 )
BC2021.r.se.tall = data.frame(temp, r, se)
BC2021.r.se.tall$isolate_id = "BC2021"

#### CNBEIJ temps 1, 5, 9, 13, 18 ####
# CNBEIJ T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t1 = summary(Testset1)$coef[1]   
CNBEIJ.se.t1 = summary(Testset1)$coef[2]

# CNBEIJ T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.1), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t3 = summary(Testset1)$coef[1]   
CNBEIJ.se.t3 = summary(Testset1)$coef[2]


# CNBEIJ T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t5 = summary(Testset1)$coef[1]   
CNBEIJ.se.t5 = summary(Testset1)$coef[2]

# CNBEIJ T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t7 = summary(Testset1)$coef[1]   
CNBEIJ.se.t7 = summary(Testset1)$coef[2]


# CNBEIJ T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t9 = summary(Testset1)$coef[1]   
CNBEIJ.se.t9 = summary(Testset1)$coef[2]

# CNBEIJ T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t11 = summary(Testset1)$coef[1]   
CNBEIJ.se.t11 = summary(Testset1)$coef[2]


# CNBEIJ T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t13 = summary(Testset1)$coef[1]   
CNBEIJ.se.t13 = summary(Testset1)$coef[2]

# CNBEIJ T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.07), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t15 = summary(Testset1)$coef[1]   
CNBEIJ.se.t15 = summary(Testset1)$coef[2]


# CNBEIJ T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t18 = summary(Testset1)$coef[1]   
CNBEIJ.se.t18 = summary(Testset1)$coef[2]

# CNBEIJ T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="CNBEIJ",]
start = cult.start[cult.start$isolate_id=="CNBEIJ",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-1), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNBEIJ.r.t21 = summary(Testset1)$coef[1]   
CNBEIJ.se.t21 = summary(Testset1)$coef[2]


# Logging CNBEIJ isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(CNBEIJ.r.t1, CNBEIJ.r.t3, CNBEIJ.r.t5, CNBEIJ.r.t7, CNBEIJ.r.t9, CNBEIJ.r.t11, CNBEIJ.r.t13, CNBEIJ.r.t15, CNBEIJ.r.t18, CNBEIJ.r.t21)
se = c(CNBEIJ.se.t1, CNBEIJ.se.t3, CNBEIJ.se.t5, CNBEIJ.se.t7, CNBEIJ.se.t9, CNBEIJ.se.t11, CNBEIJ.se.t13, CNBEIJ.se.t15, CNBEIJ.se.t18, CNBEIJ.se.t21 )
CNBEIJ.r.se.tall = data.frame(temp, r, se)
CNBEIJ.r.se.tall$isolate_id = "CNBEIJ"


#### CNGEZI temps 1, 5, 9, 13, 18 ####
# CNGEZI T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t1 = summary(Testset1)$coef[1]   
CNGEZI.se.t1 = summary(Testset1)$coef[2] 

# CNGEZI T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t3 = summary(Testset1)$coef[1]   
CNGEZI.se.t3 = summary(Testset1)$coef[2] 

# CNGEZI T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t5 = summary(Testset1)$coef[1]   
CNGEZI.se.t5 = summary(Testset1)$coef[2]

# CNGEZI T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.04), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t7 = summary(Testset1)$coef[1]   
CNGEZI.se.t7 = summary(Testset1)$coef[2] 

# CNGEZI T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t9 = summary(Testset1)$coef[1]   
CNGEZI.se.t9 = summary(Testset1)$coef[2]

# CNGEZI T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.05), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t11 = summary(Testset1)$coef[1]   
CNGEZI.se.t11 = summary(Testset1)$coef[2] 

# CNGEZI T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t13 = summary(Testset1)$coef[1]   
CNGEZI.se.t13 = summary(Testset1)$coef[2]

# CNGEZI T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.2), trace=TRUE)
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t15 = summary(Testset1)$coef[1]   
CNGEZI.se.t15 = summary(Testset1)$coef[2] 

# CNGEZI T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.31), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t18 = summary(Testset1)$coef[1]   
CNGEZI.se.t18 = summary(Testset1)$coef[2]

# CNGEZI T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="CNGEZI",]
start = cult.start[cult.start$isolate_id=="CNGEZI",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.461), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
CNGEZI.r.t21 = summary(Testset1)$coef[1]   
CNGEZI.se.t21 = summary(Testset1)$coef[2] 

# Logging CNGEZI isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(CNGEZI.r.t1, CNGEZI.r.t3, CNGEZI.r.t5, CNGEZI.r.t7, CNGEZI.r.t9, CNGEZI.r.t11, CNGEZI.r.t13, CNGEZI.r.t15, CNGEZI.r.t18, CNGEZI.r.t21)
se = c(CNGEZI.se.t1, CNGEZI.se.t3, CNGEZI.se.t5, CNGEZI.se.t7, CNGEZI.se.t9, CNGEZI.se.t11, CNGEZI.se.t13, CNGEZI.se.t15, CNGEZI.se.t18, CNGEZI.se.t21 )
CNGEZI.r.se.tall = data.frame(temp, r, se)
CNGEZI.r.se.tall$isolate_id = "CNGEZI"


#### HUNG temps 1, 5, 9, 13, 18 ####
# HUNG T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t1 = summary(Testset1)$coef[1]   
HUNG.se.t1 = summary(Testset1)$coef[2]

# HUNG T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t3 = summary(Testset1)$coef[1]   
HUNG.se.t3 = summary(Testset1)$coef[2]

# HUNG T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t5 = summary(Testset1)$coef[1]   
HUNG.se.t5 = summary(Testset1)$coef[2]

# HUNG T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t7 = summary(Testset1)$coef[1]   
HUNG.se.t7 = summary(Testset1)$coef[2]

# HUNG T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t9 = summary(Testset1)$coef[1]   
HUNG.se.t9 = summary(Testset1)$coef[2]

# HUNG T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t11 = summary(Testset1)$coef[1]   
HUNG.se.t11 = summary(Testset1)$coef[2]

# HUNG T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.8), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t13 = summary(Testset1)$coef[1]   
HUNG.se.t13 = summary(Testset1)$coef[2]

# HUNG T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t15 = summary(Testset1)$coef[1]   
HUNG.se.t15 = summary(Testset1)$coef[2]

# HUNG T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t18 = summary(Testset1)$coef[1]   
HUNG.se.t18 = summary(Testset1)$coef[2]

# HUNG T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="HUNG",]
start = cult.start[cult.start$isolate_id=="HUNG",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.9), trace=TRUE)
summary(Testset1)
# recording summary coefs #
HUNG.r.t21 = summary(Testset1)$coef[1]   
HUNG.se.t21 = summary(Testset1)$coef[2]

# Logging HUNG isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(HUNG.r.t1, HUNG.r.t3, HUNG.r.t5, HUNG.r.t7, HUNG.r.t9, HUNG.r.t11, HUNG.r.t13, HUNG.r.t15, HUNG.r.t18, HUNG.r.t21)
se = c(HUNG.se.t1, HUNG.se.t3, HUNG.se.t5, HUNG.se.t7, HUNG.se.t9, HUNG.se.t11, HUNG.se.t13, HUNG.se.t15, HUNG.se.t18, HUNG.se.t21 )
HUNG.r.se.tall = data.frame(temp, r, se)
HUNG.r.se.tall$isolate_id = "HUNG"

#### MONG1 temps 1, 5, 9, 13, 18 ####
# MONG1 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t1 = summary(Testset1)$coef[1]   
MONG1.se.t1 = summary(Testset1)$coef[2]

# MONG1 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t3 = summary(Testset1)$coef[1]   
MONG1.se.t3 = summary(Testset1)$coef[2]

# MONG1 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t5 = summary(Testset1)$coef[1]   
MONG1.se.t5 = summary(Testset1)$coef[2]

# MONG1 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t7 = summary(Testset1)$coef[1]   
MONG1.se.t7 = summary(Testset1)$coef[2]

# MONG1 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t9 = summary(Testset1)$coef[1]   
MONG1.se.t9 = summary(Testset1)$coef[2]

# MONG1 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t11 = summary(Testset1)$coef[1]   
MONG1.se.t11 = summary(Testset1)$coef[2]

# MONG1 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t13 = summary(Testset1)$coef[1]   
MONG1.se.t13 = summary(Testset1)$coef[2]

# MONG1 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t15 = summary(Testset1)$coef[1]   
MONG1.se.t15 = summary(Testset1)$coef[2]

# MONG1 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t18 = summary(Testset1)$coef[1]   
MONG1.se.t18 = summary(Testset1)$coef[2]

# MONG1 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="MONG1",]
start = cult.start[cult.start$isolate_id=="MONG1",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG1.r.t21 = summary(Testset1)$coef[1]   
MONG1.se.t21 = summary(Testset1)$coef[2]

# Logging MONG1 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(MONG1.r.t1, MONG1.r.t3, MONG1.r.t5, MONG1.r.t7, MONG1.r.t9, MONG1.r.t11, MONG1.r.t13, MONG1.r.t15, MONG1.r.t18, MONG1.r.t21)
se = c(MONG1.se.t1, MONG1.se.t3, MONG1.se.t5, MONG1.se.t7, MONG1.se.t9, MONG1.se.t11, MONG1.se.t13, MONG1.se.t15, MONG1.se.t18, MONG1.se.t21 )
MONG1.r.se.tall = data.frame(temp, r, se)
MONG1.r.se.tall$isolate_id = "MONG1"

#### MONG2 temps 1, 5, 9, 13, 18 ####
# MONG2 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG2.r.t1 = summary(Testset1)$coef[1]   
MONG2.se.t1 = summary(Testset1)$coef[2]

# MONG2 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.37), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
MONG2.r.t3 = summary(Testset1)$coef[1]   
MONG2.se.t3 = summary(Testset1)$coef[2]

# MONG2 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
MONG2.r.t5 = summary(Testset1)$coef[1]   
MONG2.se.t5 = summary(Testset1)$coef[2]

# MONG2 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG2.r.t7 = summary(Testset1)$coef[1]   
MONG2.se.t7 = summary(Testset1)$coef[2]

# MONG2 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE) 
summary(Testset1)
# recording summary coefs #
MONG2.r.t9 = summary(Testset1)$coef[1]   
MONG2.se.t9 = summary(Testset1)$coef[2]

# MONG2 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG2.r.t11 = summary(Testset1)$coef[1]   
MONG2.se.t11 = summary(Testset1)$coef[2]

# MONG2 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)  
summary(Testset1)
# recording summary coefs #
MONG2.r.t13 = summary(Testset1)$coef[1]   
MONG2.se.t13 = summary(Testset1)$coef[2]

# MONG2 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG2.r.t15 = summary(Testset1)$coef[1]   
MONG2.se.t15 = summary(Testset1)$coef[2]

# MONG2 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.6), trace=TRUE)  #WARNING
summary(Testset1)
# recording summary coefs #
MONG2.r.t18 = summary(Testset1)$coef[1]   
MONG2.se.t18 = summary(Testset1)$coef[2]

# MONG2 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="MONG2",]
start = cult.start[cult.start$isolate_id=="MONG2",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
MONG2.r.t21 = summary(Testset1)$coef[1]   
MONG2.se.t21 = summary(Testset1)$coef[2]

# Logging MONG2 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(MONG2.r.t1, MONG2.r.t3, MONG2.r.t5, MONG2.r.t7, MONG2.r.t9, MONG2.r.t11, MONG2.r.t13, MONG2.r.t15, MONG2.r.t18, MONG2.r.t21)
se = c(MONG2.se.t1, MONG2.se.t3, MONG2.se.t5, MONG2.se.t7, MONG2.se.t9, MONG2.se.t11, MONG2.se.t13, MONG2.se.t15, MONG2.se.t18, MONG2.se.t21 )
MONG2.r.se.tall = data.frame(temp, r, se)
MONG2.r.se.tall$isolate_id = "MONG2"


#### UK1 temps 1, 5, 9, 13, 18 ####
# UK1 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t1 = summary(Testset1)$coef[1]   
UK1.se.t1 = summary(Testset1)$coef[2]

# UK1 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t3 = summary(Testset1)$coef[1]   
UK1.se.t3 = summary(Testset1)$coef[2]

# UK1 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t5 = summary(Testset1)$coef[1]   
UK1.se.t5 = summary(Testset1)$coef[2]

# UK1 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t7 = summary(Testset1)$coef[1]   
UK1.se.t7 = summary(Testset1)$coef[2]

# UK1 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t9 = summary(Testset1)$coef[1]   
UK1.se.t9 = summary(Testset1)$coef[2]

# UK1 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t11 = summary(Testset1)$coef[1]   
UK1.se.t11 = summary(Testset1)$coef[2]

# UK1 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.7), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t13 = summary(Testset1)$coef[1]   
UK1.se.t13 = summary(Testset1)$coef[2]

# UK1 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t15 = summary(Testset1)$coef[1]   
UK1.se.t15 = summary(Testset1)$coef[2]

# UK1 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.1), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK1.r.t18 = summary(Testset1)$coef[1]   
UK1.se.t18 = summary(Testset1)$coef[2]

# UK1 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="UK1",]
start = cult.start[cult.start$isolate_id=="UK1",]
N0 = start$gdL  
##we probably(?) want to sub in this value for every single 0 
N0[is.na(N0)==T] = min(new.cult$gdL, na.rm = T)
N0
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.15), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
UK1.r.t21 = summary(Testset1)$coef[1]   
UK1.se.t21 = summary(Testset1)$coef[2]


# Logging UK1 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
r = c(UK1.r.t1, UK1.r.t3, UK1.r.t5, UK1.r.t7, UK1.r.t9, UK1.r.t11, UK1.r.t13, UK1.r.t15, UK1.r.t18, UK1.r.t21)
se = c(UK1.se.t1, UK1.se.t3, UK1.se.t5, UK1.se.t7, UK1.se.t9, UK1.se.t11, UK1.se.t13, UK1.se.t15, UK1.se.t18, UK1.se.t21 )
UK1.r.se.tall = data.frame(temp, r, se)
UK1.r.se.tall$isolate_id = "UK1"


#### UK2 temps 1, 5, 9, 13, 18 ####
# UK2 T1
Data<-new.cult[new.cult$day > 0 & new.cult$temp==1 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t1 = summary(Testset1)$coef[1]   
UK2.se.t1 = summary(Testset1)$coef[2]

# UK2 T3
Data<-new.cult[new.cult$day > 0 & new.cult$temp==3 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.285), trace=TRUE) #WARNINGS 
summary(Testset1)
# recording summary coefs #
UK2.r.t3 = summary(Testset1)$coef[1]   
UK2.se.t3 = summary(Testset1)$coef[2]

# UK2 T5
Data<-new.cult[new.cult$day > 0 & new.cult$temp==5 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t5 = summary(Testset1)$coef[1]   
UK2.se.t5 = summary(Testset1)$coef[2]

# UK2 T7
Data<-new.cult[new.cult$day > 0 & new.cult$temp==7 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t7 = summary(Testset1)$coef[1]   
UK2.se.t7 = summary(Testset1)$coef[2]

# UK2 T9
Data<-new.cult[new.cult$day > 0 & new.cult$temp==9 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t9 = summary(Testset1)$coef[1]   
UK2.se.t9 = summary(Testset1)$coef[2]

# UK2 T11
Data<-new.cult[new.cult$day > 0 & new.cult$temp==11 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t11 = summary(Testset1)$coef[1]   
UK2.se.t11 = summary(Testset1)$coef[2]

# UK2 T13
Data<-new.cult[new.cult$day > 0 & new.cult$temp==13 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.5), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t13 = summary(Testset1)$coef[1]   
UK2.se.t13 = summary(Testset1)$coef[2]

# UK2 T15
Data<-new.cult[new.cult$day > 0 & new.cult$temp==15 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.4), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t15 = summary(Testset1)$coef[1]   
UK2.se.t15 = summary(Testset1)$coef[2]

# UK2 T18
Data<-new.cult[new.cult$day > 0 & new.cult$temp==18 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=0.3), trace=TRUE)
summary(Testset1)
# recording summary coefs #
UK2.r.t18 = summary(Testset1)$coef[1]   
UK2.se.t18 = summary(Testset1)$coef[2]

# UK2 T21
Data<-new.cult[new.cult$day > 0 & new.cult$temp==21 & new.cult$isolate_id=="UK2",]
start = cult.start[cult.start$isolate_id=="UK2",]
N0 = start$gdL  
Testset1<-nls(gdL~K/(1+((K/N0)-1)*exp(-r*day)), 
              data=Data,
              start=list(r=-0.5), trace=TRUE) #WARNINGS
summary(Testset1)
# recording summary coefs #
UK2.r.t21 = summary(Testset1)$coef[1]   
UK2.se.t21 = summary(Testset1)$coef[2]


# Logging UK2 isolate #
#setting temp indataframe 
#temp = c(1,5,9,13,18); 
temp=c(unique(new.cult$temp))

#logging in R's and SE for R
#logging in R's and SE for R
r = c(UK2.r.t1, UK2.r.t3, UK2.r.t5, UK2.r.t7, UK2.r.t9, UK2.r.t11, UK2.r.t13, UK2.r.t15, UK2.r.t18, UK2.r.t21)
se = c(UK2.se.t1, UK2.se.t3, UK2.se.t5, UK2.se.t7, UK2.se.t9, UK2.se.t11, UK2.se.t13, UK2.se.t15, UK2.se.t18, UK2.se.t21 )
UK2.r.se.tall = data.frame(temp, r, se)
UK2.r.se.tall$isolate_id = "UK2"


#### BINDING ALL TOGETHER ####
temp.rse = bind_rows(HL2007.r.se.tall, SC2007.r.se.tall, WH2008.r.se.tall, ATCC.r.se.tall, FS2009.r.se.tall, BM2009.r.se.tall,
                     TY2012.r.se.tall, PY2019.r.se.tall, NR2019.r.se.tall, BM2019.r.se.tall, SC2020.r.se.tall, WH2020.r.se.tall,
                     WL2020.r.se.tall, FS2020.r.se.tall, BB2020.r.se.tall, BC2021.r.se.tall, CNBEIJ.r.se.tall, CNGEZI.r.se.tall,
                     HUNG.r.se.tall, MONG1.r.se.tall, MONG2.r.se.tall, UK1.r.se.tall, UK2.r.se.tall)  
View(temp.rse)
#write.csv(temp.rse, "temp.rse.csv", row.names = FALSE)
temp.rse = read.csv("temp.rse.csv")
unique(temp.rse$isolate_id)
unique(temp.rse$temp)
typeof(temp.rse$temp)

#### making plots of r's ####
#setting temp as numeric to run loess
temp.rse$temp.c = as.character(temp.rse$temp)
typeof(temp.rse$temp.c)
temp.rse$temp.n = as.numeric(temp.rse$temp.c)## this changes temps to 1-5.. 
typeof(temp.rse$temp.n)

unique(temp.rse$temp)
#make a plot with temp on x axis and r on the y axis, each facet isolate id or line color as an isolate ID
rse1=ggplot(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), alpha=0.2)+
  geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Temperature")+
  scale_color_viridis_d(name="Isolate",option="viridis")+
  theme_bw()+
  facet_wrap(~isolate_id)+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        legend.text = element_text(size=10,face="italic"))
rse1


#compare rs at temp for bats
#make a plot with temp on x axis and r on the y axis, each facet isolate id or line color as an isolate ID
temp.rse$temp.n
rse1=ggplot(data=subset(temp.rse),aes(x=isolate_id, y=r,color=isolate_id), alpha=0.2)+
  geom_point()+
  geom_smooth()+
  geom_errorbar(aes(ymax=r+se,ymin=r-se))+
  #geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  #geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Isolate")+
  #scale_color_viridis_d(name="Isolate",option="viridis")+
  theme_bw()+
  facet_wrap(~temp)+
  #facet_wrap(~isolate_id, scales = "free_y")+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12, angle=90),
        panel.grid = element_blank(),
        legend.text = element_text(size=10,face="italic"))
rse1

library(tidyverse)
head(temp.rse)
temp.rse = temp.rse %>%
  mutate(time.period = str_sub(isolate_id, 3, 6))%>%
  mutate(isolate.place = str_sub(isolate_id, 1, 2))

temp.rse$time.period.n = as.numeric(temp.rse$time.period) 

temp.rse2 = temp.rse %>%
  filter(time.period.n>3)%>%
  drop_na(time.period.n) %>%
  mutate(phase = ifelse(time.period.n>2013, "contemporary", "historic"))

temp.rse2$time.period.nf=as.factor(temp.rse2$time.period.n)

##quick plot

rse.growth=ggplot(data=subset(temp.rse2, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18),aes(x=phase, y=r), alpha=0.2)+
  geom_boxplot()+
  geom_jitter(data=subset(temp.rse2, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18), aes(x=phase, y=r,color=time.period.nf), size=2.5,alpha=0.7, width=.05, height=0.0)+
#rse.growth=ggplot(data=temp.rse2,aes(x=phase, y=r), alpha=0.2)+  geom_boxplot()+
  #   geom_jitter(data=temp.rse2, aes(x=phase, y=r,color=time.period.nf), size=3, width=.1, height=.1)+
  #  geom_smooth()+
  #geom_errorbar(aes(ymax=r+se,ymin=r-se))+
  ylim(0,0.6)+
  #geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  #geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Time period")+
  scale_color_viridis_d(name="Year",option="turbo")+
  theme_bw()+
  facet_wrap(~temp, nrow =1)+
  #facet_wrap(~isolate_id, scales = "free_y")+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12, angle=45),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=10,face="italic"))
rse.growth
ggsave(file="/Users/NicholeLaggan/Dropbox/Laggan_WNS_Project/manuscripts/thermal_isolates_expt/hist_cont_comparison.png",width=15,height=5,units="in",dpi=300)


temp.rse2$phase=as.factor(temp.rse2$phase)
temp.rse2$phase= relevel(temp.rse2$phase, ref="historic") 

#renaming isolates 
temp.rse2$isolate.place=as.character(temp.rse2$isolate.place)
unique(temp.rse2$isolate.place)
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="HL"]="Haile's, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="SC"]="Schoharie, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="WH"]="William's Hotel, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="FS"]="Fahnestock, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="BM"]="Barton Mine, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="TY"]="Tawney's, VA"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="PY"]="Peery's, VA"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="NR"]="New River, VA"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="WL"]="William's Lake, NY"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="BB"]="Blackball, IL"
temp.rse2$n_isolate.place[temp.rse2$isolate.place=="BC"]="Bear Creek, WI"
unique(temp.rse2$n_isolate.place)



# same as previous plot but colored by site
head(temp.rse2)
unique(temp.rse2$isolate.place)
rse.growth=ggplot(data=subset(temp.rse2, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18),aes(x=phase, y=r), alpha=0.2)+
  geom_boxplot()+
  geom_jitter(data=subset(temp.rse2, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18), aes(x=phase, y=r,color=n_isolate.place), size=2.5,alpha=0.7, width=.05, height=0.0)+
  #rse.growth=ggplot(data=temp.rse2,aes(x=phase, y=r), alpha=0.2)+  geom_boxplot()+
  #   geom_jitter(data=temp.rse2, aes(x=phase, y=r,color=time.period.nf), size=3, width=.1, height=.1)+
  #  geom_smooth()+
  #geom_errorbar(aes(ymax=r+se,ymin=r-se))+
  ylim(0,0.6)+
  #geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  #geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Time period")+
  scale_color_viridis_d(name="Site",option="turbo")+
  theme_bw()+
  facet_wrap(~temp, nrow =1)+
  #facet_wrap(~isolate_id, scales = "free_y")+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12, angle=45),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=10))
rse.growth
ggsave(file="/Users/NicholeLaggan/Dropbox/Laggan_WNS_Project/manuscripts/thermal_isolates_expt/hist_cont_comparison_colorsite.png",width=15,height=5,units="in",dpi=300)





# same as previous plot but colored by site and only matched sites!!!!!
unique(temp.rse2.matched$isolate.place)
temp.rse2.matched=subset(temp.rse2, isolate.place=="SC"|isolate.place=="FS"|isolate.place=="WH"|isolate.place=="BM")


head(temp.rse2)
unique(temp.rse2$isolate.place)
rse.growth=ggplot(data=subset(temp.rse2.matched, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18),aes(x=phase, y=r), alpha=0.2)+
  geom_boxplot()+
  geom_jitter(data=subset(temp.rse2.matched, temp==1|temp==3|temp==5|temp==7|temp==9|temp==11|temp==18), aes(x=phase, y=r,color=n_isolate.place), size=2.5,alpha=0.7, width=.05, height=0.0)+
  #rse.growth=ggplot(data=temp.rse2,aes(x=phase, y=r), alpha=0.2)+  geom_boxplot()+
  #   geom_jitter(data=temp.rse2, aes(x=phase, y=r,color=time.period.nf), size=3, width=.1, height=.1)+
  #  geom_smooth()+
  #geom_errorbar(aes(ymax=r+se,ymin=r-se))+
  ylim(0,0.6)+
  #geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  #geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Time period")+
  scale_color_viridis_d(name="Site",option="turbo")+
  theme_bw()+
  facet_wrap(~temp, nrow =1)+
  #facet_wrap(~isolate_id, scales = "free_y")+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12, angle=45),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size=10))
rse.growth
ggsave(file="/Users/NicholeLaggan/Dropbox/Laggan_WNS_Project/manuscripts/thermal_isolates_expt/hist_cont_comparison_colorsite.matched.png",width=15,height=5,units="in",dpi=300)








head(temp.rse2)
s1=lm(r~phase+as.factor(temp), data=subset(temp.rse2, r>0))
s1=lm(r~phase+as.factor(temp), data=subset(temp.rse2))
summary(s1)
anova(s1)

temp.rse$temp=factor(temp.rse$temp,
                 levels=c("1","3","5", "7", "9", "11","13","15","18","21"))

# looking at r's x temp, facet wrap by isolate
rse1=ggplot(data=temp.rse,aes(x=temp, y=r,color=isolate_id), alpha=0.2)+
  geom_point()+
  geom_smooth()+
  geom_errorbar(aes(ymax=r+se,ymin=r-se))+
  #geom_point(data=temp.rse,aes(x=temp.n, y=r,color=isolate_id), position=position_jitter(w=0.1), size=2.5, alpha=0.5)+ 
  #geom_smooth(data=temp.rse, aes(x=temp.n, y=r, color=isolate_id), method="loess", size=0.5, se=FALSE, span=1.5)+
  #geom_smooth(data=temp.rse, aes(x=temp, y=r), method="loess", size=0.5, se=se, span=1.5)+
  ylab("Growth rate (r)")+
  xlab("Temperature")+
  #scale_color_viridis_d(name="Isolate",option="viridis")+
  theme_bw()+
  facet_wrap(~isolate_id, scales="free_y")+
  #facet_wrap(~isolate_id)+
  theme(axis.title.y = element_text(size = 16),
        axis.title.x=element_text(size=16),
        axis.text=element_text(size=12),
        panel.grid = element_blank(),
        legend.text = element_text(size=10,face="italic"))
rse1










temp.rse

mod.test = lm(r~isolate_id+0, data=temp.rse)
summary(mod.test)
library(car)
Anova(mod.test)
emmeans(mod.test, ~ isolate_id)
prs = emmeans(mod.test, ~ isolate_id)
pairs(prs)
