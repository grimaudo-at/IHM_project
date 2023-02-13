library(tidyverse)
library(ggplot2)
library(lme4)
library(effects)

dat<-read.csv("/Users/alexg8/Dropbox/Grimaudo_WNS_Project/Data/IHM Project/ibutton_IDs_weights.csv")

dat$c<-1

agg<-dat %>%
  group_by(site) %>%
  mutate(tot.t=sum(c)) %>%
  mutate(tot.rec=sum(c[recovered=="Y"])) %>%
  mutate(tot.fun=sum(c[functional=="Y"])) %>%
  mutate(tot.22L.dep=sum(c[model=="DS1922L"])) %>%
  mutate(tot.21G.dep=sum(c[model=="DS1921G"])) %>%
  mutate(tot.22L.rec=sum(c[recovered=="Y" & model=="DS1922L"])) %>%
  mutate(tot.21G.rec=sum(c[recovered=="Y" & model=="DS1921G"])) %>%
  mutate(tot.22L.fun=sum(c[functional=="Y" & model=="DS1922L"])) %>%
  mutate(tot.21G.fun=sum(c[functional=="Y" & model=="DS1921G"]))
agg <- as.data.frame(agg)
agg$site <- as.factor(agg$site)

summ <- aggregate(tot.fun ~ site, data=agg, FUN=mean)
sum(summ$tot.fun)

### is weight predictive of retrieval? 

dat$recovered[dat$recovered=="Y"]<-1
dat$recovered[dat$recovered=="N"]<-0
dat$recovered<-as.integer(dat$recovered)

m1<-glm(recovered~weight_g, family = "binomial", data=dat);summary(m1)
plot(allEffects(m1))
#perhaps some effect of transmitter weight on recovery probability 

p1<-ggplot(aes(x=weight_g, y=recovered), data=dat)+
  geom_point();p1

