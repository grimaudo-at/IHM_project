TestTemp0.8<-nls(ColonySize~K/(1+((K/N0)-1)*exp(-(a*(((1+c*exp(-p*Temp))^-1)-exp(-(Tmax-Temp)/(Tmax-Topt))))
                                                *Day)), 
                 data=Data,
                 start=list(r=0.2), trace=TRUE)

#Above is the structure of the logistic pathogen growth model. The nls() function takes data and fist the logistic growth model 
#provided to the data using least squares method.

#Bayesian alternatives to this approach are the packages 'JAGS' and 'BRMS'

Test<-nls(r~a*(((1+c*exp(-p*Temp))^-1)-exp(-(Tmax-Temp)/(Tmax-Topt))), 
          data=JohnsonLewin,
          start=list(a=2, c=5, p=0.3, Topt=14, Tmax=19), trace=TRUE,
          control=nls.control(printEval=TRUE, maxiter=1000, warnOnly=TRUE),
          upper=c(Inf,Inf,Inf,20,20), lower=c(0,0,0,0,0), algorithm = "port")



model <- "
model {
#Priors for fixed effects - 
a ~ dunif(0, 5) ##a shape parameter
g ~ dunif(0, 30) ##a shape parameter
p ~ dunif(0, 1) ##a shape parameter
Tmax ~ dnorm(21, 2) ##Thermal maximum
Topt ~ dnorm(14, 3) ##Delta Temp

sigma ~ dunif(0, 100) # standard deviation
tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS

#Likelihoods
for (i in 1:NCoef) {
y[i] ~ dnorm(mu[i], tau) #tau is precision (1 / variance)
mu[i] <- a*((pow((1+g*exp(-p*Temp[i])), -1))-exp(-(Tmax-Temp[i])/(Tmax - Topt)))
}

}"


Tmax_verant<-as.numeric(jags.m2$BUGSoutput$mean[1]) #21.5 [20.54-22.52]
a_verant<-as.numeric(jags.m2$BUGSoutput$mean[3]) #0.719 [0.419-1.211]
g_verant<-as.numeric(jags.m2$BUGSoutput$mean[5]) #5.262 [1.995-19.459]
p_verant<-as.numeric(jags.m2$BUGSoutput$mean[6]) #0.309 [0.156-0.819]
Topt_verant<-as.numeric(jags.m2$BUGSoutput$mean[2]) #14.05 [95% CI = 12.965-15.162]
## THESE ARE THE PARAMETERS THAT SKYLAR ESTIMATED USING VERANT'S DATA AND JAGS
## THESE PARAMETERS ARE USED TO ESTIMATE r