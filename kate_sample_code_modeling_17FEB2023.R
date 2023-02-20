

##what kate thinks this function will look like...


TestTemp0.8<-nls(ColonySize~K/(1+((K/N0)-1)*exp(-(a*(((1+c*exp(-p*Temp))^-1)-exp(-(Tmax-Temp)/(Tmax-Topt))))
                                                  *Day)), 
                 data=Data,
                 start=list(r=0.2), trace=TRUE)

##skylar - nls fit of logistic model
TestTemp0.8<-nls(ColonySize~K/(1+((K/N0)-1)*exp(-r*Day)), 
                 data=Data,
                 start=list(r=0.2), trace=TRUE)


#skylar - Fit of JohnsonLewin curve
Test<-nls(r~a*(((1+c*exp(-p*Temp))^-1)-exp(-(Tmax-Temp)/(Tmax-Topt))), 
          data=JohnsonLewin,
          start=list(a=2, c=5, p=0.3, Topt=14, Tmax=19), trace=TRUE,
          control=nls.control(printEval=TRUE, maxiter=1000, warnOnly=TRUE),
          upper=c(Inf,Inf,Inf,20,20), lower=c(0,0,0,0,0), algorithm = "port")


##Define model - from JAGS
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

#Bundle data to pass to Jags
NCoef<-length(JohnsonLewin$Coef)
Temp<-JohnsonLewin$Temp
y<-JohnsonLewin$Coef
data<-list(y=y, Temp=Temp, NCoef=NCoef)

##Inits function
inits<-function(){list(a=rnorm(n = 1, mean = 2, sd = 0.5),
                       g=rnorm(n = 1, mean = 10, sd = 1),
                       p=rnorm(n = 1, mean = 0.5, sd = 0.2),
                       Tmax=rnorm(n = 1, mean = 21, sd = 0.5),
                       Topt=rnorm(n = 1, mean = 14, sd = 1))}
##Parameters to estimate
params<-c("a", "g", "p", "Tmax", "sigma", "Topt")

#Run the model and pull posterior samples
jags.m2 <- jags(data = data, inits = inits,
                parameters.to.save = params, n.chains = 3, n.iter = 60000, #niter should be much higher
                n.burnin = 30000, model.file = textConnection(model))
jags.m2
#jags.m<- as.mcmc(jags.m)
summary(jags.m2) 


##ALEX - FIND THE PARAMETERS TO SIMULATE HERE!###
#Save the params to use w/ field data
Tmax_verant<-as.numeric(jags.m2$BUGSoutput$mean[1]) #21.5 [20.54-22.52]
a_verant<-as.numeric(jags.m2$BUGSoutput$mean[3]) #0.719 [0.419-1.211]
g_verant<-as.numeric(jags.m2$BUGSoutput$mean[5]) #5.262 [1.995-19.459]
p_verant<-as.numeric(jags.m2$BUGSoutput$mean[6]) #0.309 [0.156-0.819]
Topt_verant<-as.numeric(jags.m2$BUGSoutput$mean[2]) #14.05 [95% CI = 12.965-15.162]

#Predictions
JohnsonLewin$Test2<-a_verant*(((1+g_verant*exp(-p_verant*JohnsonLewin$Temp))^-1)-exp(-(Tmax_verant-JohnsonLewin$Temp)/(Tmax_verant-Topt_verant)))



