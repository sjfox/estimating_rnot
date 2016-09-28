rm(list=ls())

library(deSolve)
sapply(c("rnot_code/network_fxns.R", "rnot_code/fitting_fxns.R", "rnot_code/plotting_fxns.R"), source)


## Set up the parameters for the run




##########################
# Setup data for fitting
# Choose data from 2008-2009 flu season prior to pandemic
# only the increasing areas
##########################
dataOriginal <- read.csv("data/ILINet.csv", skip=1, na.strings = "X")
whoOriginal <- read.csv("data/WHO_NREVSS.csv", na.strings = "X")

## Isolate data for 2009 flu season prior to pandemic
## Week 50 chosen of 2008, because it's roughly equal to the end of the season before the pandemic
# ili <- dataOriginal[dataOriginal$YEAR==2009, ]
# who <- whoOriginal[whoOriginal$YEAR==2009, ]
plot(dataOriginal$X..WEIGHTED.ILI*whoOriginal$PERCENT.POSITIVE, type="l")
##############################
# Use these datasets for fitting the seasonal 2009 epidemic
##############################
ili <- dataOriginal[(dataOriginal$YEAR==2009 & dataOriginal$WEEK <16), ]
who <- whoOriginal[(whoOriginal$YEAR==2009 & whoOriginal$WEEK<16), ]


obsFlu <- data.frame( incidence = ili$X..WEIGHTED.ILI/100 * who$PERCENT.POSITIVE/100)
obsFlu$day <- (seq_along(obsFlu$incidence)-1)*7
plot(obsFlu$day, obsFlu$incidence, type="l")

expMean <- function(rate, desiredMean, degs=seq(1,10000)){
  sum(seq_along(degs) * normalize(dexp(rate = rate, x = degs))) - desiredMean
}

meanDeg <- 16
degs <- seq(1,10000)
expRate <- uniroot(expMean, interval = c(1/1000,1), desiredMean=meanDeg, degs=degs)$root
expRate
degExp <- normalize(dexp(rate = expRate, x = degs))

degUnif <- rep(0,10000)
degUnif[meanDeg] <- 1

# source("network_fxns.R")

################################
## Fit Uniform network
################################
unif_parms <- seir_params(deg_dist=degUnif)
est_parms <- c(beta =parms$beta, 
               eta=parms$eta, 
               gamma=parms$gamma)
optim_unif <- optim(par = c(beta=0.05, eta=1/2.62, gamma=1/3.38), fn = objFXN, parms=unif_parms, obsDat = obsFlu, mod="ssf")



plot_fit(optim_unif$par, parms=unif_parms, data = obsFlu)

###############################
## Fit Exponential Network
###############################
exp_parms <- seir_params(deg_dist=degExp)

est_parms <- c(beta =parms$beta, 
               eta=parms$eta, 
               gamma=parms$gamma)
optim_exp <- optim(par = c(beta=0.05, eta=1/2.62, gamma=1/3.38), fn = objFXN, parms=exp_parms, obsDat = obsFlu, mod="ssf")



plot_fit(optim_exp$par, parms=exp_parms, data = obsFlu)


calc_fit_rnot <- function(fit_parms, ref_parms){
  parms <- subs_parms(fit_parms, ref_parms)
  calcTheoRNot(beta=parms$beta, parms$deg_dist, gamma = parms$gamma)
}
  
calc_fit_rnot(optim_unif$par, unif_parms)
calc_fit_rnot(optim_exp$par, exp_parms)
  
  
