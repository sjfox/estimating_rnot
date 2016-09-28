## Fitting functions


nllikelihood <- function(init, parms = seir_params(), obsDat) {
  # simDat <- simEpidemic(init, parms=parms, deg_dist, tseq=obsDat$day, method=method)
  #nlls <- -dnbinom(x = round(obsDat$X..WEIGHTED.ILI*10000), mu = simDat$incI*10000, size=1000, log=T)
  
  simDat <- simEpidemic(init, parms=parms, deg_dist, tseq=seq(0,obsDat[nrow(obsDat), "day"]))
  
  maxDat <- obsDat$day[match(max(obsDat$X..WEIGHTED.ILI), obsDat$X..WEIGHTED.ILI)]
  maxSim <- simDat$time[match(max(simDat$incI), simDat$incI)]
  
  nlls <- -dnorm(x = maxDat, mean = maxSim, sd = 4, log=T)
  
  return(sum(nlls))
}

sumsq <- function(init, parms = seir_params(), obsDat) {
  tseq <- seq(0, max(obsDat$day))
  simDat <- simEpidemic(init, parms=parms, tseq=tseq)
  matchedTimes <- simDat$time %in% obsDat$day
  return(sum((simDat$incI[matchedTimes] - obsDat$incidence)^2))
}

sumsqflex <- function(init, parms = seir_params(), obsDat) {
  # Sum squares, but allows for matching the time series flexibly in terms of 
  # day number
  # Will allow model to have initial epidemic begin anywhere up to a full epidemic before or after
  #browser()
  maxDay <- max(obsDat$day)
  tseq <- seq(0, maxDay*2+1)
  simDat <- simEpidemic(init, parms=parms, tseq=tseq)
  beforeEpi <- data.frame(time = seq(-maxDay,-1), incI = 0)
  
  simDat <- rbind(beforeEpi, simDat[,c("time", "incI")])
  
  sumsquares <- vector('numeric', length = 2*maxDay)
  minShift <- 0
  minVal <- 1000
  for(ii in -maxDay:maxDay){
    simDat$flextime <- simDat$time - ii
    matchedTimes <- simDat$flextime %in% obsDat$day
    sumsquares[ii+maxDay+1] <- sum((simDat$incI[matchedTimes] - obsDat$incidence)^2)
    # sumsquares[ii+maxDay+1] <- sum(abs(simDat$incI[matchedTimes] - obsDat$incidence)/obsDat$incidence) ## didn't work well
    # sumsquares[ii+maxDay+1] <- - sum(dpois(x = round(10000*simDat$incI[matchedTimes]) , lambda = round(10000*obsDat$incidence), log=T))
    if(sumsquares[ii+maxDay+1] < minVal) {
      minShift <- ii
      minVal <- sumsquares[ii+maxDay+1]
    }
  }
  
  return(list(minSS = min(sumsquares, na.rm=T), shift = minShift))
}


objFXN <- function(est_parms, obsDat, parms=seir_params(), mod) {
  parms <- subs_parms(est_parms, parms)
  # init <- with(as.list(parms), c(theta=1, thetaDot= -beta*phiI0, I=phiI0, R=phiR0, cumI=phiI0))
  init <- with(as.list(parms), c(theta=1, phiI=phiI0, I=phiI0, R=phiR0, cumI=phiI0))
  
  if(mod=="ll"){
    nllikelihood(parms = parms, obsDat = obsDat, init=init)  
  } else if (mod=="ss"){
    sumsq(parms = parms, obsDat = obsDat, init=init)  
  }else if (mod=="ssf"){
    sumsqflex(parms = parms, obsDat = obsDat, init=init)$minSS  
  }
}


subs_parms <- function(sub_parms=NULL, ## parameters to change in ref_parms
                       ref_parms) {  ## Reference params
  within(ref_parms, { ## subs fitting parameters into reference parameter vector
    
    for(nm in names(sub_parms)) {
      assign(nm, sub_parms[[nm]])
    }
    rm(nm)
  })
}
