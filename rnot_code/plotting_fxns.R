## Functions for plotting results

plot_fit <- function(fit_parms, parms, data){
  testParms <- subs_parms(sub_parms = fit_parms, ref_parms = parms)
  testInit <- with(as.list(testParms), c(theta=1, phiI= phiI0, I=phiI0, R=phiR0, cumI=phiI0))
  shift <- sumsqflex(testInit, testParms, obsFlu)$shift
  simData <- simEpidemic(testInit, 1:250, testParms)
  simData$time <- simData$time - shift
  plot(obsFlu$day, obsFlu$incidence)  
  lines(simData$time, cumToInc(simData$cumI))
}