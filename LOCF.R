#################################################################################
#  Last observation carried forward: Weibull baseline hazard and gamma frailty  #
#  h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)                        #
#################################################################################
#                                                                               #
#  The function arguments are:                                                  #
#   - data     : a data.frame containing all the variables;                     #
#   - cluster  : the name of the variable in data containing cluster(twin) IDs; #
#   - id       : the name of the variable in data containing individual IDs;    #
#   - MarkerTime  : time when the covariate measurements were taken             #
#   - status   : censoring indicator(1: non censored; 0 censored)               #
#   - SurvTime : Survival time                                                  #
#   - A0       : Age at entry                                                   #
#   - inip     : initial values for the parameters(Lambda,Rho,alpha,theta)      #
#################################################################################
#                                                                               #
#   Last modification : April  07, 2025                                         #
#                                                                               #
#################################################################################

#The log of marginal likelihood function;
Lik_LOCF <- function(p, dat1, e, di, GG) {
  # e -- Total number of non-censored observations
  # GG -- total number of clusters (twin pairs)
  cumhazStart <- as.numeric(aggregate(out ~ cluster, transform(dat1,
                                                               out = exp(p[1]) * (start^exp(p[2])) * exp(zij * p[4])), FUN = sum)[,2])
  
  cumhazEntry <- as.numeric(aggregate(out ~ cluster,
                                      transform(subset(dat1, !duplicated(id)),
                                                out = exp(p[1]) * (start^exp(p[2])) * exp(zij * p[4])), FUN = sum)[,2])
  
  cumhazStop <- as.numeric(aggregate(out ~ cluster, transform(dat1,
                                                              out = exp(p[1]) * (stop^exp(p[2])) * exp(zij * p[4])), FUN = sum)[,2])
  
  lnhaz <- as.numeric(aggregate(out ~ cluster, transform(dat1,
                                                         out = status * (zij * p[4] + log(exp(p[1]) * exp(p[2]) * stop^(exp(p[2]) - 1)))), FUN = sum)[,2])
  
  ## Likelihood Equation
  lik <- e * log(p[3]) -
    GG * log(gamma(1 / p[3])) +
    sum(log(gamma(di + 1 / p[3]))) +
    sum(lnhaz) +
    sum((1 / p[3]) * log(1 + cumhazEntry * p[3])) -
    sum((di + 1 / p[3]) * log(1 + cumhazEntry * p[3] + (cumhazStop - cumhazStart) * p[3]))
  
  return(-lik)
}

fit_LOCF <- function(data, cluster, id, inip, CovTime, status, covariate, SurvTime, A0, hess) {
  Dat <- data
  Dat$id <- Dat[[deparse(substitute(id))]] 
  Dat$start <- Dat[[deparse(substitute(CovTime))]] 
  Dat$zij <- Dat[[deparse(substitute(covariate))]] 
  Dat$status <- Dat[[deparse(substitute(status))]] 
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ATF <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  Dat <- Dat[Dat$A0 < Dat$ATF,]
  
  Dat3 <- Dat %>% group_by(id) %>%
    mutate(stop = lead(start, default = 1e6),
           stop = ifelse(stop > ATF, ATF, stop), .after = start) %>%
    filter(stop <= ATF) %>%
    ungroup()
  
  Dat4 <- Dat3[Dat3$A0 < Dat3$stop,]
  Dat4$start <- pmax(Dat4$A0, Dat4$start)
  Dat4$status[Dat4$stop < Dat4$ATF] <- 0
  dat11 <- Dat4[Dat4$stop > Dat4$start,]
  
  dat1 <- subset(dat11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
  
  Di <- aggregate(dat1$status, by = list(dat1$cluster), FUN = sum)
  e <- colSums(Di[2])
  di <- Di[2]
  GG <- length(unique(dat1$cluster))
  
  t2 <- nlminb(inip, Lik_LOCF, dat1 = dat1, e = e, di = di, GG = GG, hessian = hess, control = list(rel.tol = 1e-6))
  
  output <- c(alpha = t2$par[4], theta = t2$par[3],
              lambda = exp(t2$par[1]),
              rho = exp(t2$par[2]))
  
  if (hess == TRUE) {
    std3 <- sqrt(diag(solve(hessian(Lik_LOCF, t2$par, dat1 = dat1, e = e, di = di, GG = GG)))
    )
    
    Output2 <- round(matrix(c(t2$par[4], std3[4],
                              t2$par[3], std3[3],
                              exp(t2$par[1]), std3[1] * exp(t2$par[1]),
                              exp(t2$par[2]), std3[2] * exp(t2$par[2])), 
                            nrow = 4, byrow = TRUE), 3)
    
    colnames(Output2) <- c("Estimate", "SE")
    rownames(Output2) <- c("alpha", "theta", "Lambda", "Rho")
    return(Output2)
  } else {
    return(output)
  }
}

# Data Application step;
library(numDeriv); library(dplyr);
DatH<-read.csv( file="C:/Users/staammu/Desktop/SimulData.csv" )
attach(DatH)
fit_LOCF ( data=DatH, cluster = TwinID, id=id, inip=c(  log(0.01), log(2),0.3,2 ), MarkerTime=MarkerTime, status=status,SurvTime=SurvTime, A0=A0, Marker=Marker, hess=T)

