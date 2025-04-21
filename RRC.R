
#################################################################################
#  RiskSet regression calibration: Weibull baseline hazard and gamma frailty    #
#  h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)                        #
#################################################################################
#                                                                               #
#  The function arguments are:                                                  #
#   - data     : a data.frame containing all the variables;                     #
#   - cluster  : the name of the variable in data containing cluster(twin) IDs; #
#   - id       : the name of the variable in data containing individual IDs;    #
#   - Marker    : the longitudinal covariate measurements                       #
#   - MarkerTime  : time when the covariate measurements were taken             #
#   - status   : censoring indicator (1: non censored; 0 censored)              #
#   - SurvTime : Survival time                                                  #
#   - A0       : Age at entry                                                   #
#   - inip     : initial values for the parameters (Lambda,Rho,alpha,theta)     #
#   - hess     : hessian matrix to generate standard error                      #
#   - randslope: include id level random slope in the longitudinal model        #
#################################################################################

# ---- Prediction function for randslope = FALSE
predicted_Marker_fun_F <- function(st, Dat) {
  ds_subset <- Dat[Dat$ST >= st & Dat$A0 <= st & Dat$start <= st, ]
  fit <- lmer(yij ~ start + (1|id) + (1|cluster), data = ds_subset)
  newSubset <- data.frame(id = unique(ds_subset$id),
                          cluster = subset(ds_subset, !duplicated(id))$cluster,
                          A0 = subset(ds_subset, !duplicated(id))$A0,
                          ST = subset(ds_subset, !duplicated(id))$ST,
                          status = subset(ds_subset, !duplicated(id))$status,
                          start = st)
  newSubset$predicted_Marker <- predict(fit, newdata = newSubset)
  return(newSubset)
}

# ---- Prediction function for randslope = TRUE
predicted_Marker_fun_T <- function(st, Dat) {
  ds_subset <- Dat[Dat$ST >= st & Dat$A0 <= st & Dat$start <= st, ]
  fit <- lmer(yij ~ start + (1|id) + (0+start|id) + (1|cluster), data = ds_subset, REML = FALSE, 
              control = lmerControl(optimizer = "Nelder_Mead"))
  newSubset <- data.frame(id = unique(ds_subset$id),
                          cluster = subset(ds_subset, !duplicated(id))$cluster,
                          A0 = subset(ds_subset, !duplicated(id))$A0,
                          ST = subset(ds_subset, !duplicated(id))$ST,
                          status = subset(ds_subset, !duplicated(id))$status,
                          start = st)
  newSubset$predicted_Marker <- predict(fit, newdata = newSubset)
  return(newSubset)
}

# ---- Likelihood function for randslope = FALSE
Lik_RRC_fun_F <- function(p, dat1, predicted_Marker, di, e, GG) {
  cumhazStart <- aggregate(out ~ cluster, transform(dat1, out = exp(p[1]) * (start^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  cumhazEntry <- aggregate(out ~ cluster, transform(subset(dat1, !duplicated(id)), out = exp(p[1]) * (start^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  cumhazStop <- aggregate(out ~ cluster, transform(dat1, out = exp(p[1]) * (stop^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  lnhaz <- aggregate(out ~ cluster, transform(dat1, out = status * (predicted_Marker * p[4] + log(exp(p[1]) * exp(p[2]) * stop^(exp(p[2]) - 1)))), sum)[,2]
  
  lik <- e * log(p[3]) -
    GG * log(gamma(1 / p[3])) +
    sum(log(gamma(di + 1 / p[3]))) +
    sum(lnhaz) +
    sum((1 / p[3]) * log(1 + cumhazEntry * p[3])) -
    sum((di + 1 / p[3]) * log(1 + cumhazEntry * p[3] + (cumhazStop - cumhazStart) * p[3]))
  
  return(-lik)
}

# ---- Likelihood function for randslope = TRUE
Lik_RRC_fun_T <- function(p, dat1, predicted_Marker, di, e, GG) {
  cumhazStart <- aggregate(out ~ cluster, transform(dat1, out = exp(p[1]) * (start^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  cumhazEntry <- aggregate(out ~ cluster, transform(subset(dat1, !duplicated(id)), out = exp(p[1]) * (start^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  cumhazStop <- aggregate(out ~ cluster, transform(dat1, out = exp(p[1]) * (stop^exp(p[2])) * exp(predicted_Marker * p[4])), sum)[,2]
  lnhaz <- aggregate(out ~ cluster, transform(dat1, out = status * (predicted_Marker * p[4] + log(exp(p[1]) * exp(p[2]) * stop^(exp(p[2]) - 1)))), sum)[,2]
  
  lik <- e * log(p[3]) -
    GG * log(gamma(1 / p[3])) +
    sum(log(gamma(di + 1 / p[3]))) +
    sum(lnhaz) +
    sum((1 / p[3]) * log(1 + cumhazEntry * p[3])) -
    sum((di + 1 / p[3]) * log(1 + cumhazEntry * p[3] + (cumhazStop - cumhazStart) * p[3]))
  
  return(-lik)
}

# ---- Main function ----
fit_RRC <- function(data, cluster, id, inip, Marker, MarkerTime, status, SurvTime, A0, hess, randslope) {
  
  Dat <- data
  Dat$id <- Dat[[deparse(substitute(id))]] 
  Dat$start <- Dat[[deparse(substitute(MarkerTime))]] 
  Dat$yij <- Dat[[deparse(substitute(Marker))]] 
  Dat$status <- Dat[[deparse(substitute(status))]] 
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ST <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  Dat <- Dat[Dat$A0 < Dat$ST, ]
  Dat <- Dat[order(Dat$ST), ]
  tt <- sort(unique(c(Dat$A0, Dat$ST[Dat$status == 1])))
  nt <- length(tt)
  
  # Choose correct predicted_Marker function based on randslope
  predicted_Marker <- if (randslope) {
    function(st) predicted_Marker_fun_T(st, Dat)
  } else {
    function(st) predicted_Marker_fun_F(st, Dat)
  }
  
  # Run prediction for all unique times
  LMM_Output <- purrr::map_dfr(tt[1:(nt - 1)], predicted_Marker)
  
  # Merge and process output
  Dat2 <- LMM_Output[LMM_Output$ST > LMM_Output$start, ]
  
  # Arrange data in start-stop format
  Dat3 <- Dat2 %>%
    group_by(id) %>%
    mutate(stop = lead(start, default = 1e6),
           stop = ifelse(stop > ST, ST, stop), .after = start) %>%
    filter(stop <= ST) %>%
    ungroup()
  
  Dat3 <- orderBy(~id, data = Dat3)
  Dat3$status[Dat3$stop < Dat3$ST] <- 0
  dat11 <- Dat3[Dat3$stop > Dat3$start, ]
  dat1 <- subset(dat11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
  
  # Compute the necessary values for likelihood computation
  Di <- aggregate(dat1$status, by = list(dat1$cluster), FUN = sum)
  e <- colSums(Di[2])
  di <- Di[2]
  GG <- length(unique(dat1$cluster))
  
  # Choose correct likelihood function
  Lik_RRC <- if (randslope) {
    function(p) Lik_RRC_fun_T(p, dat1, dat1$predicted_Marker, di, e, GG)
  } else {
    function(p) Lik_RRC_fun_F(p, dat1, dat1$predicted_Marker, di, e, GG)
  }
  
  # Optimization
  system.time( t2 <- nlminb(inip, Lik_RRC, hessian = hess, control = list(rel.tol = 1e-6, trace=1))   )
  
  # If Hessian is required, return standard errors
  if (hess) {
    std3 <- sqrt(diag(solve(hessian(Lik_RRC, t2$par))))
    Output2 <- round(matrix(c(t2$par[4], std3[4], t2$par[3], std3[3], exp(t2$par[1]), std3[1] * exp(t2$par[1]), exp(t2$par[2]), std3[2] * exp(t2$par[2])),
                            nrow = 4, byrow = TRUE), 4)
    colnames(Output2) <- c("Estimate", "SE")
    rownames(Output2) <- c("alpha", "theta", "Lambda", "Rho")
    return(Output2)
  } else {
    # Output the results
    output <- round( c(alpha = t2$par[4], theta = t2$par[3],
                       lambda = exp(t2$par[1]), rho = exp(t2$par[2])) ,4)
    return(output)
  }
}


# Data Application step;
library(dplyr); library(numDeriv); library(lme4); library(doBy); library(purrr)

DatH<-read.csv( file="C:/Users/staammu/Desktop/SimulData.csv" )
attach(DatH)
fit_RRC ( data=DatH, cluster = TwinID, id=id, inip=c(  log(0.01), log(2),0.3, 2 ) ,
          Marker=Marker,  MarkerTime=MarkerTime, status=status,SurvTime=SurvTime, A0=A0,
          hess=T, randslope=T)
