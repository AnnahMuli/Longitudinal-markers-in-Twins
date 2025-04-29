############################################################################################
#  Joint Model: Weibull baseline hazard and gamma frailty                                  #
#  h0 = lambda*rho*t^(rh0-1)   and   H0 = lambda*t^(rho)                                   #
############################################################################################
#                                                                                          #
#  The function arguments are:                                                             #
#   - data     : a data.frame containing all the variables;                                #
#   - cluster  : the name of the variable in data containing cluster(twin) IDs;            #
#   - id       : the name of the variable in data containing individual IDs;               #
#   - Marker    : the longitudinal covariate measurements                                  #
#   - MarkerTime  : time when the covariate measurements were taken                        #
#   - status   : censoring indicator (1: non censored; 0 censored)                         #
#   - SurvTime : Survival time                                                             #
#   - A0       : Age at entry                                                              #
#   - inip     : initial values for the parameters (Lambda,Rho,alpha,theta,sigma_u,sigma_e)#
#   - hess     : hessian matrix to generate standard error                                 #
#   - randslope: include id level random slope in the longitudinal model                   #
############################################################################################

# Likelihood for randomslope=F
likelihood.Joint_F <- function(p, LongData, SurvData, Beta0, Beta1) {
  lambda <- exp(p[4])
  rho <- exp(p[6])
  theta <- exp(p[2])  # s=log(var.eps) <=> exp(s)=var.eps
  alpha <- p[3]
  sigma_u <- exp(p[1])
  sigma_e <- exp(p[5])
  
  n <- length(unique(LongData$cluster))
  L.ni <- as.vector(table(LongData$cluster))
  
  yit <- LongData$yij - Beta0 - LongData$bij - Beta1 * LongData$start
  L.yit <- split(yit, LongData$cluster)
  S.status <- aggregate(SurvData$status, by = list(SurvData$cluster), FUN = sum)[2]
  
  hit <- ((lambda * rho * SurvData$ST^(rho - 1)) * exp(alpha * (Beta0 + SurvData$bij + Beta1 * SurvData$ST)))^SurvData$status
  S.hit <- split(hit, SurvData$cluster)
  S.ni <- as.vector(table(SurvData$cluster))
  
  Hit <- unlist(sapply(1:nrow(SurvData), function(i) adaptIntegrate(function(s)
    exp(alpha * (Beta0 + SurvData$bij[i])) * s^(rho - 1) * exp(alpha * Beta1 * s), lowerLimit = c(0), upperLimit = c(SurvData$ST[i]))[1]), use.names = FALSE)
  S.Hit <- split(Hit, SurvData$cluster)
  
  H0it <- unlist(sapply(1:nrow(SurvData), function(i) adaptIntegrate(function(s)
    exp(alpha * (Beta0 + SurvData$bij[i])) * s^(rho - 1) * exp(alpha * Beta1 * s), lowerLimit = c(0), upperLimit = c(SurvData$A0[i]))[1]), use.names = FALSE)
  S0.Hit <- split(H0it, SurvData$cluster)
  
  lik <- sum(log(unlist(sapply(1:n, function(i)
    adaptIntegrate(function(uv)
      as.numeric(((2 * pi * (sigma_e^2))^(-L.ni[i] / 2) * exp(-1 / (2 * (sigma_e^2)) * sum((L.yit[[i]] - uv[1])^2))) *
                   (prod(S.hit[[i]]) * (exp(alpha * uv[1]))^S.status[i,]) *
                   ((1 / theta + sum(S0.Hit[[i]]) * lambda * rho * exp(alpha * uv[1]))^(1 / theta) / gamma(1 / theta)) *
                   (gamma(S.status[i,] + 1 / theta) / ((1 / theta + sum(S.Hit[[i]]) * lambda * rho * exp(alpha * uv[1]))^(S.status[i,] + 1 / theta)))) *
        dnorm(uv[1], mean = 0, sd = sigma_u), lowerLimit = c(-10), upperLimit = c(10))[1]))))
  
  return(-lik)
}

# Likelihood for randomslope=T
likelihood.Joint_T <- function(p, LongData, SurvData, Beta0, Beta1) {
  lambda <- exp(p[4])
  rho <- exp(p[6])
  theta <- exp(p[2])  # s=log(var.eps) <=> exp(s)=var.eps
  alpha <- p[3]
  sigma_u <- exp(p[1])
  sigma_e <- exp(p[5])
  
  n <- length(unique(LongData$cluster))
  L.ni <- as.vector(table(LongData$cluster))
  
  yit <- LongData$yij - Beta0 - LongData$bij0 - LongData$bij1 * LongData$start - Beta1 * LongData$start
  L.yit <- split(yit, LongData$cluster)
  S.status <- aggregate(SurvData$status, by = list(SurvData$cluster), FUN = sum)[2]
  
  hit <- ((lambda * rho * SurvData$ST^(rho - 1)) * exp(alpha * (Beta0 + SurvData$bij0 + SurvData$bij1 * SurvData$ST + Beta1 * SurvData$ST)))^SurvData$status
  S.hit <- split(hit, SurvData$cluster)
  S.ni <- as.vector(table(SurvData$cluster))
  
  Hit <- unlist(sapply(1:nrow(SurvData), function(i) adaptIntegrate(function(s)
    exp(alpha * (Beta0 + SurvData$bij0[i])) * s^(rho - 1) * exp(alpha * (Beta1 + SurvData$bij1[i]) * s), lowerLimit = c(0), upperLimit = c(SurvData$ST[i]))[1]), use.names = FALSE)
  S.Hit <- split(Hit, SurvData$cluster)
  
  H0it <- unlist(sapply(1:nrow(SurvData), function(i) adaptIntegrate(function(s)
    exp(alpha * (Beta0 + SurvData$bij0[i])) * s^(rho - 1) * exp(alpha * (Beta1 + SurvData$bij1[i]) * s), lowerLimit = c(0), upperLimit = c(SurvData$A0[i]))[1]), use.names = FALSE)
  S0.Hit <- split(H0it, SurvData$cluster)
  
   lik <- sum(log(unlist(sapply(1:n, function(i)
    adaptIntegrate(function(uv)
      as.numeric(((2 * pi * (sigma_e^2))^(-L.ni[i] / 2) * exp(-1 / (2 * (sigma_e^2)) * sum((L.yit[[i]] - uv[1])^2))) *
                   (prod(S.hit[[i]]) * (exp(alpha * uv[1]))^S.status[i,]) *
                   ((1 / theta + sum(S0.Hit[[i]]) * lambda * rho * exp(alpha * uv[1]))^(1 / theta) / gamma(1 / theta)) *
                   (gamma(S.status[i,] + 1 / theta) / ((1 / theta + sum(S.Hit[[i]]) * lambda * rho * exp(alpha * uv[1]))^(S.status[i,] + 1 / theta)))) *
        dnorm(uv[1], mean = 0, sd = sigma_u), lowerLimit = c(-10), upperLimit = c(10))[1]))))
  
  return(-lik)
}

# Main function to estimate parameters
fit_JM <- function(data, cluster, id, inip, Marker, MarkerTime, status, SurvTime, A0, hess, randslope) {
  Dat <- data
  Dat$id <- Dat[[deparse(substitute(id))]]
  Dat$start <- Dat[[deparse(substitute(MarkerTime))]]
  Dat$yij <- Dat[[deparse(substitute(Marker))]]
  Dat$status <- Dat[[deparse(substitute(status))]]
  Dat$cluster <- Dat[[deparse(substitute(cluster))]]
  Dat$ST <- Dat[[deparse(substitute(SurvTime))]]
  Dat$A0 <- Dat[[deparse(substitute(A0))]]
  
  if (randslope) {
    # Random intercept and slope model
    model <- lmer(yij ~ start + (1 | id) + (0 + start | id) + (1 | cluster), 
                  data = Dat, REML = FALSE, 
                  control = lmerControl(optimizer = "Nelder_Mead"))
    # Extract random effects: ID
    bij <- ranef(model)$id
    B_ij <- data.frame(id = rownames(bij), bij0 = bij[[1]], bij1 = bij[[2]])
    } 
  else {
    # Random intercept only model
    model <- lmer(yij ~ start + (1 | id) + (1 | cluster), data = Dat)
    # Extract random effects: ID
    bij <- ranef(model)$id
    B_ij <- data.frame(id = rownames(bij), bij = bij[[1]])
  }
  
  # Extract fixed effects
  Beta0 <- coef(summary(model))[1, 1]
  Beta1 <- coef(summary(model))[2, 1]
  
  
  # Merge random effects by ID
  Merged_data <- merge(Dat, B_ij, by = "id")
  
  # Extract and merge random effects: Cluster
  ui <- ranef(model)$cluster
  Ui_cluster <- data.frame(cluster = rownames(ui), ui = ui[[1]])
  Merged_data1 <- merge(Merged_data, Ui_cluster, by = "cluster")
  
  # Ensure bij values are numeric if randslope
  if (randslope) {
    Merged_data1$bij0 <- as.numeric(Merged_data1$bij0)
    Merged_data1$bij1 <- as.numeric(Merged_data1$bij1)
  }
  
  # Filter longitudinal data
  LongData11 <- Merged_data1[Merged_data1$ST > Merged_data1$start, ]
  LongData <- subset(LongData11, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
  
  # Filter survival data
  SurvData1 <- subset(LongData, !duplicated(id))
  SurvData <- subset(SurvData1, ave(id, cluster, FUN = function(x) length(unique(x))) > 1)
  
  
  # Choose the appropriate likelihood function based on randslope
  Lik_JM <- if (randslope) {
    function(p) likelihood.Joint_T(p, LongData, SurvData, Beta0, Beta1)
  } else {
    function(p) likelihood.Joint_F(p, LongData, SurvData, Beta0, Beta1)
  }
  
  # Optimization
  system.time( t2 <- nlminb(inip, Lik_JM, hessian = hess, control = list(rel.tol = 1e-6, step.min = 1e-3, trace = 1))   )
  
  if (hess == TRUE) {
    std3 <- sqrt(diag(solve(hessian(Lik_JM, t2$par))))
    Output2 <- round(matrix(c(t2$par[3], std3[3], exp(t2$par[2]), exp(t2$par[2]) *std3[2],
                              exp(t2$par[4]), exp(t2$par[4]) * std3[4],
                              exp(t2$par[6]), exp(t2$par[6]) * std3[6],
                              exp(t2$par[1]), exp(t2$par[1]) * std3[1],
                              exp(t2$par[5]), exp(t2$par[5]) * std3[5]),
                            nrow = 6, byrow = TRUE), 4)
    colnames(Output2) <- c("Estimate", "SE")
    rownames(Output2) <- c("alpha", "theta", "lambda", "rho", "Sigma_u", "Sigma_e")
    return(Output2)
  } else {
    output <- round( c(sigma_u_JM = exp(t2$par[1]), Theta_JM = exp(t2$par[2]),
                alpha_JM = t2$par[3], Lambda_JM = exp(t2$par[4]),
                sigma_e_JM = exp(t2$par[5]), Rho_JM = exp(t2$par[6])), 4)
    return(output)
  }
}


# Data Application step;
library(dplyr); library(numDeriv);library(lme4) ; library(doBy) ; library(cubature)
DatH<-read.csv( file="C:/Users/staammu/Desktop/SimulData.csv" )
attach(DatH)
fit_JM( data=DatH, cluster = TwinID, id=id, inip=c(  log(0.01), log(0.5), 2, log(0.05), log(0.05), log(2) ) ,
        Marker=Marker,  MarkerTime=MarkerTime, status=status,SurvTime=SurvTime, A0=A0,
        hess=T, randslope=T)
