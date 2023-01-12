#' @title Fit cp carcass persistence models
#'
#' @description Carcass persistence is modeled as survival function
#'
#' @details uses \code{survival} package to fit.
#'
#' @param cpdata matrix of interval censored carcass persistence times
#'
#' @param dist Name of the persistence distribution family: Weibull,
#'  lognormal, loglogistic, or exponential
#'
#' @return answer
#'
#' @export
#'
est_cp0 <- function(cpdata, dist){
  distn <- tolower(dist)
  if (distn == "log-logistic") distn <- "loglogistic"
  pdnms <- c("exponential", "weibull", "loglogistic", "lognormal")
  if (!distn %in% pdnms){
    msg <- "error. persistence distribution must be one of the following: "
    for (nm in pdnms) msg<-paste(msg, nm)
    stop(msg)
  }
  if (ncol(cpdata) < 2){
    stop(paste0("est_cp0 requires two columns with data for CPmin and CPmax."))
  } else {
    cpdata<-cpdata[, 1:2]
    if (!is.data.frame(cpdata)) cpdata<-data.frame(cpdata)
    names(cpdata)<-c("CPmin", "CPmax")
  }
  if (!is.numeric(cpdata[,1]) || !is.numeric(cpdata[,2]) ||
      sum(cpdata[,1] < 0) > 0 || sum(cpdata[,1] > cpdata[,2]) > 0) {
    stop("Error in cpdata. Cannot calculate.")
  }
  xind<-which(cpdata$CPmin == 0 & cpdata$CPmax == Inf)
  if (length(xind)>0){
    cpdata$CPmin<-cpdata$CPmin[-xind]
    cpdata$CPmax<-cpdata$CPmax[-xind]
  }
  cpdata$CPmin <- pmax(cpdata$CPmin, 0.001)
  event<-ifelse(cpdata$CPmin == cpdata$CPmax, 1, ifelse(cpdata$CPmax == Inf, 0, 3))
  left<-cpdata$CPmin
  right<-cpdata$CPmax
  right[event==0]<-cpdata$CPmin[event==0]

  if (sum(cpdata$CPmax==Inf) == length(cpdata$CPmax)){
    stop("No carcasses removed in persistence trials. Cannot fit model using ",
      "'survival' package. Use CP analysis submodule of Single Class module to ",
      "fit model.")
  }
  surv<- survival::Surv(time=left, time2=right, event=event, type=c('interval'))
  # fit survival models to persistence data and plot
  return(survival::survreg(surv ~ 1, dist = distn))
}
#' @title generate random cp parameters or persistence times
#'
#' @description Given a fitted cp model (\code{survreg} object), generate random
#'  \code{pda} and \code{pdb} parameters or random persistence times.

#'  NOTE: This function is likely to move to call GenEst's \code{rcp} in future.
#'  This will not change the results, but the GenEst version is more nicely
#'  coded and keeping some coherence among the models is helpful.
#'
#' @param cpmodel Fitted cp model ((\code{survreg} object))
#'
#' @param nsim Number of simulation draws
#'
#' @param option \code{option = "parms"} returns random draws of parameters from
#'  the fitted model; \code{option != "parms"} returns random draws of carcass
#'  persistence times
#'
#' @return answer
#'
#' @export
#'
sim_cp0 <- function(cpmodel, nsim, option = "parms"){
  CPab <- array(dim=c(nsim,2))
  if (cpmodel$dist == "exponential"){
    CPab[, 2] <- exp(rnorm(nsim,
      mean = cpmodel$coefficients, sd = sqrt(cpmodel$var[1])))
    CPab[, 1] <- 1/CPab[, 2]
  } else if (cpmodel$dist == "weibull"){
    CPparms <- MASS::mvrnorm(nsim, mu = cpmodel$icoef, Sigma = cpmodel$var)
    CPab[,1] <- 1/exp(CPparms[, 2]) #shape
    CPab[,2] <- exp(CPparms[, 1]) # scale
  } else if (cpmodel$dist == "loglogistic"){
    # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])
    CPparms <- MASS::mvrnorm(nsim, mu = cpmodel$icoef, Sigma = cpmodel$var)
    CPab[,1] <- 1/exp(CPparms[, 2]) #shape
    CPab[,2] <- exp(CPparms[, 1]) # scale
  } else if (cpmodel$dist == "lognormal"){
    # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
    CPparms <- MASS::mvrnorm(nsim, mu = cpmodel$icoef, Sigma = cpmodel$var)
    CPab[, 1] <- exp(CPparms[, 2])^2 #shape
    CPab[, 2] <- CPparms[, 1] # scale
  }
  if (option == "parms") {
    colnames(CPab) <- c("pda", "pdb")
    return(CPab)
  } else {
    return(switch(cpmodel$dist,
      "exponential" = rexp(nsim, CPab[,1]),
      "weibull"     = rweibull(nsim, shape = CPab[,1], scale = CPab[,2]),
      "loglogistic" = actuar::rllogis(nsim, shape = CPab[,1], scale = CPab[,2]),
      "lognormal"   = rlnorm(nsim, sdlog = sqrt(CPab[,1]), meanlog = CPab[,2])
    ))
  }
}

#' @title retrieve EoA parameterization from \code{survival} parameterization of
#'  a fitted cp model (or \code{survreg} object with exponential, weibull,
#'  lognormal, or loglogistic distribution)
#'
#' @param cpmodel fitted cp model (from \code{est_cp0} or \code{survreg} object)
#'
#' @return 2-vector of pda and pdb
#'
#' @export
#'
getab <- function(cpmodel){
  if (!"survreg" %in% class(cpmodel)) stop("getab requires survreg object")
  if (!cpmodel$dist %in% c("exponential", "weibull", "loglogistic", "lognormal"))
    stop("persistence distribution must be ",
         "exponential, weibull, loglogistic, or lognormal")
  if (cpmodel$dist == "exponential"){
      pdb0 <- exp(cpmodel$coefficients[1])
      pda0 <- 1/pdb0
  } else if (cpmodel$dist == "weibull"){
      # NOTE: for weibull...
      #  a = shape = 1/cpmodel$scale
      #  b = scale = exp(cpmodel$coef[1])
    pdb0 <- exp(cpmodel$coefficients[1])
    pda0 <- 1/cpmodel$scale
  } else if (cpmodel$dist == "loglogistic"){
      # NOTE: for loglogistic...
      #  a = shape = 1/cpmodel$scale
      #  b = scale = exp(cpmodel$coef[1])
      pdb0 <- exp(cpmodel$coefficients[1])
      pda0 <- 1/cpmodel$scale
  } else if (cpmodel$dist == "lognormal"){
      # NOTE: lognormal...
      #  a = sdlog^2 = cpmodel$scale^2
      #  b = meanlog = cpmodel$coef[1]
      pdb0 <- cpmodel$coefficients[1]
      pda0 <- cpmodel$scale^2
  }
  ans <- c(pda0, pdb0)
  names(ans) <- c("pda", "pdb")
  return(ans)
}
