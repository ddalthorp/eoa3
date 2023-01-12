#' @title estimate g from fitted pk and cp models and search schedule
#'
#' @description Given a fitted pk model (from \code{\link{est_pk0}}), a fitted
#'  cp model (from \code{\link{est_cp0}} or a \code{survreg} object), and a
#'  search schedule, estimate detection probability.
#'
#' @param pkmodel fitted pk model (from \code{est_pk0})
#' @param cpmodel fitted cp model (from \code{est_cp0} or \code{survreg} object)
#' @param days vector of days since searches begin (\code{days[1]} == 0)
#' @param a fraction of carcasses arriving in the area searched
#' @param v fraction of carcasses arriving in the period spanned by the monitoring
#' @param ... additional arguments (ignored)
#' @examples
#'  pkmodel <- est_pk0(pkdata = pkdata0)
#'  cpmodel <- est_cp0(cpdata = cpdata0, dist = "weibull")
#'  ghat <- est_g0(pkmodel = pkmodel, cpmodel = cpmodel, days = days0, a = 0.4, v = 0.75)
#'  summary(ghat)
#' @return list of parameters for a beta distributions fit to the vectors of
#'  for \eqn{\hat{g}} for the searched area within the period monitored
#'  (\code{$BabRaw}), for the whole site within the period monitored
#'  (\code{$Bab}), and for the whole site extrapolated to the whole year
#'  (\code{$BabAnn}). In addition, the models and parameters that went into the
#'  estimate are included as well (\code{pkmodel, cpmodel, a, v}).
#'
#' @export
#'
est_g0 <- function(pkmodel, cpmodel, days, a = NULL, v = NULL, ...){
  samtype <- ifelse(length(unique(diff(days))) == 1, "Formula", "Custom")
  if (samtype=='Formula'){
    Isam <- days[2] - days[1]
    nsearch <- length(days) - 1
  } else {
    Isam <- round(max(days)/(length(days) - 1), 1)
    nsearch <- length(days) - 1
  }
  if (pkmodel$type == "pOnly") stop("cannot estimate g using a k-free pkmodel")
  ###1. setting estimation control parameters
  ##  a. search limit: required number of searches after arrival to include
  ##   in estimate of searcher efficiency [when the number of searches is high,
  ##   including them all in the estimation is calculation intensive but does
  ##   not contribute signficantly to the result]
  ### Better would be to use SYcalcg with fitted pk and cp models. This would
  ### give greater stability to the calculation of control parameters as well
  ### as more accuracy (and speed) in determining the number of simulation
  ### reps required to get the desired level of precision (for eoa v3)
  sumstat <- summary(pkmodel)
  f0 <- sumstat["p", "mean"]
  k0 <- sumstat["k", "mean"]
  ind1 <- rep(1:nsearch, times = nsearch:1)
  ind2 <- ind1 + 1
  ind3 <- unlist(lapply(1:nsearch, function(x) x:nsearch)) + 1
  schedule.index <- cbind(ind1, ind2, ind3)
  schedule <- cbind(days[ind1], days[ind2], days[ind3])
  nmiss <- schedule.index[,3] - schedule.index[,2]
  maxmiss <- max(nmiss)
  powk <- cumprod(c(1, rep(k0, maxmiss))) # vector of k^i's
  notfind <- cumprod(1 - f0 * powk[-length(powk)])
  nvec <- c(1, notfind) * f0
  # conditional probability of finding a carcass on the ith search (row)
  #  after arrival for given (simulated) searcher efficiency (column)...
  pfind.si <- nvec * powk
  # persistences
  intxsearch <- unique(cbind(
    schedule[,2] - schedule[,1], schedule[,3] - schedule[,2]
  ), MARGIN = 1)
  ab0 <- getab(cpmodel)
  ppersu <- GenEst::ppersist(pda = ab0[1], pdb = ab0[2], dist = cpmodel$dist,
    t_arrive0 = 0,
    t_arrive1 = intxsearch[, 1],
    t_search = intxsearch[, 1] + intxsearch[, 2])
  if (any(is.na(ppersu))) ppersu[is.na(ppersu)] <- max(ppersu[!is.na(ppersu)])
  # fraction of total carcasses that arrive in each search interval
  # ...inference is limited to the monitoring period
  # ...assume uniform arrivals
  arrvec <- (schedule[, 2] - schedule[, 1])/max(days) #
  prob_obs <- numeric(dim(schedule)[1])
  for (i in 1:length(prob_obs)){
    prob_obs[i] <- pfind.si[nmiss[i] + 1] *
      ppersu[which(
        abs(intxsearch[, 1] - (schedule[i, 2] - schedule[i, 1])) < 0.001 &
        abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 2])) < 0.001), ] *
      arrvec[i]
  }
  prob_obs <- prob_obs * ifelse(is.null(a), 1, a)
  ggnm <- numeric(maxmiss + 1)
  for (missi in 0:maxmiss){
    ggnm[missi + 1] <- sum(prob_obs[nmiss == missi])
  }
  if (nsearch > 10){ # calculations can be slow if there are very many searches
    # ignore intervals beyond the 99th percentile of discovery intervals
    iskip <- min(which(cumsum(ggnm)/sum(ggnm) > 0.99)) + 1
    # cutting off the end introduces a bias.
    # Correct by multiplying the final g's by sum(ggnm)/ggnm[iskip]
    gadj <- sum(ggnm)/sum(ggnm[1:iskip])
  } else {
    iskip <- maxmiss
    gadj <- 1
  }
  ##  b. determine how many simulation draws are needed:
  ##   choose nsim large enough so that the estimated ghat
  ##   is known within 1% [or SE(ghat) < 0.005 * ghat]
  r0 <- GenEst::ppersist(pda = ab0[1], pdb = ab0[2],dist = cpmodel$dist,
    t_arrive0 = 0,
    t_arrive1 = round(Isam),
    t_search = round(Isam)) # need an "if"
  g0 <- f0 * r0
  # estimated variance
  nsim <- 1000
  CPab <- sim_cp0(cpmodel, nsim)
  rr <- as.vector(GenEst::ppersist(pda = CPab[, 1], pdb  =  CPab[, 2],
    dist = cpmodel$dist,
    t_arrive0 = 0, t_arrive1 = round(Isam), t_search = round(Isam)))
  rr[is.na(rr)] <- max(rr[!is.na(rr)])
  Vhatr <- var(rr)
  Vhatf <- (sumstat["p", "sd"])^2
  shatg <- sqrt(f0^2 * Vhatr + r0^2 * Vhatf + Vhatf * Vhatr)
  nsim <- min(round((shatg/(g0 * 0.003))^2), 20000)

  ###2. estimation of g
  if (nsim <= 1 | shatg < 0.001){
    ghat <- sum(ggnm)
    prob_obs <- sum(ggnm)
    Bab <- c(-1, -1)    # this is an indicator of fixed g
    BabAnn <- c(-1, -1)
    BabRaw <- c(-1, -1)
  } else {
    schedule <- cbind(days[ind1], days[ind2], days[ind3])[ind2 >= ind3 - iskip + 1,]
    schedule.index <- cbind(ind1, ind2, ind3)[ind2 >= ind3 - iskip + 1,]
    nmiss <- schedule.index[, 3] - schedule.index[, 2]
    maxmiss <- max(nmiss)
    # searcher efficiencies
    tmppk <- sim_pk0(pkmodel, nsim)
    if (maxmiss  ==  0) pfind.si <- tmppk[, "p"]
    if (maxmiss  ==  1) pfind.si <- cbind(tmppk[, "p"],
      (1 - tmppk[, "p"]) * tmppk[, "k"] * tmppk[, "p"])
    if (maxmiss > 1){
      powk <- array(dim = c(nsim, maxmiss + 1))
      powk[, 1] <- 1
      for (i in 1:maxmiss + 1) powk[,i]  <-  powk[,i - 1] * tmppk[, "k"]
      pfind.si <- tmppk[, "p"] * powk * cbind(1, t(apply( # should be able to use matrixStats here
        1 - (tmppk[, "p"] * powk)[, 1:maxmiss], FUN  =  cumprod, MARGIN  =  1)))
    }

    intxsearch <- unique(
      cbind(schedule[, 2] - schedule[, 1], schedule[, 3] - schedule[, 2]),
      MARGIN = 1)

    CPab <- sim_cp0(cpmodel, nsim)
    # deja vu...should save calcs from earlier to save time?
    ppersu <- GenEst::ppersist(pda = CPab[,1], pdb = CPab[,2],
      dist = cpmodel$dist,
      t_arrive0 = 0,
      t_arrive1 = intxsearch[, 1],
      t_search = intxsearch[, 1] + intxsearch[,2])
    ppersu[is.na(ppersu)]  <-  max(ppersu[!is.na(ppersu)])
    # arrivals
    # if uniform arrivals
    arrvec <- (schedule[, 2] - schedule[, 1])/(nsearch * Isam)

    # add the probabilities
    prob_obs <- numeric(nsim)
    if (maxmiss > 0){
      for (i in 1:dim(schedule)[1]){
        prob_obs <- prob_obs + pfind.si[, nmiss[i] + 1] *
          ppersu[which(
            abs(intxsearch[, 1] - (schedule[i, 2] - schedule[i, 1])) < 0.001 &
            abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 2])) < 0.001),
          ] * arrvec[i]
      }
    } else {
      for (i in 1:dim(schedule)[1]){
        prob_obs <- prob_obs + pfind.si[nmiss[i] + 1] *
          ppersu[which(
            abs(intxsearch[, 1] - (schedule[i, 2] - schedule[i, 1])) < 0.001 &
            abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 2])) < 0.001),
          ] * arrvec[i]
      }
    }
    # g for monitored period
    if (max(prob_obs) * gadj < 1) prob_obs  <-  prob_obs * gadj
    muB <- mean(prob_obs); sig2B <- var(prob_obs)
    Ba <- muB^2/sig2B*(1 - muB) - muB; Bb <- Ba*(1/muB - 1)
    BabRaw <- suppressWarnings(
      MASS::fitdistr(prob_obs, 'beta', # would lognormal (or logit-normal) be better?
        start = list(shape1 = Ba, shape2 = Bb),
        control = list(parscale = c(Ba, Bb))
      )$estimate)
      names(BabRaw) <- c("Ba", "Bb")
     if (!is.null(a)){
      prob_obs  <-  prob_obs * a
      muB <- mean(prob_obs); sig2B <- var(prob_obs)
      Ba <- muB^2/sig2B * (1 - muB) - muB; Bb <- Ba * (1/muB - 1)
      Bab <- suppressWarnings(
        MASS::fitdistr(prob_obs,"beta",
          start = list(shape1 = Ba, shape2 = Bb),
          control = list(parscale = c(Ba, Bb))
        )$estimate
      )
    } else {
      Bab <- c(NA, NA)
    }
    names(Bab) <- names(BabRaw)
    if (!is.null(v)){
      prob_obs <- prob_obs * v
      muB <- mean(prob_obs); sig2B <- var(prob_obs)
      Ba <- muB^2/sig2B * (1 - muB) - muB; Bb <- Ba*(1/muB - 1)
      if (sig2B > 0.00001){
        BabAnn <- suppressWarnings(MASS::fitdistr(prob_obs,'beta',
          start = list(shape1 = Ba, shape2 = Bb),
          control = list(parscale = c(Ba, Bb)))$estimate)
      } else {
        BabAnn <- c(Ba, Bb)
      }
    } else {
      BabAnn <- c(NA, NA)
    }
    names(BabAnn) <- names(Bab)
  }
  out <- list(BabRaw = BabRaw, Bab = Bab, BabAnn = BabAnn,
    pkmodel = pkmodel, cpmodel = cpmodel, a = a, v = v)
  class(out) <- "estg"
  return(out)
}

#' @title Summary statistics for estimated g
#' @param object An \code{\link[=est_g0]{estg}} object
#' @param crlev Credibility level of estimated CI to be returned
#' @param ... additional (optional) arguments passed to \code{rjags::coda.samples}
#' @return summary statistics for estimated \code{g}. \code{searched} is for the
#'  fraction of carcasses falling in the searched area during the monitored
#'  period, \code{site} is area-adjusted to account for carcasses falling outside
#'  the searched area, and \code{full} is further extrapolated to the full year.
#' @export
#'
summary.estg <- function(object, crlev = 0.95, ...){
  x <- object
  lwr <- (1 - crlev)/2
  upr <- 1 - (1 - crlev)/2
  out <- array(dim = c(3, 8), dimnames = list(
    c("searched", "site", "full"), c("Ba", "Bb", "mean",
    paste0(lwr * 100, "%"), "25%", "median", "75%", paste0(upr * 100, "%"))))
  probs <- c((1 - crlev)/2, 0.25, 0.5, 0.75, 1 - (1 - crlev)/2)
  for (i in 1:nrow(out)){
    if (!any(is.na(x[[i]]))){
      a <- x[[i]]["Ba"]
      b <- x[[i]]["Bb"]
      out[i, ] <- c(a, b, a/(a + b), qbeta(probs, shape1 = a, shape2 = b))
    } else {
      out[i, ] <- NA
    }
  }
  attr(out, "a") <- x$a
  attr(out, "v") <- x$v
  return(out)
}