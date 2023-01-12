#' @title Fit pk searcher efficiency models
#'
#' @description Searcher efficiency is modeled as a function of the number of
#'  times a carcass has been missed in previous searches and any number of
#'  covariates.
#' @details
#'  The probability of finding a carcass that is present at the time of
#'  search is \code{p} on the first search after carcass arrival and is
#'  assumed to decrease by a factor of \code{k} each time the carcass is
#'  missed in searches.
#' @param pkdata Search trial data entered in a list of N-vectors, $n and $y,
#'  indicating the number of carcasses available and the number discovered in
#'  searcher efficiency field trials in which carcasses were available for
#'  discovery. [NOTE: In earlier versions of \code{eoa}, the vectors were $M and
#'  $X. The names have been changed to avoid confusion with the \code{M} and
#'  \code{X} for total mortality and carcasses discovered carcass survey.]
#' @param kFixed If trial carcassses are available for discovery for one search
#'  and data are insufficient for estimating \code{k}, a fixed, assumed value must
#'  be entered for \code{k}.
#' @param n.iter number of iterations to use in updating the JAGS model for
#'  \eqn{p} and \eqn{k}
#' @param ... Other parameters that may be used in called functions (esp.
#'  \code{burn} for updating the JAGS function)
#' @return A list with an nsim x 2 matrix of simulated \code{p} and \code{k}
#'  values
#'  the joint posterior for SE.
#'
#' @export
#'
est_pk0 <- function(pkdata, kFixed = NULL, n.iter = 1000, ...){
  if ("list" %in% class(pkdata)){ # legacy pk format required
    ## error-checking: (NOTE: error-checking will be done outside this function in future)
    if (!all(c("n","y") %in% names(pkdata)))
      stop("pkdata must be list with elements n and y")
    if (!is.vector(pkdata$n) | !is.vector(pkdata$y))
      stop("pkdata$n and pkdata$y must be vectors")
    if (length(pkdata$n) != length(pkdata$y))
      stop("pkdata$n and pkdata$y must be vectors of with the same length")
    if (!is.numeric(pkdata$n) || !is.numeric(pkdata$y))
      stop("pkdata$n and pkdata$y must be numeric")
    with (pkdata, {
      if (sum(n < 1) > 0 || sum(abs(round(n) - n)) > 0.00001)
        stop("pkdata$n must be a vector of positive integers")
      if (sum(y < 0) > 0 || sum(abs(round(y) - y)) > 0.00001)
        warning("pkdata$y must be a vector of non-negative integers")
    })
    if (any(pkdata$y > pkdata$n)){
      stop(paste0("More carcasses found than available in search ",
        min(which(pkdata$y > pkdata$n)), " of SE trials."))
    }
    pkdata$N <- length(pkdata$n) # number of search occasions
  } else if (!"list" %in% class(pkdata)){
    stop("pkdata must be list with elements n and y")
  }
  out <- list()
  if (pkdata$N == 1 & is.null(kFixed)){
    out$type = "pOnly"
    out$model  <- c(pBa = pkdata$y + 0.5, pBb = pkdata$n - pkdata$y + 0.5, k = NA)
  } else if (!is.null(kFixed)){
    if (!is.numeric(kFixed) || !is.finite(kFixed) || kFixed > 1 || kFixed < 0)
      stop ("kFixed must be in interval [0, 1] or NULL")
    out$type <- "kFixed"
    out$model <- c(pBa = pkdata$y[1] + 0.5, pBb = pkdata$n[1] - pkdata$y[1] + 0.5, k = kFixed)
   } else {
    out$model <- rjags::jags.model(
      textConnection("
        model {
          for (i in 1:N) {
            y[i] ~ dbin(p*k^(i-1), n[i])
          }
          p ~ dbeta(1/2, 1/2)
          k ~ dbeta(1, 1)
        }
        "),
      data = pkdata,
      inits = with(pkdata, list(
        p = y[1]/n[1],
        k = max(min((y[2]/n[2])/(y[1]/n[1]), .99), .01))
      ), quiet = TRUE, ...
    )
    update(out$model, n.iter = n.iter, progress.bar = NULL) # simple & unlikely to need more burn-in
    out$type <- "pk"
  }
  out$data <- pkdata
  class(out) <- "estpk"
  return(out)
}

#' @title Summary statistics for estimated p and p parameters
#'
#' @param object An \code{\link[=est_pk0]{estpk}} object
#' @param crlev Credibility level of estimated CI to be returned
#' @param n.iter Number of iterations of the JAGS model for estimating the joint
#'  posterior distribution of \code{p} and \code{k} (relevant only if
#'  \code{object$type == "pk"}).
#' @param ... additional (optional) arguments passed to \code{rjags::coda.samples}
#' @return array of summary statistics for \code{p} and \code{k}
#'
#' @export
#'
summary.estpk <- function(object, crlev = 0.95, n.iter = 10000, ...){
  x <- object
  lwr <- (1 - crlev)/2
  upr <- 1 - (1 - crlev)/2
  out <- array(dim = c(2, 7), dimnames = list(c("p", "k"), c("mean", "sd",
    paste0(lwr * 100, "%"), "25%", "median", "75%", paste0(upr * 100, "%"))))
  if (x$type == "pk"){
    tmp <- as.matrix(rjags::coda.samples(x$model, n.iter = n.iter,
      variable.names = c("p", "k"), progress.bar = NULL, ...)[[1]])[, 2:1]
    for (i in rownames(out)){
      out[i, ] <- c(mean(tmp[, i]), sd(tmp[, i]), quantile(tmp[, i],
        probs = c((1 - crlev)/2, 0.25, 0.5, 0.75, 1 - (1 - crlev)/2),
        names = FALSE
      ))
    }
  } else {
    a <- x$model["pBa"]
    b <- x$model["pBb"]
    out["p", ] <- c(a/(a + b), sqrt(a * b/((a + b)^2 * (a + b + 1))),
      qbeta(c(lwr, 0.25, 0.50, 0.75, upr), shape1 = a, shape2 = b))
    out["k", ] <- x$model["k"]
  }
  return(out)
}
#' @title Simulate pk parameters from model
#'
#' @description Simple simulation
#' @param pkmodel A model returned from \code{est_pk0}
#' @param nsim Number of simulation reps for estimating the joint posterior
#'  distribution of \code{p} and \code{k}.
#' @return An nsim x 2 matrix of simulated \code{p} and \code{k} values
#'  the joint posterior for SE.
#'
#' @export
#'
sim_pk0 <- function(pkmodel, nsim = 1000){
  if (pkmodel$type == "pk"){
    return(as.matrix(rjags::coda.samples(pkmodel$model,
      variable.names = c('p','k'), n.iter = nsim, progress.bar = NULL)[[1]][, 2:1]
    ))
  }
  if (pkmodel$type %in% c("pOnly", "kFixed")){
    return(matrix(
      c(rbeta(nsim, shape1 = pkmodel$model["pBa"], shape2 = pkmodel$model["pBb"]),
        rep(pkmodel$model["k"], nsim)),
      ncol = 2, dimnames = list(NULL, c("p", "k"))
    ))
  }
  return("error")
}


