# Searcer Efficiency (SE0) ------------------------------
#' A template for summarized searcher efficiency data with the number of
#'  carcasses available and the number discovered for N = 12 search occasions
#'
#' @format A list with 2 numeric N-vectors with numbers of:
#' \describe{
#'   \item{n}{searcher efficiency trial carcasses available}
#'   \item{y}{carcasses discovered}
#' }
"pkdata0"

# Carcass Persistence (CP0) ------------------------------
#' A template for carcass persistence data with interval-censored carcass
#'  persistence times
#'
#' @format A matrix with 2 columns bracketing persistence times (in days since
#'  carcass placement) of each carcass:
#' \describe{
#'   \item{CPmin}{the last time the carcass was observed}
#'   \item{CPmax}{the first time the carcass was noted missing}
#' }
"cpdata0"

# Search Schedule (CP0) ------------------------------
#' A template for search schedule data
#'
#' @format A numeric vector with the times when searches were conducted, with
#'  days[1] = 0:
#' \describe{
#'   \item{days0}{numeric vector of search times}
#' }
"days0"