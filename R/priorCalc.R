# Prior calculation

#' Estimate AHG b exponent using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_b <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwsd <- apply(log(Wobs), 1, function(x) sd(x, na.rm = TRUE))

  b_hat <- 0.0569 + 0.3822 * lwsd #r2: 0.726
}

#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  logA0hat <- -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd #r2: 0.907

  logA0hat
}

#'Estimate bankful width using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)

  logWbhat <- 0.0037 + 1.0028 * lwbar #r2: 0.984
  logWbhat
}

#'Estimate bankful depth using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  logDbhat <- -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar #r2: 0.640
  logDbhat
}

#'Estimate channel shape using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  logrhat <- 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar #r2: 0.421
}

#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @export
estimate_logn <- function(Sobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)

  lognhat <- -0.1636 + 0.4077 * lsbar #r2: 0.631
  lognhat
}
