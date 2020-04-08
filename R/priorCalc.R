# Prior calculation using expert classification framework

#AHG b exponent functions-------------------------------------------------------------------
#' Estimate AHG b exponent using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_b <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwsd <- apply(log(Wobs), 1, function(x) sd(x, na.rm = TRUE))

  #expert classification
  temp <- c(0.281116498,
            0.24873163,
            0.233806573,
            0.221609934,
            0.190969495,
            0.186128473,
            0.145874141,
            0.15322105,
            0.405)

  class <- apply(Wobs, 1, classify_func)
  b_hat <- ifelse(class != 100, temp[class], 0.0569 + 0.3822 * lwsd) #repeat by sptial unit
  #global r2: 0.726
}

#' Estimate AHG b lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.01904573,
            0.016895241,
            0.009206385,
            0.009206385,
            0.009206385,
            0.009634045,
            0.008909195,
            0.000182357,
            0.029)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_b <- ifelse(class != 100, temp[class], 0.000182357)
  lowerbound_b <- min(lowerbound_b)
}

#' Estimate AHG b upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.718356399,
            0.718356399,
            0.773757585,
            0.773757585,
            0.773757585,
            0.773757585,
            0.659229465,
            0.659229465,
            0.77)

  class <- apply(Wobs, 1, classify_func)
  upperbound_b <- ifelse(class != 100, temp[class], 0.773757585)
  upperbound_b <- max(upperbound_b)
}

#' Estimate AHG b SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_bSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.127044741,
            0.121791926,
            0.116980495,
            0.120133338,
            0.11851495,
            0.131447085,
            0.123924935,
            0.117431499,
            0.11)

  class <- apply(Wobs, 1, classify_func)
  b_sd <- ifelse(class != 100, temp[class], 0.068077044)
}

#A0 functions---------------------------------------------------------------
#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #Expert classification
  temp <- c(4.235554731,
            4.890349128,
            5.036952602,
            5.347107531,
            5.768320996,
            6.488444764,
            7.222566019,
            8.496990484,
            4.394)

  class <- apply(Wobs, 1, classify_func)
  logA0_hat <- ifelse(class != 100, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' Estimate base cross-sectional area lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(-0.328504067,
            0.385262401,
            -0.192371893,
            -0.192371893,
            -0.192371893,
            0.91027266,
            0.91027266,
            1.545432582,
            0.262)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_A0 <- ifelse(class != 100, exp(temp[class]), exp(-0.328504067))
  lowerbound_A0 <- min(lowerbound_A0)
}

#' Estimate base cross-sectional area upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0 <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(10.20728901,
            10.20728901,
            8.997147152,
            10.20728901,
            10.20728901,
            10.20728901,
            10.20728901,
            11.6483301,
            11.55)

  class <- apply(Wobs, 1, classify_func)
  upperbound_A0 <- ifelse(class != 100, exp(temp[class]), exp(11.6483301))
  upperbound <- max(upperbound_A0)
}

#' Estimate base cross-sectional area SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_A0SD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #expert classification
  temp <- c(1.186970818,
            1.115671401,
            1.139077766,
            1.159659197,
            1.332424151,
            1.420396679,
            1.338002098,
            1.793626478,
            2.285)

  class <- apply(Wobs, 1, classify_func)
  logA0_sd <- ifelse(class != 100, temp[class], 0.58987527)
}

#Bankful Width---------------------------------------------------------
#'Estimate bankful width using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(2.773369667,
            3.231397105,
            3.417726684,
            3.557488626,
            3.732537423,
            4.090671157,
            4.4893103,
            5.145836126,
            3.039)

  class <- apply(Wobs, 1, classify_func)
  logWb_hat <- ifelse(class != 100, temp[class], 0.0037 + 1.0028 * lwbar) #repeat for each sptial unit
  #global r2: 0.984
}

#'Estimate bankful width lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(-0.122732765,
            -0.122732765,
            -0.122732765,
            0.211273379,
            -0.122732765,
            0.461215123,
            0.42199441,
            0.42199441,
            -0.1227)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logWb <- ifelse(class != 100, temp[class], -0.122732765)
  lowerbound_logWb <- min(lowerbound_logWb)
}

#'Estimate bankful width upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogWb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(6.372636963,
            6.372636963,
            5.967171855,
            6.372636963,
            6.372636963,
            6.372636963,
            6.540091608,
            7.006785802,
            6.917)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logWb <- ifelse(class != 100, temp[class], 7.006785802)
  upperbound_logWb <- max(upperbound_logWb)
}

#'Estimate bankful width SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWbSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(0.711826397,
            0.640918356,
            0.621565057,
            0.600591021,
            0.710380783,
            0.74356749,
            0.765907204,
            1.029163996,
            1.284)

  class <- apply(Wobs, 1, classify_func)
  logWb_sd <- ifelse(class != 100, temp[class], 0.137381044)
}

#Bankful depth-------------------------------------------------------------------------
#'Estimate bankful depth using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(-0.855744545,
            -0.746966348,
            -0.749188571,
            -0.666909064,
            -0.467863065,
            -0.036929277,
            0.252025621,
            1.264673325,
            -1.00)

  class <- apply(Wobs, 1, classify_func)
  logDb_hat <- ifelse(class != 100, temp[class], -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar) #repeat for each spatial unit
  #global r2: 0.640
}

#'Estimate bankful depth lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(-3.020024966,
            -3.020024966,
            -2.926934543,
            -2.926934543,
            -2.926934543,
            -2.094164783,
            -2.094164783,
            -2.094164783,
            -3.02)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logDb <- ifelse(class != 100, temp[class], -3.020024966)
  lowerbound_logDb <- min(lowerbound_logDb)
}

#'Estimate bankful depth upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogDb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(1.7933029,
            1.7933029,
            1.7933029,
            1.7933029,
            1.708874134,
            2.009770826,
            2.009770826,
            3.309358647,
            2.572)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 100, temp[class], 3.309358647)
  upperbound_logDb <- max(upperbound_logDb)
}

#'Estimate bankful depth SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDbSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(0.598955804,
            0.608728355,
            0.643855445,
            0.684081848,
            0.746547309,
            0.806479131,
            0.733313884,
            0.927305358,
            1.147)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 100, temp[class], 0.576212733)
}

#Channel shape----------------------------------------------
#'Estimate channel shape using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #supervised classification
  temp <- c(0.325295247,
            0.559767277,
            0.618153535,
            0.656156611,
            0.796844318,
            0.845216799,
            1.106992587,
            0.875041483,
            -0.249)

  class <- apply(Wobs, 1, classify_func)
  logr_hat <- ifelse(class != 100, temp[class], 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar) #repeat for each spatial unit
  #Global r2: 0.421
}

#'Estimate channel shape lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(-1.694462984,
            -1.694462984,
            -1.694462984,
            -1.504526738,
            -1.694462984,
            -1.59134178,
            -1.59134178,
            -2.580471126,
            -2.58)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logr <- ifelse(class != 100, temp[class], -2.580471126)
  lowerbound_logr <- min(lowerbound_logr)
}

#'Estimate channel shape upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogr <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(2.750427051,
            2.750427051,
            3.660494543,
            3.660494543,
            3.75765632,
            3.75765632,
            3.885278632,
            8.037716276,
            0)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logr <- ifelse(class != 100, temp[class], 8.037716276)
  upperbound_logr <- max(upperbound_logr)
}

#'Estimate channel shape SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logrSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #expert classification
  temp <- c(0.611235593,
            0.65753785,
            0.649325058,
            0.690260039,
            0.73068479,
            0.874593274,
            0.916480631,
            1.043382513,
            0.412)

  class <- apply(Wobs, 1, classify_func)
  logr_sd <- ifelse(class != 100, temp[class], 0.67332688)
}

#Manning's n------------------------------------------------------------------------
#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_logn <- function(Sobs, Wobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #expert classification
  temp <- c(-3.34557283,
            -3.312687279,
            -3.323656967,
            -3.46372292,
            -3.417279106,
            -3.470045772,
            -3.225956303,
            -3.400259437,
            -3.41)

  class <- apply(Wobs, 1, classify_func)
  logn_hat <- ifelse(class != 100, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-7.628253825,
            -7.369108582,
            -7.720003485,
            -7.385380143,
            -8.636508843,
            -8.741518112,
            -8.941359471,
            -12.60217578,
            log(0.01))

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logn <- ifelse(class != 100, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(0.726125533,
            -0.270337583,
            -0.013112528,
            0.030500965,
            0.207297661,
            0.298670281,
            0.611600396,
            3.047771386,
            log(0.05))

  class <- apply(Wobs, 1, classify_func)
  upperbound_logn <- ifelse(class != 100, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(1.154479143,
            1.129036566,
            1.167466126,
            1.200156724,
            1.223989206,
            1.206875933,
            1.306485481,
            1.484896588,
            1.23)

  class <- apply(Wobs, 1, classify_func)
  logn_sd <- ifelse(class != 100, temp[class], 0.761673112)
}


# Prior calculation using unsupervised framework------------------------------------------------------------------------------
#class 100 are noisy rivers

#AHG b exponent functions-------------------------------------------------------------------
#' Estimate AHG b exponent using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_b_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwsd <- apply(log(Wobs), 1, function(x) sd(x, na.rm = TRUE))

  #unsupervised classification
  temp <- c(0.281116498,
            0.24873163,
            0.233806573,
            0.221609934,
            0.190969495,
            0.186128473,
            0.145874141,
            0.15322105,
            0.405)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  b_hat <- ifelse(class != 101, temp[class], 0.0569 + 0.3822 * lwsd) #repeat by sptial unit
  #global r2: 0.726
}

#' Estimate AHG b lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.01904573,
            0.016895241,
            0.009206385,
            0.009206385,
            0.009206385,
            0.009634045,
            0.008909195,
            0.000182357,
            0.029)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_b <- ifelse(class != 101, temp[class], 0.000182357)
  lowerbound_b <- min(lowerbound_b)
}

#' Estimate AHG b upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.718356399,
            0.718356399,
            0.773757585,
            0.773757585,
            0.773757585,
            0.773757585,
            0.659229465,
            0.659229465,
            0.77)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_b <- ifelse(class != 101, temp[class], 0.773757585)
  upperbound_b <- max(upperbound_b)
}

#' Estimate AHG b SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_bSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.127044741,
            0.121791926,
            0.116980495,
            0.120133338,
            0.11851495,
            0.131447085,
            0.123924935,
            0.117431499,
            0.11)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  b_sd <- ifelse(class != 101, temp[class], 0.068077044)
}

#A0 functions---------------------------------------------------------------
#' Estimate base cross-sectional area using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(4.235554731,
            4.890349128,
            5.036952602,
            5.347107531,
            5.768320996,
            6.488444764,
            7.222566019,
            8.496990484,
            4.394)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logA0_hat <- ifelse(class != 101, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
  #global r2: 0.907
}

#' Estimate base cross-sectional area lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(-0.328504067,
            0.385262401,
            -0.192371893,
            -0.192371893,
            -0.192371893,
            0.91027266,
            0.91027266,
            1.545432582,
            0.262)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_A0 <- ifelse(class != 101, exp(temp[class]), exp(-0.328504067))
  lowerbound_A0 <- min(lowerbound_A0)
}

#' Estimate base cross-sectional area upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundA0_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(10.20728901,
            10.20728901,
            8.997147152,
            10.20728901,
            10.20728901,
            10.20728901,
            10.20728901,
            11.6483301,
            11.55)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_A0 <- ifelse(class != 101, exp(temp[class]), exp(11.6483301))
  upperbound_A0 <- max(upperbound_A0)
}

#' Estimate base cross-sectional area SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_A0SD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #unsupervised classification
  temp <- c(1.186970818,
            1.115671401,
            1.139077766,
            1.159659197,
            1.332424151,
            1.420396679,
            1.338002098,
            1.793626478,
            2.285)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logA0_sd <- ifelse(class != 101, temp[class], 0.58987527)
}

#Bankful Width---------------------------------------------------------
#'Estimate bankful width using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(2.773369667,
            3.231397105,
            3.417726684,
            3.557488626,
            3.732537423,
            4.090671157,
            4.4893103,
            5.145836126,
            3.039)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logWb_hat <- ifelse(class != 101, temp[class], 0.0037 + 1.0028 * lwbar) #repeat for each sptial unit
  #global r2: 0.984
}

#'Estimate bankful width lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-0.122732765,
            -0.122732765,
            -0.122732765,
            0.211273379,
            -0.122732765,
            0.461215123,
            0.42199441,
            0.42199441,
            -0.1227)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logWb <- ifelse(class != 101, temp[class], -0.122732765)
  lowerbound_logWb <- min(lowerbound_logWb)
}

#'Estimate bankful width upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogWb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(6.372636963,
            6.372636963,
            5.967171855,
            6.372636963,
            6.372636963,
            6.372636963,
            6.540091608,
            7.006785802,
            6.917)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logWb <- ifelse(class != 101, temp[class], 7.006785802)
  upperbound_logWb <- max(upperbound_logWb)
}

#'Estimate bankful width SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logWbSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(0.711826397,
            0.640918356,
            0.621565057,
            0.600591021,
            0.710380783,
            0.74356749,
            0.765907204,
            1.029163996,
            1.284)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logWb_sd <- ifelse(class != 101, temp[class], 0.137381044)
}

#Bankful depth-------------------------------------------------------------------------
#'Estimate bankful depth using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-0.855744545,
            -0.746966348,
            -0.749188571,
            -0.666909064,
            -0.467863065,
            -0.036929277,
            0.252025621,
            1.264673325,
            -1.00)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logDb_hat <- ifelse(class != 101, temp[class], -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar) #repeat for each spatial unit
  #global r2: 0.640
}

#'Estimate bankful depth lower bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-3.020024966,
            -3.020024966,
            -2.926934543,
            -2.926934543,
            -2.926934543,
            -2.094164783,
            -2.094164783,
            -2.094164783,
            -3.02)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logDb <- ifelse(class != 101, temp[class], -3.020024966)
  lowerbound_logDb <- min(lowerbound_logDb)
}

#'Estimate bankful depth upper bound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogDb_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(1.7933029,
            1.7933029,
            1.7933029,
            1.7933029,
            1.708874134,
            2.009770826,
            2.009770826,
            3.309358647,
            2.572)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logDb <- ifelse(class != 101, temp[class], 3.309358647)
  upperbound_logDb <- max(upperbound_logDb)
}

#'Estimate bankful depth SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logDbSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(0.598955804,
            0.608728355,
            0.643855445,
            0.684081848,
            0.746547309,
            0.806479131,
            0.733313884,
            0.927305358,
            1.147)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logDb <- ifelse(class != 101, temp[class], 0.576212733)
  upperbound_logDb <- max(upperbound_logDb)
}

#Channel shape----------------------------------------------
#'Estimate channel shape using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)

  #supervised classification
  temp <- c(0.325295247,
            0.559767277,
            0.618153535,
            0.656156611,
            0.796844318,
            0.845216799,
            1.106992587,
            0.875041483,
            -0.249)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logr_hat <- ifelse(class != 101, temp[class], 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar) #repeat for each spatial unit
  #Global r2: 0.421
}

#'Estimate channel shape lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundlogr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-1.694462984,
            -1.694462984,
            -1.694462984,
            -1.504526738,
            -1.694462984,
            -1.59134178,
            -1.59134178,
            -2.580471126,
            -2.58)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logr <- ifelse(class != 101, temp[class], -2.580471126)
  lowerbound_logr <- min(lowerbound_logr)
}

#'Estimate channel shape upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundlogr_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(2.750427051,
            2.750427051,
            3.660494543,
            3.660494543,
            3.75765632,
            3.75765632,
            3.885278632,
            8.037716276,
            0)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logr <- ifelse(class != 101, temp[class], 8.037716276)
  upperbound_logr <- max(upperbound_logr)
}

#'Estimate channel shape SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_logrSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lwbar <- apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- apply(log(Wobs), 1, sd, na.rm = TRUE)
  #unsupervised classification
  temp <- c(0.611235593,
            0.65753785,
            0.649325058,
            0.690260039,
            0.73068479,
            0.874593274,
            0.916480631,
            1.043382513,
            0.412)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logr_sd <- ifelse(class != 101, temp[class], 0.67332688)
}

#Manning's n------------------------------------------------------------------------
#'Estimate manning's n using bam data
#'
#' @param Sobs Observed S, as a space-down, time-across matrix
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_logn_unsupervised <- function(Sobs, Wobs) {
  Sobs[Sobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  lsbar <- apply(log(Sobs), 1, mean, na.rm = TRUE)
  #unsupervised classification
  temp <- c(-3.34557283,
            -3.312687279,
            -3.323656967,
            -3.46372292,
            -3.417279106,
            -3.470045772,
            -3.225956303,
            -3.400259437,
            -3.41)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logn_hat <- ifelse(class != 101, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(-7.628253825,
            -7.369108582,
            -7.720003485,
            -7.385380143,
            -8.636508843,
            -8.741518112,
            -8.941359471,
            -12.60217578,
            log(0.01))

  class <- apply(Wobs, 1, classify_func_unsupervised)
  lowerbound_logn <- ifelse(class != 101, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(0.726125533,
            -0.270337583,
            -0.013112528,
            0.030500965,
            0.207297661,
            0.298670281,
            0.611600396,
            3.047771386,
            log(0.05))

  class <- apply(Wobs, 1, classify_func_unsupervised)
  upperbound_logn <- ifelse(class != 101, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD_unsupervised <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #unsupervised classification
  temp <- c(1.154479143,
            1.129036566,
            1.167466126,
            1.200156724,
            1.223989206,
            1.206875933,
            1.306485481,
            1.484896588,
            1.23)

  class <- apply(Wobs, 1, classify_func_unsupervised)
  logn_sd <- ifelse(class != 101, temp[class], 0.761673112)
}
