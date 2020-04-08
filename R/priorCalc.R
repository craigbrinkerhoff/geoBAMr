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
  temp <- c(0.227633914,
            0.24121913,
            0.204897294,
            0.204451567,
            0.192887656,
            0.188263846,
            0.211646902,
            0.176051924,
            0.164989336,
            0.14444077,
            0.179925708,
            0.151049254,
            0.107294865,
            0.125720128,
            0.117940326,
            0.404257196)

  class <- apply(Wobs, 1, classify_func)
  b_hat <- ifelse(class != 17, temp[class], 0.0569 + 0.3822 * lwsd) #repeat by sptial unit
  #global r2: 0.726
}

#' Estimate AHG b lowerbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_lowerboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.009634045,
            0.009206385,
            0.010378626,
            0.00984031,
            0.000182357,
            0.008909195,
            0.025751932,
            0.005703884,
            0.017018127,
            0.013954247,
            0.019210795,
            0.010974071,
            0.016252382,
            0.017584834,
            0.004872167,
            0.029379248)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_b <- ifelse(class != 17, temp[class], 0.000182357)
  lowerbound_b <- min(lowerbound_b)
}

#' Estimate AHG b upperbound using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_upperboundb <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.401485695,
            0.450546257,
            0.442897059,
            0.427490668,
            0.415703708,
            0.388788209,
            0.422853699,
            0.41859221,
            0.397420318,
            0.432497404,
            0.404125905,
            0.413263659,
            0.361701335,
            0.328343393,
            0.348032765,
            0.773757585)

  class <- apply(Wobs, 1, classify_func)
  upperbound_b <- ifelse(class != 17, temp[class], 0.773757585)
  upperbound_b <- max(upperbound_b)
}

#' Estimate AHG b SD using bam data
#'
#' @param Wobs Observed W,as a space-down, time-across matrix.
#' @export
estimate_bSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  temp <- c(0.094625319,
            0.103986506,
            0.110231244,
            0.103910665,
            0.101259677,
            0.098511991,
            0.097920891,
            0.0962566,
            0.100083519,
            0.10412018,
            0.099927194,
            0.098731293,
            0.099282596,
            0.075781266,
            0.091252051,
            0.112180583)

  class <- apply(Wobs, 1, classify_func)
  b_sd <- ifelse(class != 17, temp[class], 0.068077044)
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
  temp <- c(3.723280881,
            4.497073804,
            4.753590191,
            4.990432587,
            4.911328517,
            5.350615603,
            5.422742508,
            5.523458921,
            5.774532256,
            6.446836611,
            6.527953643,
            6.873163834,
            7.102499356,
            8.007965013,
            8.937204637,
            4.432006567)

  class <- apply(Wobs, 1, classify_func)
  logA0_hat <- ifelse(class != 17, temp[class], -0.2918 + 1.6930 * lwbar - 0.1887 * lwsd) #repeat for each sptial unit
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
  temp <- c(0.732367894,
            1.508511994,
            0.91027266,
            2.517696473,
            1.199964783,
            2.681021529,
            2.148850993,
            1.545432582,
            2.415913778,
            3.106826321,
            3.874321138,
            2.694627181,
            3.696351469,
            3.593194204,
            4.043928076,
            0.262364264)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_A0 <- ifelse(class != 17, exp(temp[class]), exp(4.540631665))
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
  temp <- c(7.640123173,
            7.355641103,
            8.997147152,
            9.164296433,
            8.554488976,
            9.417354541,
            7.677863501,
            8.144679183,
            7.863266724,
            8.793308627,
            8.776475789,
            9.014325488,
            8.78186249,
            9.61580548,
            11.6483301,
            11.55214618)

  class <- apply(Wobs, 1, classify_func)
  upperbound_A0 <- ifelse(class != 17, exp(temp[class]), 1000000)
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
  temp <- c(1.446820524,
            1.350055707,
            1.363875688,
            1.260787706,
            1.306619211,
            1.315826017,
            1.190542029,
            1.337220271,
            1.21047899,
            1.359096608,
            1.245863462,
            1.287297345,
            1.08535437,
            1.154319081,
            1.5575699,
            2.272342301)

  class <- apply(Wobs, 1, classify_func)
  logA0_sd <- ifelse(class != 17, temp[class], 0.58987527)
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
  temp <- c(2.476118144,
            2.864001065,
            3.103015939,
            3.249308032,
            3.284178964,
            3.371669039,
            3.574693905,
            3.664586762,
            3.683922384,
            4.002696788,
            4.032912323,
            4.357733942,
            4.436574004,
            4.921742348,
            5.287893051,
            3.032064203)

  class <- apply(Wobs, 1, classify_func)
  logWb_hat <- ifelse(class != 17, temp[class], 0.0037 + 1.0028 * lwbar) #repeat for each sptial unit
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
  temp <- c(0.316633894,
            0.814036498,
            0.461215123,
            1.377505855,
            0.744077909,
            2.070653036,
            1.73032723,
            0.42199441,
            1.756995477,
            2.396075436,
            2.293796587,
            2.268873179,
            2.99460664,
            1.852305657,
            2.585317436,
            -0.122732765)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logWb <- ifelse(class != 17, temp[class], 3.089222617)
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
  temp <- c(4.016563185,
            4.410978457,
            4.833579847,
            5.094823245,
            5.029162599,
            5.216130696,
            4.5361416,
            5.091077922,
            5.194372515,
            5.100415058,
            5.372171735,
            5.860073719,
            5.521860838,
            5.749870579,
            7.006785802,
            6.917259966)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logWb <- ifelse(class != 17, temp[class], 10)
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
  temp <- c(0.750069561,
            0.738069364,
            0.710662053,
            0.669257108,
            0.685680677,
            0.630380021,
            0.563796126,
            0.812622517,
            0.691428554,
            0.654155675,
            0.587954845,
            0.728133987,
            0.540567805,
            0.680678526,
            0.889425931,
            1.268865384)

  class <- apply(Wobs, 1, classify_func)
  logWb_sd <- ifelse(class != 17, temp[class], 0.137381044)
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
  temp <- c(-1.035199401,
            -0.738778811,
            -0.759250143,
            -0.714358207,
            -0.834152133,
            -0.629032654,
            -0.471723514,
            -0.547540836,
            -0.480357291,
            0.130087152,
            -0.07180952,
            0.15467964,
            0.299392088,
            0.579044554,
            1.48611948,
            -0.974789558)

  class <- apply(Wobs, 1, classify_func)
  logDb_hat <- ifelse(class != 17, temp[class], -2.6189 - 0.2436 * lwsd + 0.6854 * lwbar) #repeat for each spatial unit
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
  temp <- c(-2.326877786,
            -1.788516004,
            -2.05494407,
            -1.892425141,
            -2.057110918,
            -1.947205432,
            -1.915682003,
            -1.514919292,
            -1.730529988,
            -1.638497383,
            -1.809169822,
            -1.622549192,
            -1.649478962,
            -0.83244782,
            -0.924234286,
            -2.379581849)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logDb <- ifelse(class != 17, temp[class], -0.951962577)
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
  temp <- c(1.49936159,
            1.045713084,
            1.7933029,
            1.725376856,
            1.127811443,
            1.827833367,
            1.174682984,
            1.449780646,
            1.260032945,
            1.708874134,
            1.765116599,
            1.327316111,
            1.639749354,
            2.102212505,
            3.309358647,
            2.571816749)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 17, temp[class], 5)
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
  temp <- c(0.754179034,
            0.719470145,
            0.751836561,
            0.699278157,
            0.699865238,
            0.771195464,
            0.69921735,
            0.665457667,
            0.662592139,
            0.810898004,
            0.764674068,
            0.719504763,
            0.690395904,
            0.664391944,
            0.8459272,
            1.132642792)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logDb <- ifelse(class != 17, temp[class], 0.576212733)
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

  #expert classification
  temp <- c(0.479582219,
            0.475424989,
            0.742030366,
            0.747458374,
            0.810066771,
            0.764035146,
            0.728537155,
            0.914703855,
            1.012946569,
            1.012082534,
            0.964521619,
            1.007259927,
            1.368493142,
            1.217111712,
            1.067996942,
            -0.247245593)

  class <- apply(Wobs, 1, classify_func)
  logr_hat <- ifelse(class != 17, temp[class], 1.4241 - 1.9097 * lwsd + 0.0420 * lwbar) #repeat for each spatial unit
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
  temp <- c(0.010687262,
            0.003572302,
            0.030800262,
            0.055237814,
            0.075268048,
            0.005163827,
            0.003031095,
            0.140133893,
            0.014043795,
            0.007552884,
            0.033288292,
            0.017194676,
            0.081171573,
            0.128060852,
            0.020365892,
            -2.580471126)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logr <- ifelse(class != 17, temp[class], -2.580471126)
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
  temp <- c(2.659590969,
            3.660494543,
            3.75765632,
            3.809886281,
            8.037716276,
            3.885278632,
            3.141534762,
            4.058323672,
            2.894177697,
            3.145908752,
            3.291615219,
            3.759676814,
            3.688940701,
            3.293279235,
            4.028976176,
            -0.002840899)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logr <- ifelse(class != 17, temp[class], 8.037716276)
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
  temp <- c(0.574788166,
            0.686443448,
            0.735995097,
            0.660600833,
            1.152463807,
            0.751672132,
            0.646178286,
            0.721545901,
            0.678734433,
            0.73732522,
            0.714799269,
            0.783012838,
            0.824128369,
            0.666709963,
            0.741805733,
            0.418335751)

  class <- apply(Wobs, 1, classify_func)
  logr_sd <- ifelse(class != 17, temp[class], 0.67332688)
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
  temp <- c(-2.956882373,
            -3.218074003,
            -3.073438472,
            -3.316850463,
            -3.318215717,
            -3.15335243,
            -3.456243594,
            -3.545592423,
            -3.240085716,
            -3.401877538,
            -3.269861085,
            -3.372948626,
            -3.404202242,
            -3.274729636,
            -3.405497386,
            -3.318138928)

  class <- apply(Wobs, 1, classify_func)
  logn_hat <- ifelse(class != 17, temp[class], -0.1636 + 0.4077 * lsbar) #repeat for each sptial unit
  #Global r2: 0.631
}

#'Estimate manning's n lowerbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lowerboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-5.678333858,
            -5.467593416,
            -5.67856698,
            -5.768715685,
            -5.682940688,
            -5.841321589,
            -5.56434293,
            -5.957755533,
            -5.646550723,
            -6.064509006,
            -6.029056118,
            -6.405084634,
            -6.678294272,
            -6.211213087,
            -6.06854241,
            -6.246537799)

  class <- apply(Wobs, 1, classify_func)
  lowerbound_logn <- ifelse(class != 17, temp[class], log(0.01))
  lowerbound_logn <- min(lowerbound_logn)
}

#'Estimate manning's n upperbound using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_upperboundlogn <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(-0.377797559,
            -1.024617813,
            -0.469297566,
            -1.330691863,
            -1.17053266,
            -1.181653683,
            -0.819808414,
            -1.085062606,
            -0.540120016,
            -1.349783005,
            -0.325687448,
            -0.025491812,
            -0.554602291,
            0.501441352,
            1.930473909,
            0.094113848)

  class <- apply(Wobs, 1, classify_func)
  upperbound_logn <- ifelse(class != 17, temp[class], log(0.05))
  upperbound_logn <- max(upperbound_logn)
}

#'Estimate manning's n SD using bam data
#'
#' @param Wobs Observed W, as a space-down, time-across matrix
#' @export
estimate_lognSD <- function(Wobs) {
  Wobs[Wobs <= 0] <- NA # I replaced missing values with 0 so Stan will accept
  #expert classification
  temp <- c(1.069101032,
            1.004209047,
            1.051391594,
            1.050420395,
            1.103812115,
            1.050696519,
            1.205687421,
            1.229666517,
            1.101795647,
            1.194241604,
            1.268715794,
            1.260820049,
            1.199816945,
            1.345333569,
            1.724262114,
            1.122652969)

  class <- apply(Wobs, 1, classify_func)
  logn_sd <- ifelse(class != 17, temp[class], 0.761673112)
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
