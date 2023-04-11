
# Using advice from read-and-delete-me:
# "Be sure to add useDynLib(mypackage, .registration = TRUE) to the NAMESPACE file
# which you can do by placing the line   #' @useDynLib rstanarm, .registration = TRUE
# in one of your .R files
# see rstanarm's 'rstanarm-package.R' file"

#' Preprocess data for BAM estimation
#'
#' Produces a bamdata object that can be passed to bam_estimate function
#'
#' @useDynLib geoBAMr, .registration = TRUE
#' @param w Matrix (or data frame) of widths: time as columns, space as rows
#' @param s Matrix (or data frame) of slopes: time as columns, space as rows
#' @param dA Matrix of area above base area: time as columns, space as rows
#' @param Qhat Vector of Q estimates. Needed to create prior on Q.
#' @param variant Which geoBAM variant to use. Options are "manning_amhg" (default),
#'   "manning", or "amhg".
#' @param max_xs Maximum number of cross-sections to allow in data. Used to reduce
#'   sampling time. Defaults to 30.
#' @param seed RNG seed to use for sampling cross-sections, if nx > max_xs.
#' @export

bam_data <- function(w,
                     s,
                     dA = NULL,
                     Qhat,
                     variant = c("manning_amhg", "manning", "amhg"),
                     max_xs = 30L,
                     seed = NULL) {
  force(Qhat)


  if(is.null(dA)==0 && variant=='amhg'){
    stop('AMHG-only geoBAM cannot process a dA matrix!')}
  if(is.null(dA)==1 && variant!='amhg'){
    stop('Mannings-based geoBAM requires a dA matrix!')}
 # w_test <- apply(w,2, function(x){sum(x, na.rm=T)})  #make sure at least 3 observations for width sd calculation for priors
#  if(any(w_test<2)){
#    stop('Need at least 3 non-NA width observations to calculate some priors!')}
  w_test <- apply(w,2, function(x){sum(x > 0, na.rm=T)})  #make sure at least 3 observations for width sd calculation for priors
  if(any(w_test<3)){
    stop('Need at least 3 non-NA width observations to calculate some priors!')}

  manning_ready <- !is.null(s) && !is.null(dA)
  if (!manning_ready) {
    dA <- matrix(1, nrow = nrow(w), ncol = ncol(w))
  }

  datalist <- list(Wobs = w,
                Sobs = s,
                dAobs = dA,
                logQ_hat = log(Qhat))

  datalist <- bam_check_args(datalist)
  datalist <- bam_check_nas(datalist, variant)

  nx <- nrow(datalist$Wobs)
  nt <- ncol(datalist$Wobs)

  out <- structure(c(datalist,
                        nx = nx,
                        nt = nt),
                   manning_ready = manning_ready,
                   class = c("bamdata"))

  if (nx > max_xs)
    out <- sample_xs(out, n = max_xs, seed = seed)

  out
}

#' Performs the following checks:
#' - types:
#'     - logQ_hat is numeric vector
#'     - everything else matrix
#' - dimensions:
#'     - all matrices have same dims
#'     - logQ_hat has length equal to ncol of matrices
#'
#' @param datalist A list of BAM data inputs
bam_check_args <- function(datalist) {

  dA <- datalist$dA
  logQ_hat <- datalist$logQ_hat
  matlist <- datalist[names(datalist) != "logQ_hat"]

  # Check types
  if (!(is(logQ_hat, "numeric") && is(logQ_hat, "vector")))
    stop("Qhat must be a numeric vector.\n")
  if (!all(vapply(matlist, is, logical(1), "matrix")))
    stop("All data must be a supplied as a matrix.\n")

  # Check dims
  nr <- nrow(matlist[[1]])
  nc <- ncol(matlist[[1]])
  if (!(all(vapply(matlist, nrow, 0L) == nr) &&
        all(vapply(matlist, ncol, 0L) == nc)))
    stop("All data must have same dimensions.\n")
  if (!length(logQ_hat) == nc)
    logQ_hat <- rep(logQ_hat, length.out = nc)

  out <- c(matlist, list(logQ_hat = logQ_hat))

  out
}

#' Add missing-data inputs to data list
#'
#' Binary matrices indicating where data are/aren't missing are
#' added to the data list. This is required in order to run
#' ragged-array data structures in the stanfile.
#'
#' Previously this function omitted any times with missing data,
#' but now that ragged arrays are accommodated in the stanfile the
#' operations are entirely different.
#'
#' @param datalist a list of BAM inputs
#' @param variant which geoBAM is being run
#' @importFrom stats median
bam_check_nas <- function(datalist, variant) {

  mats <- vapply(datalist, is.matrix, logical(1))
  nonas <- lapply(datalist[mats], function(x) !is.na(x))

  # AMHG has-data matrix (needs slopes and widths)
  hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
  hasdat_w <- (!is.na(datalist[["Wobs"]])) * 1
  hasdat_amhg <- hasdat_w * hasdat_s

  # Replace NA's with zeros so Stan will accept the data
  datalist[["Wobs"]][!hasdat_amhg] <- 0
  datalist[["Sobs"]][!hasdat_amhg] <- 0

  # Manning has-data matrix (only nonzero if all Manning obs present)
  if (identical(setdiff(c("Wobs", "Sobs", "dAobs"),
                        names(datalist[mats])),
                character(0))) {
    #hasdat_s <- (!is.na(datalist[["Sobs"]])) * 1
    hasdat_a <- (!is.na(datalist[["dAobs"]])) * 1

    hasdat_man <- hasdat_amhg * hasdat_a

    # Replace NA's with zeros so Stan will accept the data
    datalist[["dAobs"]][!hasdat_man] <- 0

  } else {
    hasdat_man <- matrix(0, nrow = nrow(hasdat_amhg), ncol = ncol(hasdat_amhg))
  }

  if (!is.null(datalist[["dAobs"]])) {
    dA_shift <- apply(datalist[["dAobs"]], 1, function(x) median(x) - min(x))
  } else {
    dA_shift <- rep(0, nrow(datalist[["Wobs"]]))
  }

  newbits <- list(
    hasdat_man = hasdat_man,
    hasdat_amhg = hasdat_amhg,
    ntot_man = sum(hasdat_man),
    ntot_amhg = sum(hasdat_amhg),
    dA_shift = dA_shift
  )

  out <- c(datalist, newbits)
  out
}

#' Establish prior hyperparameters for BAM estimation
#'
#' Produces a bampriors object that can be passed to bam_estimate function
#'
#' @useDynLib geoBAMr, .registration = TRUE
#' @param bamdata An object of class bamdata, as returned by \code{bam_data}
#' @param variant Which geoBAM variant to use. Options are "manning_amhg" (default),
#'   "manning", or "amhg".
#' @param classification Which classification framework to use. Options are 'expert' (default),
#'   or 'unsupervised'.
#' @param ... Optional manually set parameters. Unquoted expressions are allowed,
#'   e.g. \code{logQ_sd = cv2sigma(0.8)}. Additionally, any variables present in
#'   \code{bamdata} may be referenced, e.g. \code{lowerbound_logQ = log(mean(Wobs)) + log(5)}
#' @export

bam_priors <- function(bamdata,
                       variant = c("manning_amhg", "manning", "amhg"),
                       classification = c('expert', 'unsupervised'),
                       ...) {
  variant <- match.arg(variant)
  classification <- match.arg(classification)
  if (variant != "amhg" && !attr(bamdata, "manning_ready"))
    stop("bamdata must have dA data for non-amhg variants.")

  force(bamdata)
  paramset <- bam_settings("paramnames")

  myparams0 <- rlang::quos(..., .named = TRUE)
  myparams <- do.call(settings::clone_and_merge,
                      args = c(list(options = bam_settings), myparams0))

  quoparams <- myparams()[-1] # first one is parameter set
  params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)

  if (classification == 'unsupervised'){ #recalculate parameter set if classification switched to unsupervised
    paramset <- bam_settings_unsupervised("paramnames")

    myparams0 <- rlang::quos(..., .named = TRUE)
    myparams <- do.call(settings::clone_and_merge,
                        args = c(list(options = bam_settings_unsupervised), myparams0))

    quoparams <- myparams()[-1] # first one is parameter set
    params <- lapply(quoparams, rlang::eval_tidy, data = bamdata)
  }

  if (!length(params[["logQ_sd"]]) == bamdata$nt)
    params$logQ_sd <- rep(params$logQ_sd, length.out = bamdata$nt)

  if (!identical(dim(params[["sigma_man"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_man <- matrix(rep(params$sigma_man,
                                   length.out = bamdata$nt * bamdata$nx),
                               nrow = bamdata$nx)
  }

  if (!identical(dim(params[["sigma_amhg"]]),
                 as.integer(c(bamdata$nx, bamdata$nt)))) {
    params$sigma_amhg <- matrix(rep(params$sigma_amhg,
                                    length.out = bamdata$nt * bamdata$nx),
                                nrow = bamdata$nx)
  }

  #just priors the user wants to see
  user_paramset <- c('lowerbound_logQ', 'upperbound_logQ', 'lowerbound_logWc', 'upperbound_logWc', 'lowerbound_logQc', 'upperbound_logQc',
                     'logWc_hat',"logQc_hat",
                     'logQ_sd','logWc_sd','logQc_sd',
                     "Werr_sd", "Serr_sd", "dAerr_sd",
                     "sigma_man", "sigma_amhg")
  user_paramset <- params[user_paramset]

  #total priors needed to run geoBAM
  geoBAM_paramset <- c( "lowerbound_A0", "upperbound_A0", "lowerbound_logn", "upperbound_logn","lowerbound_b","upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb','upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                        "logA0_hat", "logn_hat","b_hat",  'logWb_hat', 'logDb_hat', 'logr_hat',
                       "logA0_sd", "logn_sd", "b_sd",'logWb_sd', 'logDb_sd', 'logr_sd')
  geoBAMparams <- params[geoBAM_paramset]

  riverType <- params[["River_Type"]]

  out <- list( 'River_Type'=riverType, 'river_type_priors'=geoBAMparams, 'other_priors'=user_paramset)
  out <- structure(out,
                   class = c("bampriors"))
  out
}

compose_bam_inputs <- function(bamdata, priors = bam_priors(bamdata)) {

  inps <- c(bamdata, priors)

  out <- inps
  out

}


#' Take a random sample of a bamdata object's cross-sections.
#'
#' @param bamdata a bamdata object, as returned by \code{bam_data()}
#' @param n Number of cross-sections to
#' @param seed option RNG seed, for reproducibility.
#' @importFrom methods is
#' @export
sample_xs <- function(bamdata, n, seed = NULL) {

  stopifnot(is(bamdata, "bamdata"))

  if (n >= bamdata$nx)
    return(bamdata)

  if (!is.null(seed))
    set.seed(seed)
  keepxs <- sort(sample(1:bamdata$nx, size = n, replace = FALSE))

  bamdata$nx <- n
  bamdata$Wobs <- bamdata$Wobs[keepxs, ]

  if (!is.null(bamdata$Sobs)) {
    bamdata$Sobs <- bamdata$Sobs[keepxs, ]
    bamdata$dAobs <- bamdata$dAobs[keepxs, ]
  }

  bamdata
}



#' Calculate lognormal moments based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
#' @importFrom stats dnorm pnorm
ln_moms <- function(obs, err_sigma, a = 0) {
  alpha <- (a - obs) / err_sigma
  Z <- 1 - pnorm(alpha)

  mean <- obs + (dnorm(alpha)) / Z * err_sigma
  sdquan <- 1 + (alpha * dnorm(alpha)) / Z -
    (dnorm(alpha) / Z)^2
  sd <- err_sigma * sqrt(sdquan)

  out <- list(mean = mean, sd = sd)
  out
}

#' Calculate lognormal sigma parameter based on truncated normal parameters
#'
#' Used to put measurement errors into original log-normal parameterization.
#'
#' @param obs A numeric vector of observations
#' @param err_sigma Standard deviation of measurement error
#' @param a zero-reference for method of moments.
ln_sigsq <- function(obs, err_sigma, a = 0) {
  moms <- ln_moms(obs = obs, err_sigma = err_sigma, a = a)
  mn <- unname(moms[["mean"]])
  sd <- unname(moms[["sd"]])
  mu <- 2 * log(mn) - 0.5 * log(sd^2 + mn^2)
  sigsq <- 2 * log(mn) - 2 * mu

  sigsq
}
