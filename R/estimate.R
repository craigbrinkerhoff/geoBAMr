
#' Estimate BAM 
#' 
#' Fits a BAM model of one of several variants using Hamiltonian Monte Carlo.
#' 
#' @param bamdata A bamdata object, as produced by \code{bam_data()}
#' @param variant Which BAM variant to use: amhg, manning_amhg, or manning
#' @param meas_error Include measurement error in inference? Setting this to TRUE 
#'   will slow down the inference by roughly an order of mangnitude.
#' @param reparam Reparameterize measurement errors to speed up sampling?
#' @param bampriors A bampriors object. If none is supplied, defaults are used 
#'   from calling \code{bam_priors(bamdata)} (with no other arguments).
#' @param cores Number of processing cores for running chains in parallel. 
#'   See \code{?rstan::sampling}. Defaults to \code{parallel::detectCores()}.
#' @param chains A positive integer specifying the number of Markov chains. 
#'   The default is 3.
#' @param iter Number of iterations per chain (including warmup). Defaults to 1000.
#' @param stanmodel A \code{stanmodel} object to use instead of one of the default 
#'   models. 
#' @param pars (passed to \code{rstan::sampling()}) A vector of character strings specifying 
#'   parameters of interest to be returned in the stanfit object. If not specified, 
#'   a default parameter set is returned.
#' @param include (passed to \code{rstan::sampling()}) Defaults to FALSE, which 
#'   excludes parameters specified in \code{pars} from the returned model.
#' @param ... Other arguments passed to rstan::sampling() for customizing the 
#'   Monte Carlo sampler
#' @import rstan
#' @export

bam_estimate <- function(bamdata, 
                         variant = c("manning", "amhg", "manning_amhg"), 
                         bampriors = NULL, 
                         meas_error = TRUE,
                         reparam = TRUE,
                         cores = getOption("mc.cores", default = parallel::detectCores()),
                         chains = 3L,
                         iter = 1000L,
                         stanmodel = NULL,
                         pars = NULL, 
                         include = FALSE,
                         ...) {
  variant <- match.arg(variant)
  stopifnot(is(bamdata, "bamdata"))
  if (is.null(bampriors))
    bampriors <- bam_priors(bamdata, variant = variant)
  stopifnot(is(bampriors, "bampriors"))
  
  baminputs <- compose_bam_inputs(bamdata, bampriors)
  
  if (!is.null(stanmodel)) {
    stopifnot(inherits(stanmodel, "stanmodel"))
    stanfit <- stanmodel
  } else {
    stanfit <- stanmodels[["master"]]
  }
  baminputs$meas_err <- ifelse(meas_error && !reparam, 1, 0)
  baminputs$inc_m <- ifelse(variant %in% c("manning", "manning_amhg"), 1, 0)
  baminputs$inc_a <- ifelse(variant %in% c("amhg", "manning_amhg"), 1, 0)
  
  if (is.null(pars)) {
    pars <- c("man_rhs", "amhg_rhs", "logWSpart", 
              "logQtn", "logQnbar",
              "Sact", "Wact", "dAact")
  }
  
  if (reparam && meas_error) {
    logS_sigsq_obs <- ln_sigsq(obs = baminputs$Sobs, err_sigma = baminputs$Serr_sd)
    logW_sigsq_obs <- ln_sigsq(obs = baminputs$Wobs, err_sigma = baminputs$Werr_sd)
    baminputs$sigma_man <- sqrt(baminputs$sigma_man^2 + 
                                  logS_sigsq_obs * (3 / 6)^2 +
                                  logW_sigsq_obs * (4 / 6)^2)
    baminputs$sigma_amhg <- sqrt(baminputs$sigma_amhg^2 +
                                   logW_sigsq_obs)
  }
  
  out <- sampling(stanfit, data = baminputs, 
                  cores = cores, chains = chains, iter = iter,  
                  pars = pars, include = include,
                  ...)
  
  out
}


