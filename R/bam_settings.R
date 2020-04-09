#Bounds on Manning's n are manually set to log(0.01) and log(0.05).
  #Functions are written to revert back to standard approach to defining priors, currently commented out.
#Mark's defualt BAM priors are also commented out.



#' Options manager for geoBAMr defaults
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
bam_settings <- settings::options_manager(
  paramnames = c("lowerbound_logQ", "upperbound_logQ", "lowerbound_A0",
                 "upperbound_A0", "lowerbound_logn", "upperbound_logn", "lowerbound_logQc",
                 "upperbound_logQc", "lowerbound_logWc", "upperbound_logWc", "lowerbound_b",
                 "upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb',
                 'upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                 "sigma_man", "sigma_amhg",
                 "logQc_hat", "logWc_hat", "b_hat", "logA0_hat", "logn_hat", 'logWb_hat', 'logDb_hat', 'logr_hat',
                 "logQ_sd", "logQc_sd", "logWc_sd", "b_sd", "logA0_sd", "logn_sd", 'logWb_sd', 'logDb_sd', 'logr_sd',
                 "Werr_sd", "Serr_sd", "dAerr_sd",
                 'river_type'),

  # Bounds on parameters
  lowerbound_logQ = rlang::quo(maxmin(log(Wobs)) + log(0.5) + log(0.5)),
  upperbound_logQ = rlang::quo(minmax(log(Wobs)) + log(40) + log(5)),

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0(Wobs)), #0.72,
  upperbound_A0 = rlang::quo(estimate_upperboundA0(Wobs)), #114500,
  lowerbound_logn = log(0.01), #rlang::quo(estimate_lowerboundlogn(Wobs)), #-4.60517,
  upperbound_logn = log(0.05), #rlang::quo(estimate_upperboundlogn(Wobs)), #-2.995732,

  lowerbound_logQc = 0.01,
  upperbound_logQc = 10,
  lowerbound_logWc = 1,
  upperbound_logWc = 8, # 3 km
  lowerbound_b = rlang::quo(estimate_lowerboundb(Wobs)), #0.000182,
  upperbound_b = rlang::quo(estimate_upperboundb(Wobs)), #0.773758,

  lowerbound_logDb = rlang::quo(estimate_lowerboundlogDb(Wobs)), #-3.02002,
  upperbound_logDb = rlang::quo(estimate_upperboundlogDb(Wobs)), #3.309359,
  lowerbound_logWb = rlang::quo(estimate_lowerboundlogWb(Wobs)),#-0.12273,
  upperbound_logWb = rlang::quo(estimate_upperboundlogWb(Wobs)), #7.006786,
  lowerbound_logr = rlang::quo(estimate_lowerboundlogr(Wobs)), #-2.58047,
  upperbound_logr = rlang::quo(estimate_upperboundlogr(Wobs)), #8.03772,


  # *Known* likelihood parameters
  sigma_man = 0.25,
  sigma_amhg = 0.22, # UPDATE THIS FROM CAITLINE'S WORK


  # Hyperparameters
  logQc_hat = rlang::quo(mean(logQ_hat)),
  logWc_hat = rlang::quo(mean(log(Wobs))),
  b_hat = rlang::quo(estimate_b(Wobs)),
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logn_hat = rlang::quo(estimate_logn(Sobs, Wobs)),
  logWb_hat = rlang::quo(estimate_logWb(Wobs)),
  logDb_hat = rlang::quo(estimate_logDb(Wobs)),
  logr_hat = rlang::quo(estimate_logr(Wobs)),

  logQ_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1
  logQc_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1; UPDATE THIS
  logWc_sd = sqrt(log(0.01)^2 + 1),

  #set from my model outputs
  b_sd = rlang::quo(estimate_bSD(Wobs)), #0.068077044,
  logA0_sd = rlang::quo(estimate_A0SD(Wobs)), #0.58987527,
  logn_sd = rlang::quo(estimate_lognSD(Wobs)), #0.761673112,
  logWb_sd = rlang::quo(estimate_logWbSD(Wobs)), #0.137381044,
  logDb_sd = rlang::quo(estimate_logDbSD(Wobs)), #0.576212733,
  logr_sd = rlang::quo(estimate_logrSD(Wobs)), #0.67332688,

  # Observation errors.
  Werr_sd = 10,
  Serr_sd = 1e-5,
  dAerr_sd = 10,

  #Classified river type
  river_type = rlang::quo(apply(Wobs, 1, classify_func_unsupervised))
)


#' Options manager for geoBAMr defaults using unsupervised classification
#'
#' @param ... (Optional) named settings to query or set.
#' @param .__defaults See \code{?settings::option_manager}
#' @param .__reset See \code{?settings::option_manager}
#' @export
bam_settings_unsupervised <- settings::options_manager(
  paramnames = c("lowerbound_logQ", "upperbound_logQ", "lowerbound_A0",
                 "upperbound_A0", "lowerbound_logn", "upperbound_logn", "lowerbound_logQc",
                 "upperbound_logQc", "lowerbound_logWc", "upperbound_logWc", "lowerbound_b",
                 "upperbound_b", "lowerbound_logWb", 'upperbound_logWb', 'lowerbound_logDb',
                 'upperbound_logDb', 'lowerbound_logr', 'upperbound_logr',
                 "sigma_man", "sigma_amhg",
                 "logQc_hat", "logWc_hat", "b_hat", "logA0_hat", "logn_hat", 'logWb_hat', 'logDb_hat', 'logr_hat',
                 "logQ_sd", "logQc_sd", "logWc_sd", "b_sd", "logA0_sd", "logn_sd", 'logWb_sd', 'logDb_sd', 'logr_sd',
                 "Werr_sd", "Serr_sd", "dAerr_sd",
                 'river_type'),

  # Bounds on parameters
  lowerbound_logQ = rlang::quo(maxmin(log(Wobs)) + log(0.5) + log(0.5)),
  upperbound_logQ = rlang::quo(minmax(log(Wobs)) + log(40) + log(5)),

  lowerbound_A0 = rlang::quo(estimate_lowerboundA0_unsupervised(Wobs)), #0.72,
  upperbound_A0 = rlang::quo(estimate_upperboundA0_unsupervised(Wobs)), #114500,
  lowerbound_logn = log(0.05), #rlang::quo(estimate_lowerboundlogn_unsupervised(Wobs)), #-4.60517,
  upperbound_logn = log(0.01), #rlang::quo(estimate_upperboundlogn_unsupervised(Wobs)), #-2.995732,

  lowerbound_logQc = 0.01,
  upperbound_logQc = 10,
  lowerbound_logWc = 1,
  upperbound_logWc = 8, # 3 km
  lowerbound_b = rlang::quo(estimate_lowerboundb_unsupervised(Wobs)), #0.000182,
  upperbound_b = rlang::quo(estimate_upperboundb_unsupervised(Wobs)), #0.773758,

  lowerbound_logDb = rlang::quo(estimate_lowerboundlogDb_unsupervised(Wobs)), #-3.02002,
  upperbound_logDb = rlang::quo(estimate_upperboundlogDb_unsupervised(Wobs)), #3.309359,
  lowerbound_logWb = rlang::quo(estimate_lowerboundlogWb_unsupervised(Wobs)),#-0.12273,
  upperbound_logWb = rlang::quo(estimate_upperboundlogWb_unsupervised(Wobs)), #7.006786,
  lowerbound_logr = rlang::quo(estimate_lowerboundlogr_unsupervised(Wobs)), #-2.58047,
  upperbound_logr = rlang::quo(estimate_upperboundlogr_unsupervised(Wobs)), #8.03772,


  # *Known* likelihood parameters
  sigma_man = 0.25,
  sigma_amhg = 0.22, # UPDATE THIS FROM CAITLINE'S WORK


  # Hyperparameters
  logQc_hat = rlang::quo(mean(logQ_hat)),
  logWc_hat = rlang::quo(mean(log(Wobs))),
  b_hat = rlang::quo(estimate_b_unsupervised(Wobs)),
  logA0_hat = rlang::quo(estimate_logA0_unsupervised(Wobs)),
  logn_hat = rlang::quo(estimate_logn_unsupervised(Sobs, Wobs)),
  logWb_hat = rlang::quo(estimate_logWb_unsupervised(Wobs)),
  logDb_hat = rlang::quo(estimate_logDb_unsupervised(Wobs)),
  logr_hat = rlang::quo(estimate_logr_unsupervised(Wobs)),

  logQ_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1
  logQc_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1; UPDATE THIS
  logWc_sd = sqrt(log(0.01)^2 + 1),

  #set from my model outputs
  b_sd = rlang::quo(estimate_bSD_unsupervised(Wobs)), #0.068077044,
  logA0_sd = rlang::quo(estimate_A0SD_unsupervised(Wobs)), #0.58987527,
  logn_sd = rlang::quo(estimate_lognSD_unsupervised(Wobs)), #0.761673112,
  logWb_sd = rlang::quo(estimate_logWbSD_unsupervised(Wobs)), #0.137381044,
  logDb_sd = rlang::quo(estimate_logDbSD_unsupervised(Wobs)), #0.576212733,
  logr_sd = rlang::quo(estimate_logrSD_unsupervised(Wobs)), #0.67332688,

  # Observation errors.
  Werr_sd = 10,
  Serr_sd = 1e-5,
  dAerr_sd = 10,

  #Classified river type
  river_type = rlang::quo(apply(Wobs, 1, classify_func_unsupervised))
)
