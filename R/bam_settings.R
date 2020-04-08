
#' Options manager for BAM defaults
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
     "Werr_sd", "Serr_sd", "dAerr_sd"),

  # Bounds on parameters
  lowerbound_logQ = rlang::quo(maxmin(log(Wobs)) + log(0.5) + log(0.5)),
  upperbound_logQ = rlang::quo(minmax(log(Wobs)) + log(40) + log(5)),

  lowerbound_A0 = 0.72,
  upperbound_A0 = 114500,
  lowerbound_logn = -4.61,
  upperbound_logn = -2.995732,

  lowerbound_logQc = 0,
  upperbound_logQc = 10,
  lowerbound_logWc = 1,
  upperbound_logWc = 8, # 3 km
  lowerbound_b = 0.000182,
  upperbound_b = 0.773758,

  lowerbound_logDb = -3.02002,
  upperbound_logDb = 3.309359,
  lowerbound_logWb = -0.12273,
  upperbound_logWb = 7.006786,
  lowerbound_logr = -2.58047,
  upperbound_logr = 8.03772,


  # *Known* likelihood parameters
  sigma_man = 0.25,
  sigma_amhg = 0.22, # UPDATE THIS FROM CAITLINE'S WORK


  # Hyperparameters
  logQc_hat = rlang::quo(mean(logQ_hat)),
  logWc_hat = rlang::quo(mean(log(Wobs))),
  b_hat = rlang::quo(estimate_b(Wobs)),
  logA0_hat = rlang::quo(estimate_logA0(Wobs)),
  logn_hat = rlang::quo(estimate_logn(Sobs)),
  logWb_hat = rlang::quo(estimate_logWb(Wobs)),
  logDb_hat = rlang::quo(estimate_logDb(Wobs)),
  logr_hat = rlang::quo(estimate_logr(Wobs)),

  logQ_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1
  logQc_sd = sqrt(log(1^2 + 1)), # CV of Q equals 1; UPDATE THIS
  logWc_sd = sqrt(log(0.01)^2 + 1),

  #set from my model outputs
  b_sd = 0.068,
  logA0_sd = 0.590,
  logn_sd = 0.762,
  logWb_sd = 0.137,
  logDb_sd = 0.576,
  logr_sd = 0.673,

  # Observation errors.
  Werr_sd = 10,
  Serr_sd = 1e-5,
  dAerr_sd = 10
)
