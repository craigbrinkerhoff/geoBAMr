//Craig BAM stan model
//Interventions from Mark's model:
  //1) Uses AMHG flow law from https://10.1029/2019GL084529, reformatted to derive a likelihood function for Bayesian inference.
  //2) Manning's n prior is space-varying
  //3) '_sd' priors are space varying

//line 227 must be set to inc_a to run amhg and to inc_m to run mannings.... I need someway to switch this

functions {
  // Conversion from array to vector takes place by row.
  // Nested elements are appended in sequence.


  // Convert an array to a vector based on a binary matrix
  // indicating non-missing data

  vector ragged_vec(vector[] x, int[,] bin) {
    vector[num_elements(x)] out;
    int ind;

    ind = 1;
    for (i in 1:size(x)) {
      for (t in 1:num_elements(x[1])) {
        if (bin[i, t] == 1) {
          out[ind] = x[i, t];
          ind += 1;
        }
      }
    }
    // print(out);
    return(out[1:(ind - 1)]);
  }

  // Repeat elements of a "row" vector to match with 2-D array vectorization
  vector ragged_row(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[t];
        }
      }
    }
    return(out[1:ind]);
  }

  // Repeat elements of a "column" vector to match with 2-D array vectorization
  vector ragged_col(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[i];
        }
      }
    }
    return(out[1:ind]);
  }

  // indices of vectorized bin1 that are also in vectorized bin2
  int[] commoninds(int[,] bin1, int[,] bin2) {
    int out[num_elements(bin1)];
    int vecinds[size(bin1), num_elements(bin1[1])];
    int ctr;
    int ind;

    ctr = 0;
    for (i in 1:size(bin1)) {
      for (t in 1:num_elements(bin1[1])) {
        if (bin1[i, t] == 1) {
          ctr += 1;
          vecinds[i, t] = ctr;
        }
      }
    }

    ind = 0;
    for (i in 1:size(vecinds)) {
      for (t in 1:size(vecinds[1])) {
        if (bin2[i, t] == 1) {
          ind += 1;
          out[ind] = vecinds[i, t];
        }
      }
    }
    return(out[1:ind]);
  }
}

data {

  // Options
  int<lower=0, upper=1> inc_m; // Include Manning? 0=no; 1=yes;
  int<lower=0, upper=1> inc_a; // Include AMHG? 0=no; 1=yes;
  int<lower=0, upper=1> meas_err; // 0=no; 1=yes;


  // Dimensions
  int<lower=0> nx; // number of reaches
  int<lower=0> nt; // number of times
  int<lower=0> ntot_man; // total number of non-missing Manning observations
  int<lower=0> ntot_amhg; // number of non-missing width observations

  // // Missing data
  int<lower=0,upper=1> hasdat_man[nx, nt]; // matrix of 0 (missing), 1 (not missing)
  int<lower=0,upper=1> hasdat_amhg[nx, nt];

  // *Actual* data
  vector[nt] Wobs[nx]; // measured widths, including placeholders for missing
  vector[nt] Sobs[nx]; // measured slopes
  vector[nt] dAobs[nx]; // measured partial area
  vector[nx] dA_shift; // adjustment from min to median


  real<lower=0> Werr_sd;
  real<lower=0> Serr_sd;
  real<lower=0> dAerr_sd;

  // Hard bounds on parameters
  real lowerbound_logQ;
  real upperbound_logQ;

  real lowerbound_A0; // These must be scalars, unfortunately.
  real upperbound_A0;
  real lowerbound_logn;
  real upperbound_logn;

  real lowerbound_logQc;
  real upperbound_logQc;
  real lowerbound_logWc;
  real upperbound_logWc;
  real lowerbound_b;
  real upperbound_b;

  real lowerbound_logWb;
  real upperbound_logWb;
  real lowerbound_logDb;
  real upperbound_logDb;
  real lowerbound_logr;
  real upperbound_logr;

  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_man[nx]; // Manning error standard deviation
  vector<lower=0>[nt] sigma_amhg[nx]; // AMHG error standard deviation


  // Hyperparameters
  vector[nt] logQ_hat; // prior mean on logQ
  real logQc_hat; // prior mean on logQc
  real logWc_hat;
  real b_hat[nx]; // ADD CHECK ON THIS FOR DATA PREP
  real logA0_hat[nx];
  real logn_hat[nx];
  real logWb_hat[nx];
  real logDb_hat[nx];
  real logr_hat[nx];

  vector<lower=0>[nt] logQ_sd;
  real<lower=0> logQc_sd;
  real<lower=0> logWc_sd;
  real<lower=0> b_sd[nx];
  real<lower=0> logA0_sd[nx];
  real<lower=0> logn_sd[nx];
  real<lower=0> logWb_sd[nx];
  real<lower=0> logDb_sd[nx];
  real<lower=0> logr_sd[nx];
}



transformed data {
  // Transformed data are *vectors*, not arrays. This to allow ragged structure

  vector[nt] dApos_array[nx];

  vector[ntot_man] Wobsvec_man;
  vector[ntot_amhg] Wobsvec_amhg;
  vector[inc_a ? ntot_amhg : ntot_man] Wobsvec;
  vector[ntot_man] Sobsvec_man;
  vector[ntot_amhg] Sobsvec_amhg;
  vector[inc_a ? ntot_amhg : ntot_man] Sobsvec;

  vector[ntot_man] logWobs_man;
  vector[ntot_amhg] logWobs_amhg;
  vector[ntot_man] logSobs_man;
  vector[ntot_amhg] logSobs_amhg;
  vector[ntot_man] dApos_obs;
  vector[ntot_man] sigmavec_man;
  vector[ntot_amhg] sigmavec_amhg;

  int maninds_amhg[ntot_man];

  int ntot_w; // how many widths in likelihood: equal to ntot_man unless inc_a
  ntot_w = inc_a ? ntot_amhg : ntot_man;

  for (i in 1:nx) {
    dApos_array[i] = dAobs[i] - min(dAobs[i]); // make all dA positive
  }

  // convert pseudo-ragged arrays to vectors
  Wobsvec_man = ragged_vec(Wobs, hasdat_man);
  Wobsvec_amhg = ragged_vec(Wobs, hasdat_amhg);
  Wobsvec = inc_a ? Wobsvec_amhg : Wobsvec_man;
  Sobsvec_man = ragged_vec(Sobs, hasdat_man);
  Sobsvec_amhg = ragged_vec(Sobs, hasdat_amhg);
  dApos_obs = ragged_vec(dApos_array, hasdat_man);

  logWobs_man = log(Wobsvec_man);
  logSobs_man = log(Sobsvec_man);
  logSobs_amhg = log(Sobsvec_amhg);

  sigmavec_man = ragged_vec(sigma_man, hasdat_man);
  sigmavec_amhg = ragged_vec(sigma_amhg, hasdat_amhg);

  maninds_amhg = commoninds(hasdat_amhg, hasdat_man);
}

parameters {
  vector<lower=lowerbound_logn,upper=upperbound_logn>[nx] logn[1]; //for reach-defined n
  vector<lower=lowerbound_logQ,upper=upperbound_logQ>[nt] logQ;
  vector<lower=lowerbound_A0,upper=upperbound_A0>[nx] A0[inc_m];

  real<lower=lowerbound_logWc,upper=upperbound_logWc> logWc[inc_a];
  real<lower=lowerbound_logQc,upper=upperbound_logQc> logQc[inc_a];
  vector<lower=lowerbound_b,upper=upperbound_b>[nx] b[inc_a];

  vector<lower=lowerbound_logWb, upper=upperbound_logWb>[nx] logWb[inc_a];
  vector<lower=lowerbound_logDb, upper=upperbound_logDb>[nx] logDb[inc_a];
  vector<lower=lowerbound_logr, upper=upperbound_logr>[nx] logr[inc_a];

  vector<lower=0>[ntot_w] Wact[meas_err];
  vector<lower=0>[ntot_man] Sact[meas_err * inc_m];
  vector[ntot_man] dApos_act[meas_err * inc_m];
}


transformed parameters {

  vector[ntot_man] man_lhs[inc_m]; // LHS for Manning likelihood
  vector[ntot_man] logA_man[inc_m]; // log area for Manning's equation
  vector[ntot_man] man_rhs[inc_m]; // RHS for Manning likelihood
  vector[ntot_man] Wact_man[inc_m * meas_err]; // subset of Wact parameter
  vector[ntot_man] logQ_man[inc_m]; // location-repeated logQ

  vector[ntot_amhg] amhg_rhs[inc_a]; // RHS for AMHG likelihood
  vector[ntot_amhg] logQ_amhg[inc_a]; // location-repeated logQ
  vector[ntot_amhg] logQc_amhg[inc_a]; //new Qc term

  // Manning params
  if (inc_m) {
    if (meas_err) {
      Wact_man[1] = Wact[1][maninds_amhg];
      logA_man[1] = log(ragged_col(A0[1], hasdat_man) + dApos_act[1]);
      man_lhs[1] = 4. * log(Wact_man[1]) - 3. * log(Sact[1]);
    }
    else{
      logA_man[1] = log(ragged_col(A0[1], hasdat_man) + dApos_obs);
      man_lhs[1] = 4. * logWobs_man - 3. * logSobs_man;
    }

    logQ_man[1] = ragged_row(logQ, hasdat_man);
    man_rhs[1] = 10. * logA_man[1] - 6. * ragged_col(logn[1], hasdat_man) - 6. * logQ_man[1];
  }

  if (inc_a) {
  // AMHG params
    logQ_amhg[1] = ragged_row(logQ, hasdat_amhg);

    //new AMHG likelihood function
    logQc_amhg[1] = (1 ./ ragged_col(b[1], hasdat_amhg)) .* (logWc[1]- (ragged_col(b[1], hasdat_amhg) .*
                                           (((-1.67) * ragged_col(logDb[1], hasdat_amhg)) +
	                                         ((-1.67) * ragged_col(logr[1], hasdat_amhg)) -
	                                         ((-1.67) * (ragged_col(logr[1], hasdat_amhg)+1)) +
	                                         ((1.67 * ragged_col(logr[1], hasdat_amhg)) .* ragged_col(logWb[1], hasdat_amhg)) +
	                                         (ragged_col(logn[1], hasdat_amhg)) +
	                                         (-(0.5)*logSobs_amhg[1]))));

    amhg_rhs[1] = ragged_col(b[1], hasdat_amhg) .* (logQ_amhg[1] - ragged_col(logQc_amhg[1], hasdat_amhg)) + logWc[1];
  }
}

model {

  // Priors
  logQ ~ normal(logQ_hat, logQ_sd);

  if (inc_m) {
    A0[1] + dA_shift[1] ~ lognormal(logA0_hat, logA0_sd);
    logn[1] ~ normal(logn_hat, logn_sd);
  }
  if (inc_a) {
    b[1] ~ normal(b_hat, b_sd);
    logWc ~ normal(logWc_hat, logWc_sd);
    logQc ~ normal(logQc_hat, logQc_sd);

    logn[1] ~ normal(logn_hat, logn_sd);

	  logDb[1] ~ normal(logDb_hat, logDb_sd);
	  logr[1] ~ normal(logr_hat, logr_sd);
	  logWb[1] ~ normal(logWb_hat, logWb_sd);
  }
  // Likelihood and observation error

  // Manning likelihood
  if (inc_m) {
    man_lhs[1] ~ normal(man_rhs[1], 6 * sigmavec_man);
  }

  // Latent vars for measurement error
  if (meas_err) {
    Wact[1] ~ normal(Wobsvec, Werr_sd); // W meas err

    if (inc_m) {
      Sact[1] ~ normal(Sobsvec_man, Serr_sd); // S meas err
      dApos_act[1] ~ normal(dApos_obs, dAerr_sd); // dA meas err
      target += -log(Wact[1]); // Jacobian adjustments
      target += -log(Sact[1]);
    }
    if (inc_a) { // AMHG likelihood (w/ meas err)
      Wact[1] ~ lognormal(amhg_rhs[1], sigmavec_amhg);
    }
  }

  else {
    if (inc_a) { // AMHG likelihood (w/o meas err)
      Wobsvec ~ lognormal(amhg_rhs[1], sigmavec_amhg);
    }
  }
}
