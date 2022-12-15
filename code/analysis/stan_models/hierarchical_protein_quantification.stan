data { 
    int<lower=1> J; // Number of unique carbon sources
    int<lower=1> K; // Number of calibration curve replicates
    int<lower=1> N_meas; // Total number of measurements 
    int<lower=1> N_calib; // Total number of calibration measurements
    array[N_meas] int<lower=1, upper=J> meas_idx; // ID vector mapping protein measurements to carbon source 
    array[N_calib] int<lower=1, upper=K> calib_idx; // ID vector for mapping calibration curves

    // Measured data for calibration
    vector<lower=0>[N_calib] concentration;
    vector<lower=0>[N_calib] od_595nm_calib;

    // Measured data for quantification
    vector<lower=0>[N_meas] od_595nm_meas;
    vector<lower=0>[N_meas] od_600nm_meas;
    vector<lower=0>[N_meas] dil_factor;
    vector<lower=0>[N_meas] ext_volume;
    vector<lower=0>[N_meas] cult_volume;
    }

transformed data {
    vector<lower=0>[N_meas] od_per_biomass = od_595nm_meas .* dil_factor .* ext_volume ./ (cult_volume .* od_600nm_meas);
}

parameters {
    // Hyperparameters for calibration
    real<lower=0> slope;
    real<lower=0> intercept;
    real tau;

    // Mid-level (uncentered) parameters for calibration
    vector[K] slope_1_tilde;
    vector[K] intercept_1_tilde;
    real<lower=0> calib_sigma;

    // Parameters for protein measurements 
    vector[J] od_per_biomass_mu;
    vector<lower=0>[J] od_per_biomass_sigma;
}

transformed parameters { 
    // Uncenter the hyperparameters
    // real slope = exp(log_slope);
    // real intercept = exp(log_intercept);
    vector<lower=0>[K] slope_1 = slope_1_tilde + tau * slope;
    vector<lower=0>[K] intercept_1 = intercept_1_tilde + tau * intercept;
}

model {
    // Priors
    slope ~ gamma(5, 5);
    intercept ~ gamma(1.5, 3.5);
    slope_1_tilde ~ std_normal();
    intercept_1_tilde ~ std_normal();
    tau ~ std_normal();
    calib_sigma ~ std_normal();
    od_per_biomass_mu ~ normal(1000, 500);
    od_per_biomass_sigma ~ normal(0, 100);

    // Likelihoods
    od_595nm_calib ~ normal(slope_1[calib_idx] .* concentration[calib_idx] + intercept_1[calib_idx], calib_sigma);
    od_per_biomass ~ normal(od_per_biomass_mu[meas_idx], od_per_biomass_sigma[meas_idx]);
}

generated quantities { 
    vector[N_meas] od_per_biomass_rep;
    vector[N_calib] od_595nm_calib_rep;
    vector[J] prot_per_biomass;
    // Posterior Predictive Checks
    for (i in 1:N_meas) {
        od_per_biomass_rep[i] = normal_rng(od_per_biomass_mu[meas_idx[i]], od_per_biomass_sigma[meas_idx[i]]);
    }
    for (i in 1:N_calib) { 
        od_595nm_calib_rep[i] = normal_rng(slope_1[calib_idx[i]] * concentration[i] + intercept_1[calib_idx[i]], calib_sigma);
    } 
    for (i in 1:J) {
        prot_per_biomass[i] = (od_per_biomass_mu[i] - intercept) ./ slope;
    }
}