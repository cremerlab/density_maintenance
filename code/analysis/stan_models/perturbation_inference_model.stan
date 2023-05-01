data { 
    // Measurement dimensions
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_cal;
    int<lower=1> N_brad;
    int<lower=1> N_flow;
    int<lower=1> N_growth;

    // Condition dimensions
    int<lower=1> J_size_cond;
    int<lower=1> J_brad_cond; 
    int<lower=1> J_growth_cond;
    int<lower=1> J_flow_cond;

    // Size measurements
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] length;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] peri_volume;

    // Growth measurements
    vector<lower=0>[N_growth] growth_rates;

    // Total protein measurements
    vector<lower=0>[N_prot] total_protein;
    vector<lower=0>[N_prot] total_protein_lam;

    // Flow measurements
    vector<lower=0>[N_flow] flow_events;
    array[N_flow] int<lower=1, upper=J_flow_cond> flow_idx;
    array[N_flow] int<lower=1, upper=J_size_cond> flow_mapper;

    // Bradford assay measurements
    vector<lower=0>[N_brad] brad_od595; 
    vector<lower=0>[N_brad] brad_od600; 
    vector<lower=0>[N_brad] conv_factor;
    vector<lower=0>[N_cal] concentration; 
    vector<lower=0>[N_cal] cal_od;

    array[N_growth] int<lower=1, upper=J_growth_cond> growth_idx;  
    array[N_size] int<lower=1, upper=J_size_cond> size_idx;  
    array[N_brad] int<lower=1, upper=J_brad_cond> brad_idx;  
    array[J_brad_cond] int<lower=1, upper=J_size_cond> brad_size_mapper; // Maps bradford conditions to size conditions
    array[J_brad_cond] int<lower=1, upper=J_growth_cond> brad_growth_mapper;
    array[J_brad_cond] int<lower=1, upper=J_flow_cond> brad_flow_mapper; // Maps bradford conditions to size conditions


}

transformed data {
    vector[N_size] aspect_ratio = length ./ width;
}

parameters {
    // Size parameters
    vector<lower=0>[J_size_cond] width_mu;
    vector<lower=0>[J_size_cond] width_sigma;
    vector<lower=0>[J_size_cond] length_mu;
    vector<lower=0>[J_size_cond] length_sigma;
    vector<lower=0>[J_size_cond] volume_mu;
    vector<lower=0>[J_size_cond] volume_sigma;
    vector<lower=0>[J_size_cond] peri_volume_mu;
    vector<lower=0>[J_size_cond] peri_volume_sigma;
    vector<lower=0>[J_size_cond] alpha_mu;
    vector<lower=0>[J_size_cond] alpha_sigma;

    // Growth parameters
    vector[J_growth_cond] log_growth_mu;
    vector<lower=0>[J_growth_cond] growth_rates_sigma;

    // Bradford assay parameters
    real<lower=0> cal_slope;
    real<lower=0> cal_intercept;
    real<lower=0> cal_sigma;
    vector[J_brad_cond] log_prot_per_biomass_mu;
    vector<lower=0>[J_brad_cond] od595_per_biomass_sigma;

    // Biomass and total protein parameters
    real<lower=0> total_prot_0;
    real total_prot_slope;
    real<lower=0> total_prot_sigma;

    // Cell count numbers
    vector[J_flow_cond] log_flow_mu;
    real<lower=0> flow_sigma;
}

transformed parameters {
    vector<lower=0>[J_brad_cond] prot_per_biomass_mu = exp(log_prot_per_biomass_mu); 
    vector<lower=0>[J_growth_cond] growth_mu = exp(log_growth_mu);
    vector<lower=0>[J_flow_cond] flow_mu = exp(log_flow_mu);

}

model {

    // Bradford assay - Calibration
    cal_slope ~ normal(0, 0.1);
    cal_intercept ~ normal(0, 0.1);
    cal_sigma ~ std_normal();
    cal_od ~ normal(cal_intercept + cal_slope .* concentration, cal_sigma);

    // Bradford assay - Measurements
    log_prot_per_biomass_mu ~ std_normal();
    od595_per_biomass_sigma ~ normal(0, 0.1);

    // Likelihood
    log((brad_od595 - cal_intercept)./brad_od600)  ~ normal(log(cal_slope .* 
            prot_per_biomass_mu[brad_idx] ./ conv_factor), 
            od595_per_biomass_sigma[brad_idx]);

    // Growth rate measurements
    log_growth_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    log(growth_rates) ~ normal(log_growth_mu[growth_idx], growth_rates_sigma[growth_idx]);

    // Size measurement
    width_mu ~ std_normal(); 
    width_sigma ~ std_normal();
    length_mu ~ normal(0, 2);
    length_sigma ~ std_normal();
    volume_mu ~ normal(0, 2);
    volume_sigma ~ std_normal();
    peri_volume_mu ~ std_normal();
    peri_volume_sigma ~ std_normal();
    alpha_mu ~ normal(4, 1);
    alpha_sigma ~ normal(0, 0.1);
    width ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    length ~ normal(length_mu[size_idx], length_sigma[size_idx]);
    volume ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);
    peri_volume ~ normal(peri_volume_mu[size_idx], peri_volume_sigma[size_idx]);
    aspect_ratio ~ normal(alpha_mu[size_idx], alpha_sigma[size_idx]);

    // Cell count measurements
    flow_sigma ~ std_normal();
    log(flow_events) ~ normal(log_flow_mu[flow_idx], flow_sigma);

    //Biomass emeasurements
    total_prot_0 ~ normal(500, 100);
    total_prot_slope ~ normal(0, 100);
    total_protein ~ normal(total_prot_0 + total_prot_slope * total_protein_lam, total_prot_sigma);
}

generated quantities {
    //PPCs for bradford calibration
    vector[N_cal] od595_cal_rep;
    for (i in 1:N_cal) {
        od595_cal_rep[i] = normal_rng(cal_intercept + cal_slope * concentration[i], cal_sigma);
    }
    vector[N_brad] od595_meas_rep;
    for (i in 1:N_brad) {
        od595_meas_rep[i] = cal_intercept + brad_od600[i] * exp(normal_rng(log(cal_slope .* 
            prot_per_biomass_mu[brad_idx[i]] ./ conv_factor[i]), 
            od595_per_biomass_sigma[brad_idx[i]]));
    }
    vector[N_flow] flow_rep;
    for (i in 1:N_flow) {
        flow_rep[i] = exp(normal_rng(log_flow_mu[flow_idx[i]], flow_sigma));
    }
    //PPCs for total protein
    vector[N_prot] total_protein_rep;
    for (i in 1:N_prot) {
        total_protein_rep[i] = normal_rng(total_prot_0 + total_prot_slope * total_protein_lam[i], total_prot_sigma); 
    }

    //PPcs for growth rate
    vector[N_growth] growth_rate_rep;
    for (i in 1:N_growth) {
        growth_rate_rep[i] = exp(normal_rng(log_growth_mu[growth_idx[i]], growth_rates_sigma[growth_idx[i]]));
    }

    vector[J_brad_cond] phi_peri = prot_per_biomass_mu ./ (total_prot_0 + total_prot_slope * growth_mu[brad_growth_mapper]); 
    vector[J_brad_cond] m_peri = 1E9 * prot_per_biomass_mu ./ flow_mu[brad_flow_mapper];
    vector[J_brad_cond] rho_peri = m_peri ./ peri_volume_mu[brad_size_mapper];   
    vector[J_size_cond] width_rep;
    vector[J_size_cond] length_rep;
    vector[J_size_cond] volume_rep;
    vector[J_size_cond] peri_volume_rep;
    vector[J_size_cond] alpha_rep;


    for (i in 1:J_size_cond) {
        width_rep[i] = normal_rng(width_mu[i], width_sigma[i]);
        length_rep[i] = normal_rng(length_mu[i], length_sigma[i]);
        volume_rep[i] = normal_rng(volume_mu[i], volume_sigma[i]);
        peri_volume_rep[i] = normal_rng(peri_volume_mu[i], peri_volume_sigma[i]);
        alpha_rep[i] = normal_rng(alpha_mu[i], alpha_sigma[i]);
    }
}