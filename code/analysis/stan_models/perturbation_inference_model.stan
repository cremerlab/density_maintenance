data { 
    // Measurement dimensions
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_biomass;
    int<lower=1> N_cal;
    int<lower=1> N_brad;
    int<lower=1> N_flow;
    int<lower=1> N_growth;

    // Condition dimensions
    int<lower=1> J_size_cond;
    int<lower=1> J_brad_cond; 
    int<lower=1> J_growth_cond;

    // ID vectors
    array[N_growth] int<lower=1, upper=J_growth_cond> growth_idx;  
    array[N_size] int<lower=1, upper=J_size_cond> size_idx;  
    array[N_brad] int<lower=1, upper=J_brad_cond> brad_idx;  
    array[J_brad_cond] int<lower=1, upper=J_size_cond> brad_mapper; // Maps bradford conditions to size conditions


    // Size measurements
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] length;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] peri_volume;

    // Growth measurements
    vector<lower=0>[N_growth] growth_rates;

    // Total protein measurements
    vector<lower=0>[N_prot] total_protein;
    // vector<lower=0>[N_prot] total_protein_lam;

    // Total biomass measurement
    vector<lower=0>[N_biomass] biomass;

    // Flow measurements
    vector<lower=0>[N_flow] flow_events;
    array[N_flow] int<lower=1, upper=J_size_cond> flow_mapper;

    // Bradford assay measurements
    vector<lower=0>[N_brad] brad_od595; 
    vector<lower=0>[N_brad] brad_od600; 
    vector<lower=0>[N_brad] conv_factor;
    vector<lower=0>[N_cal] concentration; 
    vector<lower=0>[N_cal] cal_od;


}

transformed data {
    vector[N_biomass] biomass_centered = (biomass - mean(biomass)) / sd(biomass);
    vector[N_prot] total_prot_centered = (total_protein - mean(total_protein)) / sd(total_protein);
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

    // Growth parameters
    vector<lower=0>[J_growth_cond] growth_mu;
    vector<lower=0>[J_growth_cond] growth_rates_sigma;

    // Bradford assay parameters
    real<lower=0> cal_slope;
    real<lower=0> cal_intercept;
    real<lower=0> cal_sigma;
    vector[J_brad_cond] log_prot_per_biomass_mu;
    vector<lower=0>[J_brad_cond] od595_per_biomass_sigma;
    // Biomass and total protein parameters
    real biomass_centered_mu;
    real<lower=0> biomass_sigma;
    real total_prot_centered_mu;
    real<lower=0> total_prot_sigma;

    // Cell count numbers
    real<lower=0> flow_prefactor;
    real<lower=0> flow_sigma;

}

transformed parameters {
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);
    real<lower=0> total_protein_mu = total_prot_centered_mu * sd(total_protein) + mean(total_protein);
    vector<lower=0>[J_brad_cond] prot_per_biomass_mu = exp(log_prot_per_biomass_mu);
    real<lower=0> flow_slope = flow_prefactor * biomass_mu;
    vector<lower=0>[N_flow] flow_mu = flow_slope * volume_mu[flow_mapper];
    // vector<lower=1>[J_size_cond] alpha = length_mu ./ width_mu;

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
            prot_per_biomass_mu[brad_idx] ./ conv_factor[brad_idx]), 
            od595_per_biomass_sigma[brad_idx]);

    // Growth rate measurements
    growth_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    growth_rates ~ normal(growth_mu[growth_idx], growth_rates_sigma[growth_idx]);

    // Size measurement
    width_mu ~ std_normal(); 
    width_sigma ~ std_normal();
    length_mu ~ normal(0, 2);
    length_sigma ~ std_normal();
    volume_mu ~ normal(0, 2);
    volume_sigma ~ std_normal();
    peri_volume_mu ~ std_normal();
    peri_volume_sigma ~ std_normal();
    width ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    length ~ normal(length_mu[size_idx], length_sigma[size_idx]);
    volume ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);
    peri_volume ~ normal(peri_volume_mu[size_idx], peri_volume_sigma[size_idx]);


    // Cell count measurements
    flow_sigma ~ normal(0, 0.1);
    flow_prefactor ~ std_normal();
    flow_events/1E9 ~ normal(1/flow_mu, flow_sigma);

    //Biomass emeasurements
    biomass_centered ~ normal(biomass_centered_mu, biomass_sigma);
    total_prot_centered ~ normal(total_prot_centered_mu, total_prot_sigma);
}

generated quantities {
    vector[J_size_cond] alpha = length_mu ./ width_mu;
    vector[J_brad_cond] phi_peri = prot_per_biomass_mu ./ total_protein_mu; 
    vector[J_brad_cond] N_cells = 1E9 / (flow_prefactor * biomass_mu * volume_mu[brad_mapper]);
    vector[J_brad_cond] m_peri = prot_per_biomass_mu ./ N_cells;
    vector[J_brad_cond] rho_peri = m_peri ./ peri_volume_mu[brad_mapper];
}