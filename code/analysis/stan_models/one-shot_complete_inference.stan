data {
    //--------------------------------------------------------------------------
    //  Growth Measurements
    //--------------------------------------------------------------------------
    int<lower=1> N_growth; // Number of growth measurements
    int<lower=1> J_growth_cond; // Number of unique conditions
    int<lower=1> J_growth_curves; // Number of biological replicates
    array[J_growth_curves] int<lower=1, upper=J_growth_cond> growth_cond_idx;
    array[N_growth] int<lower=1, upper=J_growth_curves> growth_curve_idx;
    vector<lower=0>[N_growth] growth_time;
    vector<lower=0>[N_growth] growth_od;

    //--------------------------------------------------------------------------
    //  Bradford Assay Calibration Curve
    //--------------------------------------------------------------------------
    int<lower=1> N_cal; // Number of Bradford calibration curve measurements
    vector<lower=0>[N_cal] concentration; // Protein concentrations for cal curve
    vector<lower=0>[N_cal] cal_od; // OD595nm measurements for calibration curve

    //--------------------------------------------------------------------------
    //  Size Measurements
    //--------------------------------------------------------------------------
    int<lower=1> N_size; // Total number of size measurements
    int<lower=1> J_size_cond; // Total number of conditions
    array[N_size] int<lower=1, upper=J_size_cond> size_cond_idx;
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] length;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] peri_volume;
    vector<lower=0>[N_size] surface_area;
    vector<lower=0>[N_size] surface_area_volume;

    //--------------------------------------------------------------------------
    //  Bradford Assay Protein Measurements
    //--------------------------------------------------------------------------
    int<lower=1> N_brad;  // Total number of bradford measurements
    int<lower=1> J_brad_cond; // Total number of unique conditions
    array[N_brad] int<lower=1, upper=J_brad_cond> brad_cond_idx; // ID vector
    array[J_brad_cond] int<lower=1, upper=J_size_cond> brad_cond_mapper; // Maps bradford conditions to size conditions
    vector<lower=0>[N_brad] brad_od595; // OD595 measurements per sample
    vector<lower=0>[N_brad] brad_od600; // OD595 measurements per sample
    vector<lower=0>[N_brad] conv_factor; // Conversion factor from protein per biomass to OD

    //--------------------------------------------------------------------------
    //  Flow measurements
    //--------------------------------------------------------------------------
    int<lower=1> N_flow;
    array[N_flow] int<lower=1, upper=J_size_cond> flow_mapper;
    vector<lower=0>[N_flow] cells_per_biomass;

    //--------------------------------------------------------------------------
    //  Literature Biomass Measurements
    //--------------------------------------------------------------------------
    int<lower=1> N_biomass; // Total number of biomass measurements
    vector<lower=0>[N_biomass] biomass; // Biomass in ug / OD ml
}

transformed data {
    // -------------------------------------------------------------------------
    // Growth measurements
    // -------------------------------------------------------------------------
    vector[N_growth] log_growth_od = log(growth_od);

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    vector[N_biomass] biomass_centered = (biomass - mean(biomass)) ./ sd(biomass);

}

parameters { 
    // -------------------------------------------------------------------------
    // Bradford Assay Calibration Curve
    // -------------------------------------------------------------------------
    // Hyperparameters 
    vector<lower=0>[J_growth_cond] growth_mu;
    real growth_tau;

    // Level 1 hyper parameters
    vector<lower=0>[J_growth_curves] growth_mu_1_tilde;

    // Singular
    real<lower=0> growth_sigma;
    vector<lower=0>[J_growth_curves] growth_od_init;

    // -------------------------------------------------------------------------
    // Bradford Assay Calibration Curve
    // -------------------------------------------------------------------------
    // Hyperparameters for calibration
    real<lower=0> cal_slope;
    real<lower=0> cal_intercept;
    real<lower=0> cal_sigma;

    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------
    vector[J_brad_cond] log_prot_per_biomass_mu;
    vector<lower=0>[J_brad_cond] od595_per_biomass_sigma;

    // -------------------------------------------------------------------------
    // Size Measurements
    // -------------------------------------------------------------------------
    vector<lower=0>[J_size_cond] width_mu;
    vector<lower=0>[J_size_cond] width_sigma;
    vector<lower=0>[J_size_cond] length_mu;
    vector<lower=0>[J_size_cond] length_sigma;
    vector<lower=0>[J_size_cond] volume_mu;
    vector<lower=0>[J_size_cond] volume_sigma;
    vector<lower=0>[J_size_cond] peri_volume_mu;
    vector<lower=0>[J_size_cond] peri_volume_sigma;
    vector<lower=0>[J_size_cond] surface_area_mu;
    vector<lower=0>[J_size_cond] surface_area_sigma;
    vector<lower=0>[J_size_cond] surface_area_vol_mu;
    vector<lower=0>[J_size_cond] surface_area_vol_sigma;

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    real biomass_centered_mu; // Centered mean biomass
    real<lower=0> biomass_sigma; // scale parameter
        
    // -------------------------------------------------------------------------
    // Flow measurements
    // -------------------------------------------------------------------------
    real flow_prefactor;
    real<lower=0> flow_sigma;

} 

transformed parameters {
    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------
    vector[J_growth_curves] growth_mu_1 = growth_mu[growth_cond_idx] + growth_tau * growth_mu_1_tilde;
    vector[J_growth_curves] log_growth_od_init = log(growth_od_init);

    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------
    vector[J_brad_cond] prot_per_biomass_mu = exp(log_prot_per_biomass_mu);

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);

    // -------------------------------------------------------------------------
    // Flow measurements
    // ------------------------------------------------------------------------- 
    real<lower=0> flow_slope = flow_prefactor * biomass_mu;
    vector<lower=0>[N_flow] flow_mu = (flow_prefactor * biomass_mu * volume_mu[flow_mapper]);
}

model { 
    // -------------------------------------------------------------------------
    // Growth curves
    // -------------------------------------------------------------------------
    // Prior
    growth_mu ~ gamma(2.85, 3.30);
    growth_tau ~ std_normal();
    growth_mu_1_tilde ~ std_normal();
    growth_sigma ~ normal(0, 0.1);
    growth_od_init ~ normal(0, 0.01);

    // Likelihood
    log_growth_od ~ normal(log_growth_od_init[growth_curve_idx] + growth_mu_1[growth_curve_idx] .* growth_time, growth_sigma);

    // -------------------------------------------------------------------------
    // Bradford Assay Calibration Curve
    // -------------------------------------------------------------------------
    // Prior
    cal_slope ~ normal(0, 0.1);
    cal_intercept ~ normal(0, 0.1);
    cal_sigma ~ std_normal();

    // Likelihood
    cal_od ~ normal(cal_intercept + cal_slope .* concentration, cal_sigma);

    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------    
    // Prior
    log_prot_per_biomass_mu ~ std_normal();
    od595_per_biomass_sigma ~ std_normal();

    // Likelihood
    log((brad_od595 - cal_intercept)./brad_od600)  ~ normal(log(cal_slope .* prot_per_biomass_mu[brad_cond_idx] ./ conv_factor[brad_cond_idx]), od595_per_biomass_sigma[brad_cond_idx]);

    // -------------------------------------------------------------------------
    // Size Measurements
    // -------------------------------------------------------------------------    
    // Priors
    width_mu ~ std_normal();
    width_sigma ~ std_normal();
    length_mu ~ normal(0, 2);
    length_sigma ~ std_normal();
    volume_mu ~ normal(0, 2);
    volume_sigma ~ std_normal();
    peri_volume_mu ~ std_normal();
    peri_volume_sigma ~ std_normal();
    surface_area_mu ~ normal(0, 10);
    surface_area_sigma ~ std_normal();
    surface_area_vol_mu ~ normal(0, 10);
    surface_area_vol_sigma ~ std_normal();

    // Likelihood
    width ~ normal(width_mu[size_cond_idx], width_sigma[size_cond_idx]);
    length ~ normal(length_mu[size_cond_idx], length_sigma[size_cond_idx]);
    volume ~ normal(volume_mu[size_cond_idx], volume_sigma[size_cond_idx]);
    peri_volume ~ normal(peri_volume_mu[size_cond_idx], peri_volume_sigma[size_cond_idx]);
    surface_area  ~ normal(surface_area_mu[size_cond_idx], surface_area_sigma[size_cond_idx]);
    surface_area_volume  ~ normal(surface_area_vol_mu[size_cond_idx], surface_area_vol_sigma[size_cond_idx]);
 
    // -------------------------------------------------------------------------
    // Literature biomass Measurements
    // -------------------------------------------------------------------------    
    // Prior
    biomass_centered_mu ~ std_normal();
    biomass_sigma ~ std_normal();

    // Likelihood
    biomass_centered ~ normal(biomass_centered_mu, biomass_sigma);

    // -------------------------------------------------------------------------
    // Flow measurments
    // -------------------------------------------------------------------------    
    flow_sigma ~ normal(0, 0.1);
    flow_prefactor ~ std_normal();
    cells_per_biomass/1E9 ~ normal(1/flow_mu, flow_sigma);
 
}

generated quantities {
    // -------------------------------------------------------------------------
    // Bradford Assay Calibration Curve
    // -------------------------------------------------------------------------
    vector[N_growth] growth_od_rep;
    for (i in 1:N_growth) {
        growth_od_rep[i] = exp(normal_rng(log_growth_od_init[growth_curve_idx[i]] + growth_mu_1[growth_curve_idx[i]] * growth_time[i], growth_sigma));
    }
    
    // -------------------------------------------------------------------------
    // Bradford Assay Calibration Curve
    // -------------------------------------------------------------------------
    vector[N_cal] od595_calib_rep; 
    for (i in 1:N_cal) { 
        od595_calib_rep[i] = normal_rng(cal_intercept + cal_slope * concentration[i] , cal_sigma);
    } 

    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------
    vector[N_brad] od595_brad_rep; 
    for (i in 1:N_brad) { 
        od595_brad_rep[i] = cal_intercept + brad_od600[i] * exp(normal_rng(log(cal_slope * prot_per_biomass_mu[brad_cond_idx[i]] ./ conv_factor[i]), od595_per_biomass_sigma[brad_cond_idx[i]]));
    }   

    // -------------------------------------------------------------------------
    // Size Measurements
    // -------------------------------------------------------------------------
    vector[N_size] width_rep;
    vector[N_size] length_rep;
    vector[N_size] volume_rep;
    vector[N_size] peri_volume_rep;
    vector[N_size] surface_area_rep;
    vector[N_size] surface_area_vol_rep;

    for (i in 1:N_size) {
        width_rep[i] = normal_rng(width_mu[size_cond_idx[i]], width_sigma[size_cond_idx[i]]);
        length_rep[i] = normal_rng(length_mu[size_cond_idx[i]], length_sigma[size_cond_idx[i]]);
        volume_rep[i] = normal_rng(volume_mu[size_cond_idx[i]], volume_sigma[size_cond_idx[i]]);
        peri_volume_rep[i] = normal_rng(peri_volume_mu[size_cond_idx[i]], peri_volume_sigma[size_cond_idx[i]]);
        surface_area_rep[i] = normal_rng(surface_area_mu[size_cond_idx[i]], surface_area_sigma[size_cond_idx[i]]);
        surface_area_vol_rep[i] = normal_rng(surface_area_vol_mu[size_cond_idx[i]], surface_area_vol_sigma[size_cond_idx[i]]); 
    }

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    vector[N_flow] cells_per_biomass_rep;
    for (i in 1:N_flow) {
        cells_per_biomass_rep[i] = 1E9 * normal_rng(1/flow_mu[i], flow_sigma);
    }

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    vector[N_biomass] biomass_rep;
    for (i in 1:N_biomass) {
        biomass_rep[i] = normal_rng(biomass_centered_mu, biomass_sigma) * sd(biomass) + mean(biomass);
    }

    // -------------------------------------------------------------------------
    // Compute quantities
    // -------------------------------------------------------------------------
    vector[J_brad_cond] N_cells = 1 ./ (flow_slope .* volume_mu[brad_cond_mapper]);
    vector[J_brad_cond] periplasmic_density = prot_per_biomass_mu ./ (N_cells .* peri_volume_mu[brad_cond_mapper]);
    vector[J_brad_cond] cytoplasmic_density = biomass_mu ./ (N_cells .* volume_mu[brad_cond_mapper] .* peri_volume_mu[brad_cond_mapper]);
    vector[J_brad_cond] rho_ratio = periplasmic_density ./ cytoplasmic_density;
    vector[J_brad_cond] periplasmic_biomass_fraction = prot_per_biomass_mu / biomass_mu;
}