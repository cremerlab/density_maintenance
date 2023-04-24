data {
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
    vector<lower=1>[N_size] aspect_ratio;

    //--------------------------------------------------------------------------
    //  Aggregated growth rates to infer minimum width and minimum length.
    //--------------------------------------------------------------------------
    int<lower=1> J_growth;
    int<lower=1> N_growth;
    array[N_growth] int<lower=1, upper=J_growth> growth_cond_idx;
    vector<lower=0>[N_growth] growth_rates;
   array[N_growth] int<lower=1, upper=J_size_cond> growth_size_idx;

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
    int<lower=1> N_prot;
    vector<lower=0>[N_biomass] biomass; // Biomass in ug / OD ml
    vector<lower=0>[N_prot] total_protein;
    vector<lower=0>[N_prot] total_protein_growth_rates;
}

transformed data {
    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    vector[N_biomass] biomass_centered = (biomass - mean(biomass)) ./ sd(biomass);

}

parameters { 

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
    vector<lower=1>[J_size_cond] aspect_ratio_mu;
    vector<lower=0>[J_size_cond] aspect_ratio_sigma;

    // Size parameters for growth rate dependences
    real<lower=0> width_min;
    real<lower=0> width_slope;
    real<lower=0> width_lam_sigma; 
    real<lower=0> sav_min;
    real<lower=0> sav_slope;
    real<lower=0> sav_sigma;
    real<lower=0> alpha_min;
    real alpha_slope;
    real<lower=0> alpha_sigma;
    real<lower=0> length_min;
    real<lower=0> length_slope;
    real<lower=0> length_lam_sigma;

    // -------------------------------------------------------------------------
    // Growth rate measurements
    // -------------------------------------------------------------------------
    vector<lower=0>[J_growth] growth_mu;
    vector<lower=0>[J_growth] growth_rates_sigma;

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    real biomass_centered_mu; // Centered mean biomass
    real<lower=0> biomass_sigma; // scale parameter
    real<lower=0> total_protein_mu;
    real<lower=0> total_protein_sigma;
        
    // -------------------------------------------------------------------------
    // Flow measurements
    // -------------------------------------------------------------------------
    real<lower=0> flow_prefactor;
    real<lower=0> flow_sigma;
} 

transformed parameters {

    // -------------------------------------------------------------------------
    // Literature Biomass Measurements
    // -------------------------------------------------------------------------
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);

    // -------------------------------------------------------------------------
    // Bradford Assay Protein Measurements
    // -------------------------------------------------------------------------
    vector<lower=0>[J_brad_cond] prot_per_biomass_mu = exp(log_prot_per_biomass_mu);

    // -------------------------------------------------------------------------
    // Flow measurements
    // ------------------------------------------------------------------------- 
    real<lower=0> flow_slope = flow_prefactor * biomass_mu;
    vector<lower=0>[N_flow] flow_mu = (flow_prefactor * biomass_mu * volume_mu[flow_mapper]);

       }



model { 
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
    od595_per_biomass_sigma ~ normal(0, 0.1);

    // Likelihood
    log((brad_od595 - cal_intercept)./brad_od600)  ~ normal(log(cal_slope .* 
            prot_per_biomass_mu[brad_cond_idx] ./ conv_factor[brad_cond_idx]), 
            od595_per_biomass_sigma[brad_cond_idx]);

    // -------------------------------------------------------------------------
    // Growth rate measurements
    // -------------------------------------------------------------------------    
    growth_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    growth_rates ~ normal(growth_mu[growth_cond_idx], growth_rates_sigma[growth_cond_idx]);

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
    aspect_ratio_mu ~ normal(1, 4);
    aspect_ratio_sigma ~ std_normal();
    width_min ~ normal(0, 0.5);
    width_slope ~ normal(0, 0.5);
    width_lam_sigma ~ normal(0, 0.5);
    length_min ~ normal(0, 2);
    length_slope ~ normal(0, 4);
    length_lam_sigma ~ normal(0, 0.5);
    alpha_min ~ normal(1, 2);
    alpha_slope ~ std_normal();
    alpha_sigma ~ normal(0, 0.5);

    // Likelihood
    width ~ normal(width_mu[size_cond_idx], width_sigma[size_cond_idx]);
    length ~ normal(length_mu[size_cond_idx], length_sigma[size_cond_idx]);
    volume ~ normal(volume_mu[size_cond_idx], volume_sigma[size_cond_idx]);
    peri_volume ~ normal(peri_volume_mu[size_cond_idx], peri_volume_sigma[size_cond_idx]);
    surface_area  ~ normal(surface_area_mu[size_cond_idx], surface_area_sigma[size_cond_idx]);
    surface_area_volume  ~ normal(surface_area_vol_mu[size_cond_idx], surface_area_vol_sigma[size_cond_idx]);
    aspect_ratio  ~ normal(aspect_ratio_mu[size_cond_idx], aspect_ratio_sigma[size_cond_idx]);

    // -------------------------------------------------------------------------
    // Literature biomass Measurements
    // -------------------------------------------------------------------------    
    // Prior
    biomass_centered_mu ~ std_normal();
    biomass_sigma ~ std_normal();
    total_protein_mu ~ normal(300, 100);
    total_protein_sigma ~ normal(0, 10);

    // Likelihood
    biomass_centered ~ normal(biomass_centered_mu, biomass_sigma);
    total_protein ~ normal(total_protein_mu, total_protein_sigma);

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
    vector[N_size] aspect_ratio_rep;

    for (i in 1:N_size) {
        width_rep[i] = normal_rng(width_mu[size_cond_idx[i]], width_sigma[size_cond_idx[i]]);
        length_rep[i] = normal_rng(length_mu[size_cond_idx[i]], length_sigma[size_cond_idx[i]]);
        volume_rep[i] = normal_rng(volume_mu[size_cond_idx[i]], volume_sigma[size_cond_idx[i]]);
        peri_volume_rep[i] = normal_rng(peri_volume_mu[size_cond_idx[i]], peri_volume_sigma[size_cond_idx[i]]);
        surface_area_rep[i] = normal_rng(surface_area_mu[size_cond_idx[i]], surface_area_sigma[size_cond_idx[i]]);
        surface_area_vol_rep[i] = normal_rng(surface_area_vol_mu[size_cond_idx[i]], surface_area_vol_sigma[size_cond_idx[i]]); 
        aspect_ratio_rep[i] = normal_rng(aspect_ratio_mu[size_cond_idx[i]], aspect_ratio_sigma[size_cond_idx[i]]); 
    }

    // -------------------------------------------------------------------------
    // Flow cytometry cell measurements
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

    vector[N_prot] total_protein_rep;
    for (i in 1:N_prot) {
        total_protein_rep[i] = normal_rng(total_protein_mu, total_protein_sigma);
    }

    // -------------------------------------------------------------------------
    // Compute quantities
    // -------------------------------------------------------------------------
    vector<lower=0>[J_brad_cond] N_cells = 1E9 ./ (flow_slope .* volume_mu[brad_cond_mapper]);
    vector<lower=0>[J_brad_cond] phi_peri = prot_per_biomass_mu ./ total_protein_mu;
    vector<lower=0>[J_brad_cond] rho_peri = prot_per_biomass_mu ./ (N_cells .* peri_volume_mu[brad_cond_mapper]); 
    vector<lower=0>[J_brad_cond] peri_prot_per_cell = prot_per_biomass_mu ./ N_cells;
}