data {
    // Periplasmic protein measurements
    int<lower=1>N_peri;
    int<lower=1>J_peri;
    array[N_peri] int<lower=1, upper=J_peri> peri_idx;
    vector<lower=0>[N_peri] peri_per_biomass;

    // Growth rate measurement
    int<lower=1> N_growth;
    int<lower=1> J_growth;
    array[N_growth] int<lower=1, upper=J_growth> growth_idx;
    vector<lower=0>[N_growth] growth_rates;

    // Size measurements
    int<lower=1> N_size;
    int<lower=1> J_size;
    array[N_size] int<lower=1, upper=J_size> size_idx;
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] aspect_ratios;
    vector<lower=0>[N_size] surface_areas;
    vector<lower=0>[N_size] surface_to_volumes;
    vector<lower=0>[N_size] volumes;

    // Flow cytometry measurements
    int<lower=1> N_flow;
    int<lower=1> J_flow;
    array[N_flow] int<lower=1, upper=J_flow> flow_idx;
    vector<lower=0>[N_flow] cells_per_biomass;
}

transformed data {
    vector[N_size] aspect_ratio_zeroed = aspect_ratios - 1;
}

parameters {
    // Periplasmic protein parameters
    vector<lower=0>[J_peri] peri_per_biomass_mu;
    vector<lower=0>[J_peri] peri_per_biomass_sigma;

    // Growth parameters
    vector<lower=0>[J_growth] growth_rate_mu;
    vector<lower=0>[J_growth] growth_rates_sigma;

    // Flow cytometry parameters
    vector<lower=0>[J_flow] cells_per_biomass_mu;
    vector<lower=0>[J_flow] cells_per_biomass_sigma;

    // Size parameters
    vector<lower=0>[J_size] width_mu;
    vector<lower=0>[J_size] width_sigma;
    vector<lower=0>[J_size] length_mu;
    vector<lower=0>[J_size] length_sigma;
    vector<lower=0>[J_size] surface_area_mu;
    vector<lower=0>[J_size] surface_area_sigma;
    vector<lower=0>[J_size] surface_to_volume_mu;
    vector<lower=0>[J_size] surface_to_volume_sigma;
    vector<lower=0>[J_size] volume_mu;
    vector<lower=0>[J_size] volume_sigma;
    vector<lower=0>[J_size] aspect_ratio_zeroed_mu;
    vector<lower=0>[J_size] alpha_zeroed_sigma;

}

transformed parameters {
    vector<lower=0>[J_size] aspect_ratio_mu = aspect_ratio_zeroed_mu + 1;
    
}

model {
    // Periplasmic protein inference model
    peri_per_biomass_mu ~ normal(0, 100);
    peri_per_biomass_sigma ~ normal(0, 10);
    peri_per_biomass ~ normal(peri_per_biomass_mu[peri_idx], peri_per_biomass_sigma[peri_idx]);

    // Growth rate inference model
    growth_rate_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    growth_rates ~ normal(growth_rate_mu[growth_idx], growth_rates_sigma[growth_idx]);

    // FLow cytometry inference model
    cells_per_biomass_mu ~ std_normal();
    cells_per_biomass_sigma ~ normal(0, 0.1);
    cells_per_biomass / 1E9 ~ normal(cells_per_biomass_mu[flow_idx], cells_per_biomass_sigma[flow_idx]);

    // Size inference model
    width_mu ~ std_normal();
    width_sigma ~ normal(0, 0.5);
    widths ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    
    length_mu ~ normal(0, 3);
    length_sigma ~ std_normal();
    lengths ~ normal(length_mu[size_idx], length_sigma[size_idx]);

    volume_mu ~ normal(0, 3);
    volume_sigma ~ normal(0, 0.5);
    volumes ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);

    surface_area_mu ~ normal(0, 3);
    surface_area_sigma ~ std_normal();
    surface_areas ~ normal(surface_area_mu[size_idx], surface_area_sigma[size_idx]);

    surface_to_volume_mu ~ normal(0, 3);
    surface_to_volume_sigma ~ std_normal();
    surface_to_volumes ~ normal(surface_to_volume_mu[size_idx], surface_to_volume_sigma[size_idx]);
    
    aspect_ratio_zeroed_mu ~ std_normal();
    alpha_zeroed_sigma ~ normal(0, 0.1);
    aspect_ratio_zeroed ~ normal(aspect_ratio_zeroed_mu[size_idx], alpha_zeroed_sigma[size_idx]);
}

generated quantities {
    // Calculate relevant properties
    vector<lower=0>[J_flow] N_cells = 1E9 * cells_per_biomass_mu;
    vector<lower=0>[J_peri] m_peri_mu ; 
    vector<lower=0>[J_peri] phi_peri_mu;
    vector<lower=0>[J_peri] rho_peri_mu;

    //PPCs
    vector<lower=0>[J_peri] peri_per_biomass_ppc;
    vector<lower=0>[J_flow] cells_per_biomass_ppc;
    vector<lower=0>[J_growth] growth_rates_ppc;
    vector<lower=0>[J_size] width_ppc;
    vector<lower=0>[J_size] length_ppc;
    vector<lower=0>[J_size] surface_area_ppc;
    vector<lower=0>[J_size] surface_to_volume_ppc;
    vector<lower=0>[J_size] volume_ppc;
    vector<lower=0>[J_size] aspect_ratio_ppc; 


    for (i in 1:J_peri) {
        m_peri_mu[i] = 1E9 * peri_per_biomass_mu[i] / N_cells[i];
        rho_peri_mu[i] = 1E9 * peri_per_biomass_mu[i] / (N_cells[i] * surface_area_mu[i] * 0.0276);
        peri_per_biomass_ppc[i] = normal_rng(peri_per_biomass_mu[i], peri_per_biomass_sigma[i]);
    }

    for (i in 1:J_flow) {
        cells_per_biomass_ppc[i] = 1E9 * normal_rng(cells_per_biomass_mu[i], cells_per_biomass_sigma[i]);
    }

    for (i in 1:J_growth) {
        growth_rates_ppc[i] = normal_rng(growth_rate_mu[i], growth_rates_sigma[i]);
    }

    for (i in 1:J_size) {
        width_ppc[i] = normal_rng(width_mu[i], width_sigma[i]);
        length_ppc[i] = normal_rng(length_mu[i], length_sigma[i]);
        volume_ppc[i] = normal_rng(volume_mu[i], volume_sigma[i]);
        surface_area_ppc[i] = normal_rng(surface_area_mu[i], surface_area_sigma[i]);
        surface_to_volume_ppc[i] = normal_rng(surface_to_volume_mu[i], surface_to_volume_sigma[i]);
        aspect_ratio_ppc[i] = 1 + normal_rng(aspect_ratio_zeroed_mu[i], alpha_zeroed_sigma[i]);
    }

}