data {
    int<lower=1> J_cond;

    // Growth rate data
    int<lower=1> N_growth;
    array[N_growth] int<lower=1, upper=J_cond> growth_idx;
    vector<lower=0>[N_growth] growth_rates;

    // Total protein data
    int<lower=1> N_phi;
    array[N_phi] int<lower=1, upper=J_cond> phi_idx;
    vector<lower=0>[N_phi] prot_per_biomass;
    vector<lower=0>[N_phi] rna_per_biomass;
    vector<lower=0>[N_phi] phiRb;
    

    // Size data
    int<lower=1> N_size;
    array[N_size] int<lower=1, upper=J_cond> size_idx;
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] volumes;
    vector<lower=0>[N_size] surface_areas;
    vector<lower=0>[N_size] surface_to_volume;
    vector<lower=0>[N_size] aspect_ratios;
}

transformed data{
    vector<lower=0>[N_size] aspect_ratios_zeroed = aspect_ratios - 1;
}

parameters {
    // Growth rate params
    vector<lower=0>[J_cond] width_mu;
    vector<lower=0>[J_cond] width_sigma;
    vector<lower=0>[J_cond] length_mu;
    vector<lower=0>[J_cond] length_sigma;
    vector<lower=0>[J_cond] surface_area_mu;
    vector<lower=0>[J_cond] surface_area_sigma;
    vector<lower=0>[J_cond] surface_to_volume_mu;
    vector<lower=0>[J_cond] surface_to_volume_sigma; 
    vector<lower=0>[J_cond] volume_mu;
    vector<lower=0>[J_cond] volume_sigma;
    vector<lower=0>[J_cond] aspect_ratio_zeroed_mu;
    vector<lower=0>[J_cond] aspect_ratio_zeroed_sigma;

    // Total RNA and protein parameters
    vector<lower=0>[J_cond] prot_per_biomass_mu;
    vector<lower=0>[J_cond] prot_per_biomass_sigma;
    vector<lower=0>[J_cond] rna_per_biomass_mu;
    vector<lower=0>[J_cond] rna_per_biomass_sigma;
    vector<lower=0, upper=1>[J_cond] phiRb_mu;
    vector<lower=0>[J_cond] phiRb_sigma;

    // Growth parameters
    vector<lower=0>[J_cond] growth_rate_mu;
    vector<lower=0>[J_cond] growth_rate_sigma;
}

transformed parameters {
    vector<lower=1>[J_cond] aspect_ratio_mu = 1 + aspect_ratio_zeroed_mu;
}

model { 

    // Total protein inference model
    prot_per_biomass_mu ~ normal(0, 300);
    prot_per_biomass_sigma ~ std_normal();
    prot_per_biomass ~ normal(prot_per_biomass_mu[phi_idx], prot_per_biomass_sigma[phi_idx]);

    // RNA inference model
    rna_per_biomass_mu ~ normal(0, 200);
    rna_per_biomass_sigma ~ std_normal();
    rna_per_biomass ~ normal(rna_per_biomass_mu[phi_idx], rna_per_biomass_sigma[phi_idx]);

    // Ribosomal allocation
    phiRb_mu  ~ normal(0.15, 0.1);
    phiRb_sigma ~ normal(0, 0.01);
    phiRb ~ normal(phiRb_mu[phi_idx], phiRb_sigma[phi_idx]);
    
    // Growth rate inference model
    growth_rate_mu ~ std_normal();
    growth_rate_sigma ~ normal(0, 0.1);
    growth_rates ~ normal(growth_rate_mu[growth_idx], growth_rate_sigma[growth_idx]);

    // Size inference model
    width_mu ~ std_normal();
    width_sigma ~ normal(0, 0.1);
    widths ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    
    length_mu ~ normal(0, 3);
    length_sigma ~ std_normal();
    lengths ~ normal(length_mu[size_idx], length_sigma[size_idx]);

    volume_mu ~ normal(0, 3);
    volume_sigma ~ normal(0, 0.1);
    volumes ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);

    surface_area_mu ~ normal(0, 3);
    surface_area_sigma ~ normal(0, 0.1);
    surface_areas ~ normal(surface_area_mu[size_idx], surface_area_sigma[size_idx]);

    surface_to_volume_mu ~ normal(0, 10);
    surface_to_volume_sigma ~ normal(0, 0.1);
    surface_to_volume ~ normal(surface_to_volume_mu[size_idx], surface_to_volume_sigma[size_idx]);

    aspect_ratio_zeroed_mu ~ std_normal();
    aspect_ratio_zeroed_sigma ~ normal(0, 0.1);
    aspect_ratios_zeroed ~ normal(aspect_ratio_zeroed_mu[size_idx], aspect_ratio_zeroed_sigma[size_idx]);
}

generated quantities {
    vector<lower=0>[J_cond] width_ppc;
    vector<lower=0>[J_cond] length_ppc;
    vector<lower=0>[J_cond] volume_ppc;
    vector<lower=0>[J_cond] aspect_ratio_ppc;
    vector<lower=0>[J_cond] surface_area_ppc;
    vector<lower=0>[J_cond] surface_to_volume_ppc;
    vector<lower=0>[J_cond] prot_per_biomass_ppc;
    vector<lower=0>[J_cond] rna_per_biomass_ppc;
    vector<lower=0>[J_cond] phiRb_ppc;
    vector<lower=0>[J_cond] growth_rates_ppc;

    for (i in 1:J_cond) {
        width_ppc[i] = normal_rng(width_mu[i], width_sigma[i]);
        length_ppc[i] = normal_rng(length_mu[i], length_sigma[i]);
        volume_ppc[i] = normal_rng(volume_mu[i], volume_sigma[i]);
        surface_area_ppc[i] = normal_rng(surface_area_mu[i], surface_area_sigma[i]);
        surface_to_volume_ppc[i] = normal_rng(surface_to_volume_mu[i], surface_to_volume_sigma[i]);
        aspect_ratio_ppc[i] = 1 + normal_rng(aspect_ratio_zeroed_mu[i], aspect_ratio_zeroed_sigma[i]);
        prot_per_biomass_ppc[i] = normal_rng(prot_per_biomass_mu[i], prot_per_biomass_sigma[i]);
        rna_per_biomass_ppc[i] = normal_rng(rna_per_biomass_mu[i], rna_per_biomass_sigma[i]);
        phiRb_ppc[i] = normal_rng(phiRb_mu[i], phiRb_sigma[i]);
        growth_rates_ppc[i] = normal_rng(growth_rate_mu[i], growth_rate_sigma[i]);
    }

}