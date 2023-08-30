data {

    int<lower=1> J_pert;

    // Literature data for biomass density
    int<lower=0> N_drymass_density;
    vector<lower=0>[N_drymass_density] drymass_density;

    // Literature data for total biomass
    int<lower=0> N_drymass_total;
    vector<lower=0>[N_drymass_total] total_drymass;

    // Perturbation sizes
    int<lower=1> N_sizes;
    array[N_sizes] int<lower=1, upper=J_pert> size_idx;
    vector<lower=0>[N_sizes] widths;
    vector<lower=0>[N_sizes] lengths;
    vector<lower=0>[N_sizes] volumes;
    vector<lower=0>[N_sizes] surface_areas;
    vector<lower=0>[N_sizes] aspect_ratios;

    // Perturbation total protein quantification
    int<lower=1> N_prot;
    array[N_prot] int<lower=1, upper=J_pert> prot_idx;
    vector<lower=0>[N_prot] prot_per_biomass;

    // Perturbation membrane protein quantification
    int<lower=1> N_mem;
    array[N_mem] int<lower=1, upper=J_pert> mem_idx;
    vector<lower=0>[N_mem] mem_per_biomass;

    // Perturbation rna quantification
    int<lower=1> N_rna;
    array[N_rna] int<lower=1, upper=J_pert> rna_idx;
    vector<lower=0>[N_rna] rna_per_biomass;

    // Perturbation growth rate quantification
    int<lower=1> N_growth;
    array[N_growth] int<lower=1, upper=J_pert> growth_idx;
    vector<lower=0>[N_growth] growth_rates;

}

transformed data {
    vector[N_sizes] aspect_ratios_zeroed = aspect_ratios - 1;
}


parameters {

    // Size parameters
    vector<lower=0>[J_pert] width_mu;
    vector<lower=0>[J_pert] width_sigma;
    vector<lower=0>[J_pert] length_mu;
    vector<lower=0>[J_pert] length_sigma;
    vector<lower=0>[J_pert] surface_area_mu;
    vector<lower=0>[J_pert] surface_area_sigma;
    vector<lower=0>[J_pert] volume_mu;
    vector<lower=0>[J_pert] volume_sigma;
    vector<lower=0>[J_pert] aspect_ratio_zeroed_mu;
    vector<lower=0>[J_pert] aspect_ratio_zeroed_sigma;

    // Total protein parameters
    vector<lower=0>[J_pert] prot_per_biomass_mu;
    vector<lower=0>[J_pert] prot_per_biomass_sigma;

    // Membrane protein parameters
    vector<lower=0>[J_pert] mem_per_biomass_mu;
    vector<lower=0>[J_pert] mem_per_biomass_sigma;

    // RNA parameters
    vector<lower=0>[J_pert] rna_per_biomass_mu;
    vector<lower=0>[J_pert] rna_per_biomass_sigma;

    // Growth parameters
    vector<lower=0>[J_pert] growth_rate_mu;
    vector<lower=0>[J_pert] growth_rate_sigma;

    // Drymass density parameters
    real<lower=0> drymass_density_mu;
    real<lower=0> drymass_density_sigma;

    // Total biomass parameters
    real<lower=0> total_drymass_mu;
    real<lower=0> total_drymass_sigma;


}

transformed parameters {
    vector[J_pert] aspect_ratio_mu = 1 + aspect_ratio_zeroed_mu;
}

model {

    // Drymass density model
    drymass_density_mu ~ normal(300, 100);
    drymass_density_sigma ~ normal(0, 10);
    drymass_density ~ normal(drymass_density_mu, drymass_density_sigma);

    // Total drymass density
    total_drymass_mu ~ normal(500, 100);
    total_drymass_sigma ~ normal(0, 10);
    total_drymass ~ normal(total_drymass_mu, total_drymass_sigma);

    // Total protein inference model
    prot_per_biomass_mu ~ normal(0, 300);
    prot_per_biomass_sigma ~ std_normal();
    prot_per_biomass ~ normal(prot_per_biomass_mu[prot_idx], prot_per_biomass_sigma[prot_idx]);

    // Membrane protein inference model
    mem_per_biomass_mu ~ normal(0, 20);
    mem_per_biomass_sigma ~ std_normal();
    mem_per_biomass ~ normal(mem_per_biomass_mu[mem_idx], mem_per_biomass_sigma[mem_idx]);

    // RNA inference model
    rna_per_biomass_mu ~ normal(0, 200);
    rna_per_biomass_sigma ~ std_normal();
    rna_per_biomass ~ normal(rna_per_biomass_mu[rna_idx], rna_per_biomass_sigma[rna_idx]);

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

    aspect_ratio_zeroed_mu ~ std_normal();
    aspect_ratio_zeroed_sigma ~ normal(0, 0.1);
    aspect_ratios_zeroed ~ normal(aspect_ratio_zeroed_mu[size_idx], aspect_ratio_zeroed_sigma[size_idx]);

}

generated quantities {
    vector<lower=0>[J_pert] N_cells = total_drymass_mu ./ (volume_mu .* drymass_density_mu./1E9);
    vector<lower=0>[J_pert] phi_mem_mu = mem_per_biomass_mu ./ prot_per_biomass_mu; 
    vector<lower=0>[J_pert] rho_mem_mu = (1E9 .* mem_per_biomass_mu ./ (N_cells .* 2 .* surface_area_mu));
    vector<lower=0>[J_pert] phi_Rb_mu = 0.4558 .* rna_per_biomass_mu ./ prot_per_biomass_mu;
    vector<lower=0>[J_pert] surface_to_volume_mu = surface_area_mu ./ volume_mu;

    // PPCs
    vector<lower=0>[J_pert] width_ppc;
    vector<lower=0>[J_pert] length_ppc;
    vector<lower=0>[J_pert] volume_ppc;
    vector<lower=0>[J_pert] aspect_ratio_ppc;
    vector<lower=0>[J_pert] surface_area_ppc;
    vector<lower=0>[J_pert] prot_per_biomass_ppc;
    vector<lower=0>[J_pert] rna_per_biomass_ppc;
    vector<lower=0>[J_pert] mem_per_biomass_ppc;
    vector<lower=0>[J_pert] growth_rates_ppc;
    real<lower=0> total_drymass_ppc = normal_rng(total_drymass_mu, total_drymass_sigma);
    real<lower=0> drymass_density_ppc = normal_rng(drymass_density_mu, drymass_density_sigma); 
    
    for (i in 1:J_pert) {
        width_ppc[i] = normal_rng(width_mu[i], width_sigma[i]);
        length_ppc[i] = normal_rng(length_mu[i], length_sigma[i]);
        volume_ppc[i] = normal_rng(volume_mu[i], volume_sigma[i]);
        surface_area_ppc[i] = normal_rng(surface_area_mu[i], surface_area_sigma[i]);
        aspect_ratio_ppc[i] = 1 + normal_rng(aspect_ratio_zeroed_mu[i], aspect_ratio_zeroed_sigma[i]);
        growth_rates_ppc[i] = normal_rng(growth_rate_mu[i], growth_rate_sigma[i]);
        prot_per_biomass_ppc[i] = normal_rng(prot_per_biomass_mu[i], prot_per_biomass_sigma[i]);
        mem_per_biomass_ppc[i] = normal_rng(mem_per_biomass_mu[i], mem_per_biomass_sigma[i]);
        rna_per_biomass_ppc[i] = normal_rng(rna_per_biomass_mu[i], rna_per_biomass_sigma[i]);
    }
}