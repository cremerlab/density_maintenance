data {
    // Total protein measurements
    int<lower=1>N_prot;
    int<lower=1>J_prot;
    array[N_prot] int<lower=1, upper=J_prot> prot_idx;
    vector<lower=0>[N_prot] prot_per_biomass;

    // Periplasmic protein measurements
    int<lower=1>N_peri;
    int<lower=1>J_peri;
    array[N_peri] int<lower=1, upper=J_peri> peri_idx;
    vector<lower=0>[N_peri] peri_per_biomass;

    // Membrane protein measurements
    int<lower=1>N_mem;
    int<lower=1>J_mem;
    array[N_mem] int<lower=1, upper=J_mem> mem_idx;
    vector<lower=0>[N_mem] mem_per_biomass;

    // Total RNA measurements
    int<lower=1>N_rna;
    int<lower=1>J_rna;
    array[N_rna] int<lower=1, upper=J_rna> rna_idx;
    vector<lower=0>[N_rna] rna_per_biomass;

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
    vector<lower=0>[N_size] volumes;

    // Flow cytometry measurements
    int<lower=1> N_flow;
    int<lower=1> J_flow;
    array[N_flow] int<lower=1, upper=J_flow> flow_idx;
    vector<lower=0>[N_flow] cells_per_biomass;
}

transformed data {
    // Perform centering operatiors to simplify prior choice
    vector[N_mem] mem_per_biomass_centered = (mem_per_biomass - mean(mem_per_biomass)) ./ sd(mem_per_biomass);    
    vector[N_peri] peri_per_biomass_centered = (peri_per_biomass - mean(peri_per_biomass)) ./ sd(peri_per_biomass);    
    vector[N_prot] prot_per_biomass_centered = (prot_per_biomass - mean(prot_per_biomass)) ./ sd(prot_per_biomass);    
    vector[N_rna] rna_per_biomass_centered = (rna_per_biomass - mean(rna_per_biomass)) ./ sd(rna_per_biomass);        
    vector[N_size] aspect_ratio_zeroed = aspect_ratios - 1;
}

parameters {
    // Total protein parameters
    vector[J_prot] prot_per_biomass_centered_mu;
    vector<lower=0>[J_prot] prot_per_biomass_centered_sigma;

    // Membrane protein parameters
    vector[J_mem] mem_per_biomass_centered_mu;
    vector<lower=0>[J_mem] mem_per_biomass_centered_sigma;

    // Periplasmic protein parameters
    vector[J_peri] peri_per_biomass_centered_mu;
    vector<lower=0>[J_peri] peri_per_biomass_centered_sigma;

    // RNA parameters
    vector[J_rna] rna_per_biomass_centered_mu;
    vector<lower=0>[J_rna] rna_per_biomass_centered_sigma;

    // Growth parameters
    vector<lower=0>[J_growth] growth_rates_mu;
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
    vector<lower=0>[J_size] volume_mu;
    vector<lower=0>[J_size] volume_sigma;
    vector<lower=0>[J_size] alpha_zeroed_mu;
    vector<lower=0>[J_size] alpha_zeroed_sigma;

}

transformed parameters {
    // Uncentering operations 
    vector[J_mem] mem_per_biomass_mu = mem_per_biomass_centered_mu .* sd(mem_per_biomass) + mean(mem_per_biomass);
    vector[J_peri] peri_per_biomass_mu = peri_per_biomass_centered_mu .* sd(peri_per_biomass) + mean(peri_per_biomass);
    vector[J_prot] prot_per_biomass_mu = prot_per_biomass_centered_mu .* sd(prot_per_biomass) + mean(prot_per_biomass);
    vector[J_rna] rna_per_biomass_mu = rna_per_biomass_centered_mu .* sd(rna_per_biomass) + mean(rna_per_biomass);
    vector[J_size] alpha_mu = alpha_zeroed_mu + 1;
}

model {
    // Total protein inference model
    prot_per_biomass_centered_mu ~ std_normal();
    prot_per_biomass_centered_sigma ~ std_normal();
    prot_per_biomass_centered ~ normal(prot_per_biomass_centered_mu[prot_idx], prot_per_biomass_centered_sigma[prot_idx]);

    // Membrane protein inference model
    mem_per_biomass_centered_mu ~ std_normal();
    mem_per_biomass_centered_sigma ~ std_normal();
    mem_per_biomass_centered ~ normal(mem_per_biomass_centered_mu[mem_idx], mem_per_biomass_centered_sigma[mem_idx]);

    // Periplasmic protein inference model
    peri_per_biomass_centered_mu ~ std_normal();
    peri_per_biomass_centered_sigma ~ std_normal();
    peri_per_biomass_centered ~ normal(peri_per_biomass_centered_mu[peri_idx], peri_per_biomass_centered_sigma[peri_idx]);

    // RNA inference model
    rna_per_biomass_centered_mu ~ std_normal();
    rna_per_biomass_centered_sigma ~ std_normal();
    rna_per_biomass_centered ~ normal(rna_per_biomass_centered_mu[rna_idx], rna_per_biomass_centered_sigma[rna_idx]);

    // Growth rate inference model
    growth_rates_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    growth_rates ~ normal(growth_rates_mu[growth_idx], growth_rates_sigma[growth_idx]);

    // FLow cytometry inference model
    cells_per_biomass_mu ~ std_normal();
    cells_per_biomass_sigma ~ std_normal();
    cells_per_biomass / 1E9 ~ normal(cells_per_biomass_mu[flow_idx], cells_per_biomass_sigma[flow_idx]);

    // Size inference model
    width_mu ~ std_normal();
    width_sigma ~ std_normal();
    widths ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    
    length_mu ~ normal(0, 3);
    length_sigma ~ std_normal();
    lengths ~ normal(length_mu[size_idx], length_sigma[size_idx]);

    volume_mu ~ normal(0, 3);
    volume_sigma ~ std_normal();
    volumes ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);

    surface_area_mu ~ normal(0, 3);
    surface_area_sigma ~ std_normal();
    surface_areas ~ normal(surface_area_mu[size_idx], surface_area_sigma[size_idx]);

    alpha_zeroed_mu ~ std_normal();
    alpha_zeroed_sigma ~ std_normal();
    aspect_ratios ~ normal(alpha_zeroed_mu[size_idx], alpha_zeroed_sigma[size_idx]);
}

generated quantities {
    // Calculate relevant properties
    vector<lower=0>[J_flow] N_cells = 1E9 * cells_per_biomass_mu;
    vector<lower=0>[J_mem] phi_mem = mem_per_biomass_mu ./ (N_cells .* prot_per_biomass_mu);
    vector<lower=0>[J_peri] phi_peri;
    vector<lower=0>[J_prot] phi_Rb = 0.4558 .* (rna_per_biomass_mu ./ prot_per_biomass_mu);
    vector<lower=0>[J_mem] rho_mem = mem_per_biomass_mu ./ (N_cells .* 2 .* surface_area_mu);
    vector<lower=0>[J_peri] rho_peri;

    //PPCs
    vector<lower=0>[J_prot] prot_per_biomass_ppc;
    vector<lower=0>[J_mem] mem_per_biomass_ppc;
    vector<lower=0>[J_peri] peri_per_biomass_ppc;
    vector<lower=0>[J_rna] rna_per_biomass_ppc;
    vector<lower=0>[J_flow] cells_per_biomass_ppc;
    vector<lower=0>[J_growth] growth_rates_ppc;
    vector<lower=0>[J_size] width_ppc;
    vector<lower=0>[J_size] length_ppc;
    vector<lower=0>[J_size] surface_area_ppc;
    vector<lower=0>[J_size] volume_ppc;
    vector<lower=0>[J_size] aspect_ratio_ppc; 

    for (i in 1:J_peri) {
        phi_peri[i] = peri_per_biomass_mu[i] / (N_cells[i] * prot_per_biomass_mu[i]);
        rho_peri[i] = peri_per_biomass_mu[i] / (N_cells[i] * 2 * surface_area_mu[i]);
        peri_per_biomass_ppc[i] = mean(peri_per_biomass) + sd(peri_per_biomass) * normal_rng(peri_per_biomass_centered_mu[peri_idx[i]], peri_per_biomass_centered_sigma[peri_idx[i]]);
    }
    for (i in 1:J_prot) {
        prot_per_biomass_ppc[i] = mean(prot_per_biomass) + sd(prot_per_biomass) * normal_rng(prot_per_biomass_centered_mu[prot_idx[i]], prot_per_biomass_centered_sigma[prot_idx[i]]);
    }

    for (i in 1:J_mem) {
        mem_per_biomass_ppc[i] = mean(mem_per_biomass) + sd(mem_per_biomass) * normal_rng(mem_per_biomass_centered_mu[mem_idx[i]], mem_per_biomass_centered_sigma[mem_idx[i]]);
    }

    for (i in 1:J_rna) {
        rna_per_biomass_ppc[i] = mean(rna_per_biomass) + sd(rna_per_biomass) * normal_rng(rna_per_biomass_centered_mu[rna_idx[i]], rna_per_biomass_centered_sigma[rna_idx[i]]);
    }

    for (i in 1:J_flow) {
        cells_per_biomass_ppc[i] = 1E9 * normal_rng(cells_per_biomass_mu[flow_idx[i]], cells_per_biomass_sigma[flow_idx[i]]);
    }

    for (i in 1:J_growth) {
        growth_rates_ppc[i] = normal_rng(growth_rates_mu[growth_idx[i]], growth_rates_sigma[growth_idx[i]]);
    }

    for (i in 1:J_size) {
        width_ppc[i] = normal_rng(width_mu[size_idx[i]], width_sigma[size_idx[i]]);
        length_ppc[i] = normal_rng(length_mu[size_idx[i]], length_sigma[size_idx[i]]);
        volume_ppc[i] = normal_rng(volume_mu[size_idx[i]], volume_sigma[size_idx[i]]);
        surface_area_ppc[i] = normal_rng(surface_area_mu[size_idx[i]], surface_area_sigma[size_idx[i]]);
        aspect_ratio_ppc[i] = 1 + normal_rng(alpha_zeroed_mu[size_idx[i]], alpha_zeroed_sigma[size_idx[i]]);
    }

}