data {

    // Protein calibration curve data 
    int<lower=1> N_bca_cal;
    int<lower=1> N_brad_cal;
    int<lower=1> N_biuret_cal;
    vector<lower=0>[N_bca_cal] bca_std_conc;
    vector<lower=0>[N_bca_cal] bca_cal;
    vector<lower=0>[N_brad_cal] brad_std_conc;
    vector<lower=0>[N_brad_cal] brad_cal;
    vector<lower=0>[N_biuret_cal] biuret_std_conc;
    vector<lower=0>[N_biuret_cal] biuret_cal;

    // Membrane protein measurement data
    int<lower=1> N_mem;
    int<lower=1> J_mem;
    array[N_mem] int<lower=1, upper=J_mem> mem_idx;
    vector<lower=0>[N_mem] mem_conv_factor;
    vector<lower=0>[N_mem] mem_od562nm;
    vector<lower=0>[N_mem] mem_od600nm;

    // Periplasmic protein measurement data
    int<lower=1> N_peri;
    int<lower=1> J_peri;
    array[N_peri] int<lower=1, upper=J_peri> peri_idx;
    vector<lower=0>[N_peri] peri_conv_factor;
    vector[N_peri] peri_od595nm;
    vector[N_peri] peri_od600nm;

    // Total protein measurement data
    int<lower=1> N_prot;
    int<lower=1> J_prot;
    array[N_prot] int<lower=1, upper=J_prot> prot_idx;
    vector<lower=0>[N_prot] prot_conv_factor;
    vector<lower=0>[N_prot] prot_od555nm;
    vector<lower=0>[N_prot] prot_od600nm;

    // Flow data
    int<lower=1> N_flow;
    int<lower=1> J_flow;
    array[N_flow] int<lower=1, upper=J_flow> flow_idx;
    vector<lower=0>[N_flow] cell_count;

    // Growth data
    int<lower=1> N_growth;
    int<lower=1> J_growth;
    array[N_growth] int<lower=1, upper=J_growth> growth_idx;
    vector<lower=0>[N_growth] growth_rates;

    // Size data
    int<lower=1> N_size;
    int<lower=1> J_size;
    array[N_size] int<lower=1, upper=J_size> size_idx;
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] surface_area;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] aspect_ratio;

    // // Protein per cell data;
    // int<lower=1> N_prot;
    // vector<lower=0>[N_prot] prot_growth_rate;
    // vector<lower=0>[N_prot] prot_per_cell;
}

parameters { 
    // Calibration curve parameters
    real bca_intercept;
    real bca_slope;
    real<lower=0> bca_sigma;
    real brad_intercept;
    real brad_slope;
    real<lower=0> brad_sigma;
    real biuret_intercept;
    real biuret_slope;
    real<lower=0> biuret_sigma;

    // Membrane protein measurement parameters
    vector[J_mem] log_od562nm_mem_mu;
    vector<lower=0>[J_mem] log_od562nm_mem_sigma;

    // Periplasmic protein measurement parameters
    vector<lower=0>[J_peri] peri_per_biomass_mu;
    vector<lower=0>[J_peri] peri_per_biomass_sigma; 

    // Total protein measurement parameters;
    vector[J_prot]log_od555nm_prot_mu;
    vector<lower=0>[J_prot] log_od555nm_prot_sigma;

    // Flow parameters
    vector<lower=0>[J_flow] cells_per_biomass_mu; 
    vector<lower=0>[J_flow] cells_per_biomass_sigma;

    // Growth parameters
    vector<lower=0>[J_growth] lam_mu;
    vector<lower=0>[J_growth] lam_sigma;

    // Size parameters
    vector<lower=0>[J_size] width_mu;
    vector<lower=0>[J_size] width_sigma;
    vector<lower=0>[J_size] length_mu;
    vector<lower=0>[J_size] length_sigma;
    vector<lower=0>[J_size] surface_area_mu;
    vector<lower=0>[J_size] surface_area_sigma;
    vector<lower=0>[J_size] volume_mu;
    vector<lower=0>[J_size] volume_sigma;
    vector<lower=1>[J_size] alpha_mu;
    vector<lower=0>[J_size] alpha_sigma;


    // // Protein per cell parameters
    // real log_protein_intercept;
    // real log_protein_slope;
    // real log_protein_sigma;

}

transformed parameters {

    // Various protein per biomass conversions
    vector<lower=0>[J_mem] mem_per_biomass_mu = exp(log_od562nm_mem_mu);
    // vector<lower=0>[J_peri] peri_per_biomass_mu = exp(log_od595nm_peri_mu);
    vector<lower=0>[J_prot] prot_per_biomass_mu = exp(log_od555nm_prot_mu);

    // Size parameters
    vector[J_size] alpha = length_mu ./ width_mu;

    // // Membrane protein parameters
    // vector[J_mem] mem_conc = (exp(log_od562nm_mem_mu) - bca_intercept) / bca_slope;
    // vector[J_mem] phi_mem = (mem_conc  * 1E9) ./ (exp(log_protein_intercept + log_protein_slope .* lam_mu) .* cells_per_biomass_mu .* 1E9);
    // vector[J_mem] rho_mem =  (mem_conc  * 1E9) ./ (cells_per_biomass_mu .* 1E9 .* 2 .* surface_area_mu);

    // // Periplasmic protein parameters
    // vector[J_peri] peri_conc = (exp(log_od595nm_peri_mu) - brad_intercept) / brad_slope;
    // vector[J_peri] phi_peri;
    // vector[J_peri] rho_peri;
    // vector[J_peri] m_peri;
    // for (i in 1:J_peri) { 
    //     phi_peri[i] = peri_conc[i] * 1E9 / (exp(log_protein_intercept + log_protein_slope * lam_mu[i]) * cells_per_biomass_mu[i] * 1E9);
    //     rho_peri[i] = peri_conc[i] * 1E9 / (cells_per_biomass_mu[i] * 1E9 * surface_area_mu[i] * 0.0246); // 0.0246 is periplasmic widt spacing, around 25 nm.
    //     m_peri[i] = peri_conc[i] * 1E9 / (cells_per_biomass_mu[i] * 1E9);
    // }
}

model { 

    // BCA calibration curve
    bca_intercept ~ std_normal();
    bca_slope ~ std_normal();
    bca_sigma ~ std_normal();
    bca_cal ~ normal(bca_intercept + bca_slope .* bca_std_conc, bca_sigma);

    // Bradford calibration curve
    brad_intercept ~ std_normal();
    brad_slope ~ std_normal();
    brad_sigma ~ std_normal();
    brad_cal ~ normal(brad_intercept + brad_slope .* brad_std_conc, brad_sigma);

    // Biuret calibration curve
    biuret_intercept ~ std_normal();
    biuret_slope ~ std_normal();
    biuret_sigma ~ std_normal();
    biuret_cal ~ normal(biuret_intercept + biuret_slope .* biuret_std_conc, biuret_sigma);

    // Membrane measurement model
    log_od562nm_mem_mu ~ std_normal();
    log_od562nm_mem_sigma ~ normal(0, 0.1);
    // log((mem_od562nm - bca_intercept) ./ mem_od600nm) ~ normal(log(bca_slope .* mem_per_biomass_mu[mem_idx] ./ mem_conv_factor[mem_idx]), log_od562nm_mem_sigma[mem_idx]);

    // Periplasmic measurement modelo
    // log_od595nm_peri_mu ~ std_normal();
    // log_od595nm_peri_sigma ~ normal(0, 0.1);
    peri_per_biomass_mu ~ std_normal();
    peri_per_biomass_sigma ~ std_normal();
    peri_conv_factor .* (peri_od595nm - brad_intercept) ./ (peri_od600nm .* brad_slope) ~ normal(peri_per_biomass_mu[peri_idx], peri_per_biomass_sigma[peri_idx]);
    // log((prot_od555nm - biuret_intercept) ./ prot_od600nm) ~ normal(log(biuret_slope * prot_per_biomass_mu[prot_idx] ./ prot_conv_factor[prot_idx]), log_od555nm_prot_sigma[prot_idx]);


    // Total protein measurement model
    log_od555nm_prot_mu ~ std_normal();
    log_od555nm_prot_sigma ~ std_normal();
    // peri_conv_factor .* (peri_od595nm - brad_intercept) ./ (peri_od600nm * brad_slope) ~ normal(peri_)
    // log((prot_od555nm - biuret_intercept) ./ prot_od600nm) ~ normal(log(biuret_slope * prot_per_biomass_mu[prot_idx] ./ prot_conv_factor[prot_idx]), log_od555nm_prot_sigma[prot_idx]);

    // Flow model
    cells_per_biomass_mu ~ std_normal();
    cells_per_biomass_sigma ~ normal(0, 0.1);
    cell_count / 1E9 ~ normal(cells_per_biomass_mu[flow_idx],  cells_per_biomass_sigma[flow_idx]);

    // Growth model
    lam_mu ~ std_normal();
    lam_sigma ~ normal(0, 0.1);
    growth_rates ~ normal(lam_mu[growth_idx], lam_sigma[growth_idx]);

    // Size model
    width_mu ~ std_normal();
    length_mu ~ normal(0, 2);
    surface_area_mu ~ normal(0, 3);
    volume_mu ~ normal(0, 3);
    alpha_mu ~ normal(1, 2);
    width_sigma ~ std_normal();
    length_sigma ~ std_normal();
    volume_sigma ~ std_normal();
    surface_area_sigma ~ std_normal();
    alpha_sigma ~ std_normal();
    widths ~ normal(width_mu[size_idx], width_sigma[size_idx]);
    lengths ~ normal(length_mu[size_idx], length_sigma[size_idx]);
    volume ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);
    aspect_ratio ~ normal(alpha_mu[size_idx], alpha_sigma[size_idx]);
    surface_area ~ normal(surface_area_mu[size_idx], surface_area_sigma[size_idx]);

    // // Protein per cell model
    // log_protein_intercept ~ std_normal();
    // log_protein_slope ~ std_normal();
    // log_protein_sigma ~ std_normal();
    // log(prot_per_cell) ~ normal(log_protein_intercept + log_protein_slope .* prot_growth_rate, log_protein_sigma);
}

generated quantities {
    vector<lower=0>[N_bca_cal] bca_cal_ppc;
    vector<lower=0>[N_brad_cal] brad_cal_ppc;
    vector<lower=0>[N_biuret_cal] biuret_cal_ppc;
    vector<lower=0>[N_mem] mem_ppc;
    vector<lower=0>[N_peri] peri_ppc;
    vector<lower=0>[N_prot] prot_ppc;

    vector<lower=0>[N_growth] growth_rate_ppc;
    vector<lower=0>[N_flow] cell_count_ppc;
    vector<lower=0>[N_size] width_ppc;
    vector<lower=0>[N_size] length_ppc;
    vector<lower=0>[N_size] surface_area_ppc;
    vector<lower=0>[N_size] volume_ppc;
    vector<lower=0>[N_size] alpha_ppc;

    // Membrane protein ppc
    for (i in 1:N_bca_cal) {
        bca_cal_ppc[i] = normal_rng(bca_intercept + bca_slope * bca_std_conc[i], bca_sigma);
    }

    for (i in 1:N_mem) {
        mem_ppc[i] = bca_intercept + mem_od600nm[i] * exp(normal_rng(log(bca_slope * mem_per_biomass_mu[mem_idx[i]] ./ mem_conv_factor[i]), log_od562nm_mem_sigma[mem_idx[i]]));
    }

    // Periplasmic protein ppc
    for (i in 1:N_brad_cal) {
        brad_cal_ppc[i] = normal_rng(brad_intercept + brad_slope * brad_std_conc[i], brad_sigma);
    }
    for (i in 1:N_peri) {
        peri_ppc[i] = brad_intercept + (normal_rng(peri_per_biomass_mu[peri_idx[i]], peri_per_biomass_sigma[peri_idx[i]]) * peri_od600nm[i] * brad_slope  / peri_conv_factor[i]);
    }

    // Total protein ppc
    for (i in 1:N_biuret_cal) { 
        biuret_cal_ppc[i] = normal_rng(biuret_intercept + biuret_slope * biuret_std_conc[i], biuret_sigma);
    }
    for (i in 1:N_prot) {
        prot_ppc[i] = biuret_intercept + prot_od600nm[i] * exp(normal_rng(log(biuret_slope * prot_per_biomass_mu[prot_idx[i]] ./ prot_conv_factor[i]), log_od555nm_prot_sigma[prot_idx[i]]));
    }

    // Growth rate ppc
    for (i in 1:N_growth) {
        growth_rate_ppc[i] = normal_rng(lam_mu[growth_idx[i]], lam_sigma[growth_idx[i]]);
    }

    // Flow ppc
    for (i in 1:N_flow) {
        cell_count_ppc[i] = 1E9 * normal_rng(cells_per_biomass_mu[flow_idx[i]], cells_per_biomass_sigma[flow_idx[i]]);
    }

    // Cell shape ppcs
    for (i in 1:N_size) {
        width_ppc[i] = normal_rng(width_mu[size_idx[i]], width_sigma[size_idx[i]]);
        length_ppc[i] = normal_rng(length_mu[size_idx[i]], length_sigma[size_idx][i]);
        alpha_ppc[i] = normal_rng(alpha_mu[size_idx[i]], alpha_sigma[size_idx][i]);
        surface_area_ppc[i] = normal_rng(surface_area_mu[size_idx[i]], surface_area_sigma[size_idx][i]);
        volume_ppc[i] = normal_rng(volume_mu[size_idx[i]], volume_sigma[size_idx][i]); 
    }
}