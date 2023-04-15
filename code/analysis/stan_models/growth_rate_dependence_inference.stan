data { 
    int<lower=1> N_size;
    int<lower=1> N_flow;
    int<lower=1> N_prot;
    int<lower=1> N_biomass;
    int<lower=1> N_mass_spec;
    int<lower=0> const_phi_mem; // Factor for deciding which model to consider . If 1, phi_mem is taken to be constant. if 0, rho_mem is taken to be constant.    
    // Mass spec measurements
    vector<lower=0, upper=1>[N_mass_spec] phi_peri;
    vector<lower=0, upper=1>[N_mass_spec] phi_mem;
    vector<lower=0>[N_mass_spec] mass_spec_lambda;

    // Size Measurements
    vector<lower=0>[N_size] size_lambda; // Growth rates for size measurement
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] length;

    // Cells per biomass
    vector<lower=0>[N_flow] flow_meas; // Cells detected per biomass
    vector<lower=0>[N_flow] flow_lambda; // Inferred growth rate per condition 

    // Total protein
    vector<lower=0>[N_prot] prot_meas; // Protein mass per OD ml
    vector<lower=0>[N_prot] prot_lambda; // Condition growth rates

    // Total biomass
    vector<lower=0>[N_biomass] biomass;

}

transformed data {
    vector[N_biomass] biomass_centered = (biomass - mean(biomass)) ./ sd(biomass);
    // vector[N_mass_spec] rel_phi = phi_per ./ phi_mem;

}

parameters { 
    real<lower=1> alpha;
    real<lower=0> m_peri;
    real<lower=0> phi_mem_sigma;
    real<lower=0> phi_peri_sigma; 
    real<lower=0> w_min;
    real<lower=0> k_w;
    real<lower=0> w_sigma;
    real<lower=0> ell_sigma;
    real<lower=0> m_min;
    real<lower=0> k_m;
    real<lower=0> m_sigma;

    // Flow measurements 
    real<lower=0> flow_prefactor;
    real<lower=0> flow_sigma;

    // Biomass
    real biomass_centered_mu;
    real<lower=0> biomass_sigma; 
    array[const_phi_mem] real<lower=0, upper=1> phi_mem_mu; 
    array[1 - const_phi_mem] real<lower=0> rho_mem;
}

transformed parameters { 
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);
    real<lower=0> flow_slope = flow_prefactor * biomass_mu;
    vector<lower=0>[N_flow] flow_width = w_min + k_w .* flow_lambda;
    vector<lower=0>[N_flow] flow_length = alpha .* (flow_width); 
    vector<lower=0>[N_flow] flow_volume = (pi()/12) * flow_width.^2 .* (3 * flow_length - flow_width);
    vector<lower=0>[N_flow] flow_mu = (flow_slope .* flow_volume); 
    vector<lower=0>[N_prot] prot_width = w_min + k_w .* prot_lambda;
    vector<lower=0>[N_prot] prot_length = alpha .* prot_width; 
    vector<lower=0>[N_prot] prot_volume = (pi()/12) .* prot_width.^2 .* (3 .* prot_length - prot_width);
    vector<lower=0>[N_mass_spec] ms_width = w_min + k_w .* mass_spec_lambda;
    vector<lower=0>[N_mass_spec] ms_length = alpha .* ms_width;
    vector<lower=0>[N_mass_spec] ms_volume = (pi() / 12) .* ms_width^2 .* (3 .* ms_length - ms_width);

}

model { 
    alpha ~ normal(1, 2);
    phi_peri_sigma ~ std_normal();
    phi_mem_sigma ~ std_normal();
    m_peri ~ normal(0, 10);
    flow_prefactor ~ std_normal();    
    w_min ~ std_normal();
    m_min ~ normal(0, 200);
    k_w ~ std_normal();
    k_m ~ normal(100, 50);
    biomass_centered_mu ~ std_normal();
    biomass_sigma ~ std_normal();
    w_sigma ~ std_normal();
    m_sigma ~ std_normal();

    if (const_phi_mem) {
        phi_mem_mu ~ beta(2, 10);
        phi_mem ~ normal(phi_mem_mu[1], phi_mem_sigma);
    }

    else { 
        rho_mem ~ normal(0, 10);
        phi_mem ~ normal(rho_mem[1] * pi() * alpha .* ms_width.^2 ./ (ms_volume .* k_m), phi_mem_sigma);

    }

    width ~ normal(w_min + k_w * size_lambda, w_sigma);
    length ~ normal(alpha .* (w_min + k_w * size_lambda), ell_sigma);//ell_min + k_ell * size_lambda, ell_sigma);
    prot_meas ~ normal((m_min + k_m .* prot_volume)./(flow_slope .* prot_volume), m_sigma);
    flow_meas./1E9 ~ normal(1/flow_mu, flow_sigma);
    biomass_centered ~ normal(biomass_centered_mu, biomass_sigma);
    phi_peri ~ normal(m_peri ./ (ms_volume .* k_m), phi_peri_sigma);
    // rel_phi ~ normal

}

generated quantities { 
    vector[N_size] width_rep;
    vector[N_size] length_rep;
    vector[N_size] volume_rep;
    vector[N_prot] prot_rep;
    vector[N_mass_spec] rho_mem_rep;
    vector[N_mass_spec] m_peri_rep;
    vector[N_mass_spec] phi_peri_rep;
    vector[N_mass_spec] phi_mem_rep;
    vector[N_mass_spec] rel_phi_rep;

    for (i in 1:N_size) {
        width_rep[i] = normal_rng(w_min + k_w * size_lambda[i], w_sigma);
        length_rep[i] = normal_rng(alpha * (w_min + k_w * size_lambda[i]), ell_sigma);
        volume_rep[i] = (pi()/12) * width_rep[i]^2 * (3 * length_rep[i] - width_rep[i]);
    }
    for (i in 1:N_prot) {
        prot_rep[i] = normal_rng((m_min + k_m .* prot_volume[i]) / (flow_slope * prot_volume[i]), m_sigma);
    }

    for (i in 1:N_mass_spec) {
        m_peri_rep[i] = phi_peri[i] * k_m * ms_volume[i];
        rho_mem_rep[i] = phi_mem[i] * k_m * ms_volume[i] / (pi() * ms_width[i] * ms_length[i]);
        phi_peri_rep[i] = normal_rng(m_peri / (ms_volume[i] * k_m), phi_peri_sigma);
        if (const_phi_mem) {
            phi_mem_rep[i] = normal_rng(phi_mem_mu[1], phi_mem_sigma);
        }
        else {
            phi_mem_rep[i] = normal_rng(rho_mem[1] * pi() * alpha * ms_width[i]^2 / (ms_volume[i] * k_m), phi_mem_sigma);
        }
        // rel_phi_rep[i] =  
    }
}