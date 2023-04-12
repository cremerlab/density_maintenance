data {
    // Parameters
    int<lower=1> N_mass_spec; // Total number of mass spec conditions
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_flow; // Number of flow ctyometry measurements
    int<lower=1> N_prot; // Number of total protein measurements
    int<lower=1> N_biomass; // Number of biomass measurements   

    // Mass Spec Measurements
    vector<lower=0>[N_mass_spec] phi_peri; // allocation towards periplasmic proteins
    vector<lower=0>[N_mass_spec] phi_memb; // allocation towards periplasmic proteins
    vector<lower=0>[N_mass_spec] mass_spec_lambda; // Growth rates for mass spec conditions

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
    vector[N_mass_spec] rel_phi = phi_peri ./ phi_memb;
}


parameters { 
    // Model parameters
    real<lower=0> m_peri; // Per cell periplasmic biomass
    real<lower=0> rho_mem; // Areal density of membrane protein mass
    real<lower=0> phi_peri_sigma;
    real<lower=0> phi_memb_sigma;
    real<lower=0> rel_sigma;

    // Growth rate dependence parameters
    real<lower=0> w_min;
    real<lower=0> ell_min;
    real<lower=0> M_min;
    real k_w;
    real k_ell;
    real<lower=0> k_m;
    real<lower=0> w_sigma;
    real<lower=0> ell_sigma;
    real<lower=0> M_sigma;

    // Flow measurements 
    real<lower=0> flow_prefactor;
    real<lower=0> flow_sigma;

    // Biomass
    real biomass_centered_mu;
    real<lower=0> biomass_sigma;
}

transformed parameters {
    real<lower=0> alpha = ell_min / w_min; //mean((ell_min + k_ell .* size_lambda) ./ (w_min + k_w .* size_lambda));
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);
    real<lower=0> flow_slope = flow_prefactor * biomass_mu;

    vector<lower=0>[N_flow] flow_width = w_min + k_w .* flow_lambda;
    vector<lower=0>[N_flow] flow_length = ell_min + k_ell .* flow_lambda;
    vector<lower=0>[N_flow] flow_volume = (pi()/12) * flow_width.^3 .* (3 * flow_length - flow_width);
    vector<lower=0>[N_flow] flow_mu = (flow_slope .* flow_volume); 
    // vector<lower=0>[N_mass_spec] M_prot_ms = (M_min + k_m .* mass_spec_lambda) * 1E9;
    vector<lower=0>[N_mass_spec] ms_width = w_min + k_w .* mass_spec_lambda;
    vector<lower=0>[N_mass_spec] ms_length = ell_min + k_ell .* mass_spec_lambda;
    vector<lower=0>[N_mass_spec] ms_volume = (pi()/12) * ms_width.^2 .* (3 * ms_length - ms_width); 
    vector<lower=0>[N_prot] prot_width = w_min + k_w .* prot_lambda;
    vector<lower=0>[N_prot] prot_length = ell_min + k_ell .* prot_lambda;
    vector<lower=0>[N_prot] prot_volume = (pi()/12) .* prot_width.^2 .* (3 .* prot_length - prot_width);
    // vector<lower=0>[N_prot] prot_cells = 1E9 ./ (flow_slope * prot_volume);
    // vector<lower=0>[N_mass_spec] N_cells_ms = 1E9 ./ (flow_slope .* ms_volume);
    real kappa = m_peri / (alpha * rho_mem);
}

model {
    // Define the model priors
    m_peri ~ normal(0, 20);
    rho_mem ~ normal(0, 20);
    phi_peri_sigma ~ normal(0, 0.1);
    phi_memb_sigma ~ normal(0, 0.1);
    rel_sigma ~ std_normal();
    flow_prefactor ~ std_normal();
    biomass_centered_mu ~ std_normal();
    biomass_sigma ~ std_normal();

    // Define the growth rate dependence priors
    w_min ~ std_normal(); 
    k_w ~ std_normal();
    w_sigma ~ std_normal();
    ell_min ~ normal(0, 2);
    k_ell ~ std_normal();
    ell_sigma ~ std_normal();
    M_min ~ normal(300, 100);
    k_m ~ normal(100, 50);
    M_sigma ~ normal(0, 10);

    // Lambda-dependence likelihoods
    width ~ normal(w_min + k_w .* size_lambda, w_sigma);
    length ~ normal(ell_min + k_ell .* size_lambda, ell_sigma);
    prot_meas ~ normal((M_min + k_m .* prot_volume)./(flow_slope .* prot_volume), M_sigma);

    // Biomass likelihood
    biomass_centered ~ normal(biomass_centered_mu, biomass_sigma);

    // Cells per biomass likelihood
    flow_meas/1E9 ~ normal(1/flow_mu, flow_sigma);

    // Allocation likelihood
    phi_peri ~ normal(m_peri / (ms_volume * k_m), phi_peri_sigma);
    phi_memb ~ normal((rho_mem .* pi() .* ms_length .* ms_width)./(ms_volume .* k_m), phi_memb_sigma);
    rel_phi ~ normal(kappa / ms_width.^2, rel_sigma);
}

generated quantities {
    // Model PPC
    // vector[N_mass_spec] phi_peri_rep;
    // vector[N_mass_spec] phi_memb_rep;
    vector[N_mass_spec] rel_phi_rep;
    vector[N_mass_spec] rho_mem_empirical;
    vector[N_mass_spec] m_peri_empirical;
    // vector[N_mass_spec] prot


    for (i in 1:N_mass_spec) {
        // phi_peri_rep[i] = normal_rng(m_peri * N_cells_ms[i] ./ (M_prot_ms[i] * 1E9), phi_peri_sigma);
        // phi_memb_rep[i] = normal_rng((rho_mem * pi() * alpha * ms_width[i]^2 * N_cells_ms[i]) ./ (M_prot_ms[i] * 1E9),phi_memb_sigma);
        
        rel_phi_rep[i] = normal_rng(kappa / ms_width[i]^2, rel_sigma);
        rho_mem_empirical[i] = (phi_memb[i] * ms_volume[i] * k_m)/(pi() * ms_length[i] * ms_width[i]);
        m_peri_empirical[i] = phi_peri[i] * ms_volume[i] * k_m;
    }

}