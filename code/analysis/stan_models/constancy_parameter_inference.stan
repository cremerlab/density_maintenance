data {
    // Parameters
    int<lower=1> N_mass_spec; // Total number of mass spec conditions
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_flow; // Number of flow ctyometry measurements
    int<lower=1> N_prot; // Number of total protein measurements
    int<lower=1> N_biomass; // Number of biomass measurements   

    // Mass Spec Measurements
    vector<lower=0>[N_mass_spec] phi_peri; // allocation towards periplasmic proteins
    vector<lower=0>[N_mass_spec] phi_mem; // allocation towards periplasmic proteins
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

tranformed data {
    vector<lower=0>[N_biomass] biomass_centered = (biomass - mean(biomass)) ./ sd(biomass);
}


parameters { 
    // Model parameters
    real<lower=0> m_peri; // Per cell periplasmic biomass
    real<lower=0> rho_mem; // Areal density of membrane protein mass
    real<lower=0> phi_peri_sigma;
    real<lower=0> phi_memb_sigma;
    // Growth rate dependence parameters
    real<lower=0> w_min;
    real k_w;
    real<lower=0> w_sigma;
    real<lower=0> ell_min;
    real k_ell;
    real<lower=0> ell_sigma;
    real<lower=0> M_min;
    real k_p;
    real<lower=0> M_sigma;

    // Flow measurements 
    real<lower=0> flow_prefactor;
    real<lower=0> flow_sigma;

    // Biomass
    real<lower=0> biomass_centered_mu;
    real<lower=0> biomass_sigma;
}

transformed parameters {
    real<lower=0> alpha = mean((ell_min + k_ell .* size_lambda) ./ (w_min + k_w .* size_lambda));
    real<lower=0> biomass_mu = biomass_centered_mu * sd(biomass) + mean(biomass);

    real<lower=0> flow_slope = flow_prefactor * biomass_mu;
    vector<lower=0>[N_flow] flow_volume = (pi()/12) * (w_min + k_w .* flow_lambda).^3 .* (3 * alpha - 1);
    vector<lower=0>[N_flow] flow_mu = (flow_prefactor * biomass_mu * flow_volume);

    
    vector<lower=0>[N_mass_spec] ms_volume = (pi()/12) * (w_min + k_w .* mass_spec_lambda).^3 .* (3 * alpha - 1);
    vector<lower=0>[N_mass_spec] M_prot_ms = M_min + k_p .* mass_spec_lambda;
    vector<lower=0>[N_mass_spec] ms_width = w_min + k_w .* mass_spec_lambda;
    vector<lower=0>[N_prot] N_cells_ms = 1E9 ./ (flow_slope .* mass_spec_volume);
}

model {
    // Define the priors
    m_peri ~ normal(0, 20);
    rho_mem ~ normal(0, 20);
    w_min ~ std_normal(); 
    k_w ~ std_normal();
    ell_min ~ normal(0, 2);
    k_ell ~ std_normal();
    M_min ~ normal(300, 100);
    k_m ~ normal(0, 100);


    // Lambda-dependence likelihoods

    // Biomass likelihood

    // Cells per biomass likelihood

    // Allocation likelihood
    phi_peri ~ normal(m_peri * N_cells_ms ./ M_prot_ms, phi_peri_sigma);
    phi_memb ~ normal((rho_mem * pi() * alpha * ms_width.^2) ./ (M_prot_ms ./ N_cells_ms));
}