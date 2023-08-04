 data { 
    // Dimensional information
    int<lower=1> N_ms; // Number of mass spectrometry measurements
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_ppc;  // Number of elements desired in ppcs. 
    int<lower=1> N_biomass;
    int<lower=0, upper=1> linear_alpha;

    // Observed data
    vector<lower=0>[N_ms] rho_mem; // Areal membrane protein density
    vector<lower=0>[N_ms] phi_mem; // Allocation towards membrane protein
    vector<lower=0>[N_ms] m_peri; // Mass of periplasmic protein per cell
    vector<lower=1>[N_size] alpha;  // Length-to-width aspect ratio
    vector<lower=1>[N_biomass] rho_biomass;
    vector<lower=0>[N_size] size_lam;

    // Variable for ppc simulation
    vector<lower=0>[N_ppc] phiRb; // Allocation towards ribosomes
    vector<lower=0>[N_ppc] w_range;
    vector<lower=0>[N_ppc] lam;


    // Constants for calculations
    real<lower=0> beta_rp; // Conversion factor from R/P to phiRb
    real<lower=0> delta; // Periplasmic width
    real<lower=0> SA_prefactor;
}


parameters {
    real<lower=0> rho_mem_mu;
    real<lower=0> rho_mem_sigma;
    real<lower=0> rho_biomass_mu;
    real<lower=0> rho_biomass_sigma;
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    array[1 - linear_alpha] real<lower=1> alpha_mu;
    array[linear_alpha] real<lower=0> alpha_int;
    array[linear_alpha] real<lower=0> alpha_slope;
    real<lower=0> alpha_sigma;
}


model {
    rho_mem_mu ~ std_normal();
    rho_mem_sigma ~ normal(0, 0.1);
    rho_biomass_mu ~ normal(300, 100);
    rho_biomass_sigma ~ normal(0, 10);
    phi_mem_mu ~ std_normal();
    phi_mem_sigma ~ normal(0, 0.1);
    alpha_sigma ~ normal(0, 0.1);
    m_peri_mu ~ std_normal();
    m_peri_sigma ~ normal(0, 0.1);

    rho_mem ~ normal(rho_mem_mu, rho_mem_sigma);
    rho_biomass ~ normal(rho_biomass_mu, rho_biomass_sigma);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    if (linear_alpha) { 
        alpha_slope ~ std_normal();
        alpha_int ~ std_normal();
        alpha ~ normal(alpha_int[1] + alpha_slope[1] .* size_lam, alpha_sigma);
    }
    else { 
        alpha_mu ~ std_normal();
        alpha ~ normal(alpha_mu[1], alpha_sigma);
    }

    m_peri ~ normal(m_peri_mu, m_peri_sigma);
}

generated quantities {
   // Scalar PPC declaration
   real<lower=0> rho_mem_rep; 
   real<lower=0> rho_biomass_rep;
   real<lower=0> m_peri_rep;
   vector<lower=0>[N_ppc] alpha_rep;
   vector<lower=0>[N_ppc] alpha_mu_0;
   real<lower=0> phi_mem_rep;

   // Vector PPC declaration
   vector<lower=0>[N_ppc] kappa; // Scaling constant for width
   vector<lower=0>[N_ppc] width_rep; // Average cell width
   vector[N_ppc] rib_rep; // Average cell width
   vector<lower=0>[N_ppc] ell_rep; // Average cell length
   vector<lower=0>[N_ppc] vol_rep; // Average cell volume
   vector<lower=0>[N_ppc] rho_peri; // Average periplasmic protein density 

   // Scalar PPC generation
   rho_mem_rep = normal_rng(rho_mem_mu, rho_mem_sigma);
   m_peri_rep = normal_rng(m_peri_mu, m_peri_sigma);     
   phi_mem_rep = normal_rng(phi_mem_mu, phi_mem_sigma);
   rho_biomass_rep = normal_rng(rho_biomass_mu, rho_biomass_sigma);
   if (linear_alpha) {
       for (i in 1:N_ppc) {
        alpha_rep[i] = normal_rng(alpha_int[1] + alpha_slope[1] * lam[i], alpha_sigma);
        alpha_mu_0[i] = alpha_int[1] + alpha_slope[1] * lam[i];
            }
        }
   else {
        alpha_rep[1] = normal_rng(alpha_mu[1], alpha_sigma);
        alpha_mu_0[1] = alpha_mu[1];
        for (i in 2:N_ppc) {
            alpha_rep[i] = alpha_rep[1];
            alpha_mu_0[i] = alpha_mu_0[1];
        }
   }

   for (i in 1:N_ppc) { 
        kappa[i] = (SA_prefactor * alpha_mu_0[i] / (3 * alpha_mu_0[i] - 1)) * (rho_mem_mu / rho_biomass_mu);
        width_rep[i] = kappa[i] * (1 + phiRb[i] / beta_rp) / phi_mem_mu;
        rib_rep[i] = (w_range[i] * phi_mem_mu / kappa[i]) - 1; 
        rho_peri[i] = m_peri_mu / (pi() * alpha_mu_0[i] * width_rep[i]^2 * delta);
        ell_rep[i] = alpha_mu_0[i] * width_rep[i];
        vol_rep[i] = (pi() / 12) * width_rep[i]^3 * (3 * alpha_mu_0[i]  - 1);
   }
}