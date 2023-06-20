 data { 
    // Dimensional information
    int<lower=1> N_ms; // Number of mass spectrometry measurements
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_ppc;  // Number of elements desired in ppcs. 

    // Observed data
    vector<lower=0>[N_ms] rho_mem; // Areal membrane protein density
    vector<lower=0>[N_ms] phi_mem; // Allocation towards membrane protein
    vector<lower=0>[N_ms] m_peri; // Mass of periplasmic protein per cell
    vector<lower=1>[N_size] alpha;  // Length-to-width aspect ratio

    // Variable for ppc simulation
    vector<lower=0>[N_ppc] phiRb; // Allocation towards ribosomes

    // Constants for calculations
    real<lower=0> rho_biomass; // Drymass density in units of fg / fL
    real<lower=0> beta_rp; // Conversion factor from R/P to phiRb
    real<lower=0> delta; // Periplasmic width
}

parameters {
    real<lower=0> rho_mem_mu;
    real<lower=0> rho_mem_sigma;
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=1> alpha_mu;
    real<lower=0> alpha_sigma;
}

model {
    rho_mem_mu ~ std_normal();
    rho_mem_sigma ~ normal(0, 0.1);
    phi_mem_mu ~ std_normal();
    phi_mem_sigma ~ normal(0, 0.1);
    alpha_mu ~ std_normal();
    alpha_sigma ~ normal(0, 0.1);
    m_peri_mu ~ std_normal();
    m_peri_sigma ~ normal(0, 0.1);
    rho_mem ~ normal(rho_mem_mu, rho_mem_sigma);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    alpha ~ normal(alpha_mu, alpha_sigma);
    m_peri ~ normal(m_peri_mu, m_peri_sigma);
}

generated quantities {
   // Scalar PPC declaration
   real<lower=0> rho_mem_rep; 
   real<lower=0> m_peri_rep;
   real<lower=1> alpha_rep;
   real<lower=0> phi_mem_rep;

   // Vector PPC declaration
   vector<lower=0>[N_ppc] kappa; // Scaling constant for width
   vector<lower=0>[N_ppc] width_rep; // Average cell width
   vector<lower=0>[N_ppc] ell_rep; // Average cell length
   vector<lower=0>[N_ppc] vol_rep; // Average cell volume
   vector<lower=0>[N_ppc] rho_peri; // Average periplasmic protein density 

   // Scalar PPC generation
   rho_mem_rep = normal_rng(rho_mem_mu, rho_mem_sigma);
   m_peri_rep = normal_rng(m_peri_mu, m_peri_sigma);     
   alpha_rep = normal_rng(alpha_mu, alpha_sigma);
   phi_mem_rep = normal_rng(phi_mem_mu, phi_mem_sigma);
   for (i in 1:N_ppc) { 
        kappa[i] = (24 * alpha_mu / (3 * alpha_mu - 1)) * (rho_mem_mu / rho_biomass);
        width_rep[i] = kappa[i] * (1 + phiRb[i] / beta_rp) / phi_mem_mu;
        rho_peri[i] = m_peri_mu / (pi() * alpha_mu * width_rep[i]^2 * delta);
        ell_rep[i] = alpha_mu * width_rep[i];
        vol_rep[i] = (pi() / 12) * width_rep[i]^3 * (3 * alpha_mu  - 1);
   }
}