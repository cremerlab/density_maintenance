data {
    // Dimensional information
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_mass_spec;
    int<lower=1> N_sim;
    real<lower=0> delta;
    vector<lower=0>[N_sim] lam_sim;
    vector<lower=0>[N_sim] width_sim;
    int<lower=0> const_phi_mem;

    // Observed data
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] volumes;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_prot] rho_prot_meas;
    vector<lower=0>[N_mass_spec] phi_mem;
    vector<lower=0>[N_mass_spec] phi_peri;
    vector<lower=0>[N_mass_spec] rho_mem_meas;
    vector<lower=0>[N_mass_spec] rho_prot_meas_ms;
    vector<lower=0>[N_mass_spec] m_peri_meas;
    vector<lower=0>[N_mass_spec] ms_lam;
}

parameters { 
    real<lower=0> w_min;
    real<lower=0> w_slope;
    real<lower=0> w_sigma;
    real<lower=0> ell_sigma;
    real<lower=0> vol_sigma;
    real<lower=0> prot_sigma;
    real<lower=1> alpha;
    real<lower=0> rho_prot_min;
    real<upper=0> rho_prot_slope;
    real<lower=0> rho_prot_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=0> phi_peri_sigma;
    array[const_phi_mem] real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    array[1 - const_phi_mem] real<lower=0> rho_mem_mu;  
    array[1 - const_phi_mem] real<lower=0> rho_mem_sigma;  
}

transformed parameters {
    vector<lower=0>[N_prot] rho_prot = rho_prot_min + rho_prot_slope .* prot_lam;
    vector<lower=0>[N_mass_spec] rho_prot_ms = rho_prot_min + rho_prot_slope .* ms_lam;
}

model {
   w_min ~ normal(0, 0.1);
   w_slope ~ std_normal();
   w_sigma ~ normal(0, 0.01);
   vol_sigma ~ normal(0, 0.1);
   ell_sigma ~ normal(0, 0.1);
   rho_prot_min ~ normal(200, 100);
   rho_prot_slope ~ normal(0, 100);
   rho_prot_sigma ~ normal(0, 10);
   phi_peri_sigma ~ normal(0, 0.1);
   phi_mem_sigma ~ normal(0, 0.1);
   alpha ~ normal(1, 3);
   m_peri_mu ~ normal(0, 100);
   if (const_phi_mem) {
        phi_mem_mu ~ beta(2, 10);
        phi_mem ~ normal(phi_mem_mu[1], phi_mem_sigma);
   }
   else {
        rho_mem_mu ~ std_normal();
        rho_mem_sigma ~ std_normal();
        rho_mem_meas ~ normal(rho_mem_mu[1], rho_mem_sigma[1]);
        phi_mem ~ normal(rho_mem_mu[1] * pi() * alpha .* (w_min + w_slope * ms_lam).^2 ./ ((pi()/12) * (w_min + w_slope * ms_lam).^3 * (3 * alpha - 1) .* (rho_prot_min + rho_prot_slope .* ms_lam)), phi_mem_sigma);
   }

   // Likelihoods for size measurements
   widths ~ normal(w_min + w_slope .* size_lam, w_sigma);
   lengths ~ normal(alpha .* (w_min + w_slope .* size_lam), ell_sigma);
   log(volumes) ~ normal(log((pi()/12) * (w_min + w_slope .* size_lam).^3 * (3 * alpha - 1)), vol_sigma);

   // Likelihoods for protein measurements
   rho_prot_meas ~ normal(rho_prot_min + rho_prot_slope .* prot_lam, rho_prot_sigma);
   log(prot_per_cell) ~ normal(log(rho_prot * (pi()/12) .* (w_min + w_slope .* prot_lam).^3 * (3 * alpha - 1)), prot_sigma);

   // Likelihoods based on protein measurements
   phi_peri ~ normal(m_peri_mu / (rho_prot_ms * (pi()/12) .* (w_min + w_slope .* ms_lam).^3 * (3 * alpha - 1)), phi_peri_sigma);
   m_peri_meas ~ normal(m_peri_mu, m_peri_sigma);

}

generated quantities {
    // PPC for fit quantities
    vector[N_sim] alpha_rep;
    vector[N_sim] alpha_sim;
    vector[N_sim] w_rep;
    vector[N_sim] ell_rep;
    vector[N_sim] vol_rep;
    vector[N_sim] rho_prot_rep;
    vector[N_sim] rho_prot_sim;
    vector[N_sim] prot_per_cell_rep;
    vector[N_sim] phi_mem_rep;
    vector[N_sim] phi_peri_rep;
    vector[N_sim] rho_peri_rep;
    vector[N_sim] rho_mem_rep;
    vector[N_sim] m_peri_rep;
    vector[N_sim] w_sim;
    vector[N_sim] ell_sim;
    vector[N_sim] vol_sim;
    vector[N_sim] prot_per_cell_sim;
    vector[N_sim] phi_mem_sim;
    vector[N_sim] phi_peri_sim;
    vector[N_sim] rho_peri_sim;
    vector[N_sim] rho_mem_sim;
    vector[N_sim] m_peri_sim;
    vector[N_sim] rel_phi_rep;
    vector[N_sim] rel_phi_sim;

    for (i in 1:N_sim) {
        alpha_rep[i] = alpha; 
        alpha_sim[i] = alpha;
        m_peri_sim[i] = m_peri_mu;
        m_peri_rep[i] = normal_rng(m_peri_mu, m_peri_sigma);
        rho_prot_rep[i] = normal_rng(rho_prot_min + rho_prot_slope * lam_sim[i], rho_prot_sigma);
        rho_prot_sim[i] = rho_prot_min + rho_prot_slope * lam_sim[i];
        w_sim[i] = w_min + w_slope * lam_sim[i];
        w_rep[i] = normal_rng(w_min + w_slope * lam_sim[i], w_sigma);
        ell_sim[i] = alpha * w_sim[i];
        ell_rep[i] = normal_rng(alpha * w_rep[i], ell_sigma);
        vol_rep[i] = exp(normal_rng(log((pi()/12) * (w_min + w_slope * lam_sim[i])^3 * (3 * alpha - 1)), vol_sigma));
        vol_sim[i] = (pi()/12) * w_sim[i]^3 * (3 * alpha - 1);
        prot_per_cell_rep[i] = normal_rng(rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1), prot_sigma);    
        prot_per_cell_sim[i] = rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha - 1);    
        if (const_phi_mem) {
            phi_mem_rep[i] = normal_rng(phi_mem_mu[1], phi_mem_sigma);
            phi_mem_sim[i] = phi_mem_mu[1];
            rho_mem_rep[i] = phi_mem_rep[i] * rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1) / (2 * pi() * alpha * w_rep[i]^2);
            rho_mem_sim[i] = phi_mem_sim[i] * rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha - 1) / (2 * pi() * alpha * w_sim[i]^2);

        }
        else {
            phi_mem_rep[i] = normal_rng(rho_mem_mu[1] * (2 * pi() * alpha * w_rep[i]^2) / (vol_rep[i] * rho_prot_rep[i]), phi_mem_sigma);
            phi_mem_sim[i] =  rho_mem_mu[1] * 2 * pi() * alpha * w_sim[i]^2 / (vol_sim[i] * rho_prot_sim[i]);
            rho_mem_rep[i] = normal_rng(rho_mem_mu[1], rho_mem_sigma[1]);
            rho_mem_sim[i] = rho_mem_mu[1];
        }
        phi_peri_rep[i] = normal_rng(m_peri_rep[i] / (rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1)), phi_peri_sigma);
        phi_peri_sim[i] = m_peri_sim[i] / (rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha - 1));
        rho_peri_rep[i] = m_peri_mu / (pi() * alpha * delta * w_rep[i]^2);
        rho_peri_sim[i] = m_peri_mu / (pi() * alpha * delta * w_sim[i]^2);

        if (const_phi_mem) {
            rel_phi_rep[i] = 12 * m_peri_rep[i] / (pi() * (3 * alpha - 1) * rho_prot_rep[i] * width_sim[i]^3 * phi_mem_rep[i]);
            rel_phi_sim[i] = 12 * m_peri_sim[i] / (pi() * (3 * alpha - 1) * rho_prot_sim[i] * width_sim[i]^3 * phi_mem_sim[i]);

        } 
        else {
            rel_phi_rep[i] = (m_peri_rep[i]) / (pi() * alpha * rho_mem_rep[i] * width_sim[i]^2);
            rel_phi_sim[i] = (m_peri_sim[i]) / (pi() * alpha * rho_mem_sim[i] * width_sim[i]^2); 

        }
        
    }
}

