//
// Literature Based Model Inference 
// ============================================================================ 
// This model defines a Bayesian inferential model to estimate different statistical
// models of protein density and cell size control. 
//
// There are two integer (max, 1) flags that differentiate the sinferential models.
//  `smk` - if 1, models are inferred assuming an SMK growth law for the volume. 
//   width is computed assuming a spherocylinder with an aspect ratio alpha. 
//   If 0, models are inferred assuming a linear depenendce of the width on the 
//   growth rate. Volume is computed assuming a spherocylinder with an asptect 
//   ratio alpha.
//
//  `const_phi_mem` - If 1, the model which assumes a constant allocation towards
//   membrane proteins is implemented. If 0, a model which assumes a constant 
//   membrane protein density is applied.
//
//
data {
    // Flags to differential inferrential modes
    int<lower=0> smk;
    int<lower=0> const_phi_mem;

    // Dimensional information
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_mass_spec;
    int<lower=1> N_sim;
    real<lower=0> delta;
    vector<lower=0>[N_sim] lam_sim;
    vector<lower=0>[N_sim] width_sim;

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

transformed data {
    vector<lower=1>[N_size] aspect_ratio = lengths ./ widths;
}

parameters { 
    array[1 - smk] real<lower=0> w_min;
    array[1 - smk] real<lower=0> w_slope;
    array[smk] real<lower=0> tau;
    array[smk] real<lower=0> C_0;
    real<lower=0> w_sigma;
    real<lower=0> ell_sigma;
    real<lower=0> vol_sigma;
    real<lower=0> prot_sigma;
    real<lower=1> alpha;
    real<lower=0> alpha_sigma;
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
   if (smk) {
        tau ~ std_normal();
        C_0 ~ std_normal();
    }
   else { 
        w_min ~ normal(0, 0.1);
        w_slope ~ std_normal();
        w_sigma ~ normal(0, 0.01);
   }

   vol_sigma ~ normal(0, 0.1);
   ell_sigma ~ normal(0, 0.1);
   rho_prot_min ~ normal(200, 100);
   rho_prot_slope ~ normal(0, 100);
   rho_prot_sigma ~ normal(0, 10);
   phi_peri_sigma ~ normal(0, 0.1);
   phi_mem_sigma ~ normal(0, 0.1);
   alpha ~ normal(1, 3);
   alpha_sigma ~ std_normal();
   m_peri_mu ~ normal(0, 100);
   if (const_phi_mem) {
        phi_mem_mu ~ beta(2, 10);
        phi_mem ~ normal(phi_mem_mu[1], phi_mem_sigma);
   }
   else {
        rho_mem_mu ~ std_normal();
        rho_mem_sigma ~ std_normal();
        rho_mem_meas ~ normal(rho_mem_mu[1], rho_mem_sigma[1]);
   }

   // Likelihoods for size measurements
   if (smk) { 
        widths ~ normal(cbrt(12 * C_0[1] .* exp(size_lam .* tau[1]) ./ (pi() * (3 * alpha - 1))), w_sigma);
        lengths ~ normal(alpha * cbrt(12 * C_0[1] .* exp(size_lam .* tau[1]) ./ (pi() * (3 * alpha - 1))), ell_sigma);
        log(volumes) ~ normal(log(C_0[1]) + tau[1] * size_lam, vol_sigma);  
        log(prot_per_cell) ~ normal(log(rho_prot * C_0[1] .* exp(tau[1] .* prot_lam)), prot_sigma);
        phi_peri ~ normal(m_peri_mu / (rho_prot_ms * C_0[1] .* exp(tau[1] .* ms_lam)), phi_peri_sigma);
        if (const_phi_mem == 0) {
            phi_mem ~ normal(rho_mem_mu[1] *  pi() * alpha .* cbrt(12 * C_0[1] .* exp(ms_lam .* tau[1]) ./ (pi() * (3 * alpha - 1))).^2 ./ ((rho_prot_min + rho_prot_slope .* ms_lam) .* C_0[1] .* exp(ms_lam * tau[1])), phi_mem_sigma);
        }
        else {
            phi_mem ~ normal(phi_mem_mu[1], phi_mem_sigma);
        }

   }

   else { 
        widths ~ normal(w_min[1] + w_slope[1] .* size_lam, w_sigma);
        lengths ~ normal(alpha .* (w_min[1] + w_slope[1] .* size_lam), ell_sigma);
        log(volumes) ~ normal(log((pi()/12) * (w_min[1] + w_slope[1] .* size_lam).^3 * (3 * alpha - 1)), vol_sigma);
        log(prot_per_cell) ~ normal(log(rho_prot * (pi()/12) .* (w_min[1] + w_slope[1] .* prot_lam).^3 * (3 * alpha - 1)), prot_sigma);
        phi_peri ~ normal(m_peri_mu / (rho_prot_ms * (pi()/12) .* (w_min[1] + w_slope[1] .* ms_lam).^3 * (3 * alpha - 1)), phi_peri_sigma);
        if (const_phi_mem == 0) {
            phi_mem ~ normal(rho_mem_mu[1] * pi() * alpha .* (w_min[1] + w_slope[1] * ms_lam).^2 ./ ((pi()/12) * (w_min[1] + w_slope[1] * ms_lam).^3 * (3 * alpha - 1) .* (rho_prot_min + rho_prot_slope .* ms_lam)), phi_mem_sigma);   
        }
        else {
            phi_mem ~ normal(phi_mem_mu[1], phi_mem_sigma);
        }
   }

   aspect_ratio ~ normal(alpha, alpha_sigma);
   
   // Likelihoods for protein measurements
   rho_prot_meas ~ normal(rho_prot_min + rho_prot_slope .* prot_lam, rho_prot_sigma);

   // Likelihoods based on protein measurements
   m_peri_meas ~ normal(m_peri_mu, m_peri_sigma);

}

generated quantities {
    // Vectors of simulated distribution of measurement
    vector[N_sim] alpha_rep;
    vector[N_sim] w_rep;
    vector[N_sim] ell_rep;
    vector[N_sim] vol_rep;
    vector[N_sim] rho_prot_rep;
    vector[N_sim] prot_per_cell_rep;
    vector[N_sim] phi_peri_rep;
    vector[N_sim] rho_peri_rep;
    vector[N_sim] rho_mem_rep;
    vector[N_sim] m_peri_rep;
    vector[N_sim] phi_mem_rep;
    vector[N_sim] rel_phi_rep;
    vector[N_sim] calc_lam_rep;
    vector[N_sim] rel_phi_rho_prot_rep;
    vector[N_sim] sav_rep;
    vector[N_sim] m_mem_rep;

    // Vectors for calculated distribution of true value
    vector[N_sim] alpha_sim;
    vector[N_sim] rho_prot_sim;
    vector[N_sim] w_sim;
    vector[N_sim] ell_sim;
    vector[N_sim] vol_sim;
    vector[N_sim] prot_per_cell_sim;
    vector[N_sim] phi_mem_sim;
    vector[N_sim] phi_peri_sim;
    vector[N_sim] rho_peri_sim;
    vector[N_sim] rho_mem_sim;
    vector[N_sim] m_peri_sim;
    vector[N_sim] rel_phi_sim;
    vector[N_sim] calc_lam_sim;
    vector[N_sim] rel_phi_rho_prot_sim;
    vector[N_sim] sav_sim;
    vector[N_sim] m_mem_sim;

    // Iterate through each simulation value
    for (i in 1:N_sim) {
        // Quantities for volume-scaling option and model independent parameters
        m_peri_sim[i] = m_peri_mu;
        rho_prot_sim[i] = rho_prot_min + rho_prot_slope * lam_sim[i];
        alpha_sim[i] = alpha;

        // Expected measurement distribution
        alpha_rep[i] = normal_rng(alpha, alpha_sigma);
        rho_prot_rep[i] = normal_rng(rho_prot_min + rho_prot_slope * lam_sim[i], rho_prot_sigma);
        m_peri_rep[i] = normal_rng(m_peri_mu, m_peri_sigma);

        // Quantities for model independent, but volume scaling *dependent* parameters.
        if (smk) { 
            w_sim[i] = cbrt(12 * C_0[1] * exp(tau[1] * lam_sim[i]) / (pi() * (3 * alpha_sim[i] - 1)));
            vol_sim[i] = C_0[1] * exp(tau[1] * lam_sim[i]);

            // Expected measurement distribution
            w_rep[i] = normal_rng(cbrt(12 * C_0[1] * exp(tau[1] * lam_sim[i]) / (pi() * (3 * alpha_rep[i] - 1))), w_sigma);
            vol_rep[i] = exp(normal_rng(log(vol_sim[i]), vol_sigma));
        }

        else { 
            w_sim[i] = w_min[1] + w_slope[1] * lam_sim[i];
            vol_sim[i] = (pi()/12) * w_sim[i]^3 * (3 * alpha_sim[i] - 1);

            // Expected measurement distribution
            w_rep[i] = normal_rng(w_sim[i], w_sigma);
            vol_rep[i] = exp(normal_rng(log((pi()/12) * w_rep[i]^3 * (3 * alpha_rep[i] - 1)), vol_sigma));

        }
        // Length, SAV, and protein quantities, conditioned on scaling and model-dependent parameters.
        ell_sim[i] = alpha_sim[i] * w_sim[i];
        ell_rep[i] = normal_rng(alpha_rep[i] * w_rep[i], ell_sigma);  
        sav_sim[i] = pi() * ell_sim[i] * w_sim[i] / vol_sim[i];
        sav_rep[i] = pi() * ell_rep[i] * w_rep[i] / vol_rep[i];

        prot_per_cell_rep[i] = normal_rng(rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1), prot_sigma);    
        prot_per_cell_sim[i] = rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha - 1);    

        // Quantities for model-dependent allocation parameters
        if (const_phi_mem) {
            phi_mem_rep[i] = normal_rng(phi_mem_mu[1], phi_mem_sigma);
            phi_mem_sim[i] = phi_mem_mu[1];
            rho_mem_rep[i] = phi_mem_rep[i] * rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha_rep[i] - 1) / (2 * pi() * alpha_rep[i] * w_rep[i]^2);
            rho_mem_sim[i] = phi_mem_sim[i] * rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha_sim[i] - 1) / (2 * pi() * alpha_sim[i] * w_sim[i]^2);

        }
        else {
            phi_mem_rep[i] = normal_rng(rho_mem_mu[1] * (2 * pi() * alpha * w_rep[i]^2) / (vol_rep[i] * rho_prot_rep[i]), phi_mem_sigma);
            phi_mem_sim[i] =  rho_mem_mu[1] * 2 * pi() * alpha * w_sim[i]^2 / (vol_sim[i] * rho_prot_sim[i]);
            rho_mem_rep[i] = normal_rng(rho_mem_mu[1], rho_mem_sigma[1]);
            rho_mem_sim[i] = rho_mem_mu[1];
        }
        m_mem_rep[i] = phi_mem_rep[i] * prot_per_cell_rep[i];
        m_mem_sim[i] = phi_mem_sim[i] * prot_per_cell_sim[i];
        phi_peri_rep[i] = normal_rng(m_peri_rep[i] / (rho_prot_rep[i] * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1)), phi_peri_sigma);
        phi_peri_sim[i] = m_peri_sim[i] / (rho_prot_sim[i] * (pi()/12) * w_sim[i]^3 * (3 * alpha - 1));
        rho_peri_rep[i] = m_peri_mu / (pi() * alpha_rep[i] * delta * w_rep[i]^2);
        rho_peri_sim[i] = m_peri_mu / (pi() * alpha_sim[i] * delta * w_sim[i]^2);

        // Quantities that are *width* dependent, meaning we have to calculate what the 
        // growth rate is given the volume-scaling relationship
        if (smk) {
            calc_lam_sim[i] = (1/tau[1]) * log(pi() * width_sim[i]^3 * (3 * alpha_sim[i] - 1) / (12 * C_0[1]));
            rel_phi_rho_prot_sim[i] = rho_prot_min + rho_prot_slope * calc_lam_sim[i];

            calc_lam_rep[i] = (1/tau[1]) * log(pi() * width_sim[i]^3 * (3 * alpha_rep[i] - 1) / (12 * C_0[1]));
            rel_phi_rho_prot_rep[i] = rho_prot_min + rho_prot_slope * calc_lam_rep[i];
        }
        else {
            calc_lam_sim[i] = (width_sim[i] - w_min[1]) / w_slope[1];
            calc_lam_rep[i] = calc_lam_sim[i];  
            rel_phi_rho_prot_sim[i] = rho_prot_min + rho_prot_slope * calc_lam_sim[i];
            rel_phi_rho_prot_rep[i] = rho_prot_min + rho_prot_slope * calc_lam_rep[i];
        }

        // Calculation of relative allocation prediction conditioned on width
        if (const_phi_mem) {
                rel_phi_rep[i] = 12 * m_peri_rep[i] / (phi_mem_rep[i] * rel_phi_rho_prot_rep[i] * pi() * (3 * alpha_rep[i] - 1) * width_sim[i]^3);
                rel_phi_sim[i] = 12 * m_peri_sim[i] / (phi_mem_sim[i] * rel_phi_rho_prot_sim[i] * pi() * (3 * alpha_sim[i] - 1) * width_sim[i]^3);
            } 

        else {
                rel_phi_rep[i] = (m_peri_rep[i]) / (pi() * alpha_rep[i] * rho_mem_rep[i] * width_sim[i]^2);
                rel_phi_sim[i] = (m_peri_sim[i]) / (pi() * alpha_sim[i] * rho_mem_sim[i] * width_sim[i]^2); 
                }

        }
} 


