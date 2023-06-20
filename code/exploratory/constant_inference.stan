 data { 
         int<lower=1> N_ms;
         int<lower=1> N_size;
         vector<lower=0>[N_ms] rho_mem;
         vector<lower=0>[N_ms] phi_mem;
         vector<lower=1>[N_size] alpha;
    }

    parameters {
        real<lower=0> rho_mem_mu;
        real<lower=0> rho_mem_sigma;
        real<lower=0> phi_mem_mu;
        real<lower=0> phi_mem_sigma;
        real<lower=1> alpha_mu;
        real<lower=1> alpha_sigma;
    }

    model {
    rho_mem_mu ~ std_normal();
    rho_mem_sigma ~ std_normal();
    phi_mem_mu ~ std_normal();
    phi_mem_sigma ~ std_normal();
    alpha_mu ~ std_normal();
    alpha_sigma ~ std_normal();

    rho_mem ~ normal(rho_mem_mu, rho_mem_sigma);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    alpha ~ normal(alpha_mu, alpha_sigma);
}