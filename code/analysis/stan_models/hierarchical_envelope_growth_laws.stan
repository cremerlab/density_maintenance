data { 
    int<lower=1> N_size;
    int<lower=1> J_size;
    array[N_size] int<lower=1, upper=J_size> size_idx;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_size] aspect_ratio;
}

parameters {
    real tau; 
    real<lower=1> alpha_mu;
    vector[J_size] alpha_mu_1_tilde; 
    vector<lower=0>[J_size] alpha_sigma;
}

transformed parameters {
    vector<lower=1>[J_size] alpha_mu_1 = alpha_mu + tau .* alpha_mu_1_tilde;
}

model {
    tau ~ std_normal();

    // Aspect Ratio Measurement 
    alpha_mu ~ normal(4, 1);
    alpha_mu_1_tilde ~ normal(0, 0.5);
    alpha_sigma ~ std_normal();
    aspect_ratio ~ normal(alpha_mu_1[size_idx], alpha_sigma[size_idx]);
}


// data {
//     // Protein per cell information
//     int<lower=1> J_prot_per_cell;
//     int<lower=1> N_prot_per_cell;
//     array[N_prot_per_cell] int<lower=1, upper=J_prot_per_cell> prot_per_cell_idx;
//     vector<lower=0>[N_prot_per_cell] prot_per_cell_growth_rate;
//     vector<lower=0>[N_prot_per_cell] prot_per_cell;
// }

// transformed data {
//     vector<lower=0>[N_prot_per_cell] log_prot_per_cell = log(prot_per_cell);
// } 

// parameters {
//     real tau; // Hyper parameter variance

//     // Hierarchical protein per cell inference
//     real<lower=0> log_prot_per_cell_slope;
//     real log_prot_per_cell_intercept; 
//     vector<lower=0>[J_prot_per_cell] log_prot_per_cell_slope_1_tilde;
//     vector<lower=0>[J_prot_per_cell] log_prot_per_cell_intercept_1_tilde;
//     vector<lower=0>[J_prot_per_cell] prot_per_cell_sigma; 

// }

// transformed parameters {
//     vector[J_prot_per_cell] log_prot_per_cell_slope_1 = log_prot_per_cell_slope + tau * log_prot_per_cell_slope_1_tilde;
//     vector[J_prot_per_cell] log_prot_per_cell_intercept_1 = log_prot_per_cell_intercept + tau * log_prot_per_cell_intercept_1_tilde;
// }

// model {
//     tau ~ std_normal();
//     log_prot_per_cell_slope ~ std_normal();
//     log_prot_per_cell_intercept ~ std_normal();
//     log_prot_per_cell_slope_1_tilde ~ std_normal();
//     log_prot_per_cell_slope_1_tilde ~ std_normal();

//     log_prot_per_cell = log_prot_per_cell_intercept_1[prot_per_cell_idx] + log_prot_per_cell_slope_1[prot_per_cell_idx] .* prot_per_cell_growth_rate
// }