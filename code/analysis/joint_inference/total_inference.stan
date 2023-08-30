data {
    // Test data
    int<lower=1> N_pred;
    vector<lower=0>[N_pred] pred_phiRb;

    //--------------------------------------------------------------------------
    // LITERATURE DATA
    //--------------------------------------------------------------------------

    // Literature size data.
    int<lower=1> N_size_lit;
    vector<lower=0>[N_size_lit] lit_size_lam;
    vector<lower=0>[N_size_lit] lit_SAV;

    // Protein-per-cell data 
    int<lower=0> N_prot_lit;
    vector<lower=0>[N_prot_lit] lit_prot_lam;
    vector<lower=0>[N_prot_lit] lit_prot_per_cell;

    // Drymass density data
    int<lower=0> N_drymass_lit;
    vector<lower=0>[N_drymass_lit] lit_drymass_lam;
    vector<lower=0>[N_drymass_lit] lit_drymass_density;

    // Mass spec data
    int<lower=0> N_ms;
    vector<lower=0>[N_ms] lit_ms_lam;
    vector<lower=0>[N_ms] lit_phi_mem;
    vector<lower=0>[N_ms] lit_phi_peri;


    // PhiRb data
    int<lower=1> N_phiRb_lit;
    vector<lower=0>[N_phiRb_lit] lit_phiRb;
    vector<lower=0>[N_phiRb_lit] lit_phiRb_lam;

    -------------------------------------------------------------------------- 
    //  OUR MEASUREMENTS
    //-------------------------------------------------------------------------- 

    // Total protein measurements
    int<lower=1>N_prot;
    int<lower=1>J_prot;
    array[N_prot] int<lower=1, upper=J_prot> prot_idx;
    vector<lower=0>[N_prot] prot_per_biomass;

    // Periplasmic protein measurements
    int<lower=1>N_peri;
    int<lower=1>J_peri;
    array[N_peri] int<lower=1, upper=J_peri> peri_idx;
    vector<lower=0>[N_peri] peri_per_biomass;

    // Membrane protein measurements
    int<lower=1>N_mem;
    int<lower=1>J_mem;
    array[N_mem] int<lower=1, upper=J_mem> mem_idx;
    vector<lower=0>[N_mem] mem_per_biomass;

    // Total RNA measurements
    int<lower=1>N_rna;
    int<lower=1>J_rna;
    array[N_rna] int<lower=1, upper=J_rna> rna_idx;
    vector<lower=0>[N_rna] rna_per_biomass;


    // Growth rate measurement
    int<lower=1> N_growth;
    int<lower=1> J_growth;
    array[N_growth] int<lower=1, upper=J_growth> growth_idx;
    vector<lower=0>[N_growth] growth_rates;

    // Size measurements
    int<lower=1> N_size;
    int<lower=1> J_size;
    array[N_size] int<lower=1, upper=J_size> size_idx;
    // vector<lower=0>[N_size]  SAV;
    // vector<lower=0>[N_size] widths;
    // vector<lower=0>[N_size] lengths;
    // vector<lower=0>[N_size] aspect_ratios;
    vector<lower=0>[N_size] surface_areas;
    vector<lower=0>[N_size] volumes;

    // Flow cytometry measurements
    int<lower=1> N_flow;
    int<lower=1> J_flow;
    array[N_flow] int<lower=1, upper=J_flow> flow_idx;
    vector<lower=0>[N_flow] cells_per_biomass;
}


parameters {
    //-------------------------------------------------------------------------- 
    //  LITERATURE PARAMETERS
    //-------------------------------------------------------------------------- 
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> phi_peri_sigma;
    real<lower=0> kappa;
    real<lower=0> lit_drymass_mu;
    real<lower=0> lit_drymass_sigma;
    real<lower=0> rho_mem_sigma;
    real<lower=0> SAV_intercept;
    real<upper=0> SAV_slope;
    real<lower=0> SAV_sigma;
    real<lower=0> phiRb_intercept;
    real<lower=0> phiRb_slope;
    real<lower=0> phiRb_sigma;
    real log_prot_intercept;
    real<lower=0> log_prot_slope;
    real<lower=0> log_prot_sigma;

    //-------------------------------------------------------------------------- 
    //  OUR PARAMETERS
    //-------------------------------------------------------------------------- 

    // Total protein parameters
    vector<lower=0>[J_prot] prot_per_biomass_mu;
    vector<lower=0>[J_prot] prot_per_biomass_sigma;

    // Membrane protein parameters
    vector<lower=0>[J_mem] mem_per_biomass_mu;
    vector<lower=0>[J_mem] mem_per_biomass_sigma;

    // Periplasmic protein parameters
    vector<lower=0>[J_peri] peri_per_biomass_sigma;

    // RNA parameters
    vector<lower=0>[J_rna] rna_per_biomass_mu;
    vector<lower=0>[J_rna] rna_per_biomass_sigma;

    // Growth parameters
    vector<lower=0>[J_growth] growth_rate_mu;
    vector<lower=0>[J_growth] growth_rates_sigma;

    // Flow cytometry parameters
    vector<lower=0>[J_flow] cells_per_biomass_mu;
    vector<lower=0>[J_flow] cells_per_biomass_sigma;

    // Size parameters
    vector<lower=0>[J_size] surface_area_mu;
    vector<lower=0>[J_size] surface_area_sigma;
    vector<lower=0>[J_size] volume_mu;
    vector<lower=0>[J_size] volume_sigma;

}

transformed parameters {
    vector<lower=0>[N_ms] rho_mem_ms = lit_phi_mem .* exp(log_prot_intercept + log_prot_slope .* lit_ms_lam) ./ (2 .* (SAV_intercept + SAV_slope .* lit_ms_lam)); 
    vector<lower=0>[J_size] SAV = surface_area_mu ./ volume_mu;
    vector<lower=0>[J_mem] phi_mem = mem_per_biomass_mu ./ prot_per_biomass_mu;
    vector<lower=0>[J_mem] rho_mem = 1E9 .* mem_per_biomass_mu ./ (2 .* surface_area_mu .* 1E9 .* cells_per_biomass_mu);

}

model {

    //-------------------------------------------------------------------------- 
    //  LITERATURE DATA MODEL
    //-------------------------------------------------------------------------- 
    phi_mem_mu ~ beta(2.5, 8);
    phi_mem_sigma ~ normal(0, 0.1);
    lit_phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma); 

    m_peri_mu ~ normal(10, 2);
    phi_peri_sigma ~ std_normal();
    lit_phi_peri ~ normal(m_peri_mu ./ exp(log_prot_intercept + log_prot_slope .* lit_ms_lam), phi_peri_sigma);    

    kappa ~ normal(100, 20);
    lit_drymass_mu ~ normal(300, 20);
    lit_drymass_sigma ~ normal(0, 10);
    lit_drymass_density ~ normal(lit_drymass_mu, lit_drymass_sigma);
    rho_mem_sigma ~ normal(0, 1);
    rho_mem_ms ~ normal(lit_drymass_mu ./ kappa, rho_mem_sigma);
    // lit_drymass_density ~ normal(rho_mem_mu .* kappa, lit_drymass_sigma);

    SAV_intercept ~ normal(0, 10);
    SAV_slope ~ std_normal();
    SAV_sigma ~ std_normal();
    lit_SAV ~ normal(SAV_intercept + SAV_slope .* lit_size_lam, SAV_sigma);

    phiRb_slope ~ std_normal();
    phiRb_intercept ~ std_normal();
    phiRb_sigma ~ normal(0, 0.1);
    lit_phiRb ~ normal(phiRb_intercept + phiRb_slope .* lit_phiRb_lam, phiRb_sigma);

    log_prot_intercept ~ std_normal();
    log_prot_slope ~ std_normal();
    log_prot_sigma ~ std_normal();
    log(lit_prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* lit_prot_lam, log_prot_sigma);

    //-------------------------------------------------------------------------- 
    //  OUR DATA MODEL
    //-------------------------------------------------------------------------- 

    // Total protein inference model
    prot_per_biomass_mu ~ normal(0, 300);
    prot_per_biomass_sigma ~ std_normal();
    prot_per_biomass ~ normal(prot_per_biomass_mu[prot_idx], prot_per_biomass_sigma[prot_idx]);

    // Membrane protein inference model
    mem_per_biomass_mu ~ normal(0, 20);
    mem_per_biomass_sigma ~ std_normal();
    mem_per_biomass ~ normal(mem_per_biomass_mu[mem_idx], mem_per_biomass_sigma[mem_idx]);

    // Periplasmic protein inference model
    peri_per_biomass_sigma ~ normal(0, 10);
    peri_per_biomass ~ normal(m_peri_mu .* cells_per_biomass_mu[peri_idx], peri_per_biomass_sigma[peri_idx]);

    // RNA inference model
    rna_per_biomass_mu ~ normal(0, 200);
    rna_per_biomass_sigma ~ std_normal();
    rna_per_biomass ~ normal(rna_per_biomass_mu[rna_idx], rna_per_biomass_sigma[rna_idx]);

    // Growth rate inference model
    growth_rate_mu ~ std_normal();
    growth_rates_sigma ~ std_normal();
    growth_rates ~ normal(growth_rate_mu[growth_idx], growth_rates_sigma[growth_idx]);

    // FLow cytometry inference model
    cells_per_biomass_mu ~ std_normal();
    cells_per_biomass_sigma ~ normal(0, 0.1);
    cells_per_biomass / 1E9 ~ normal(cells_per_biomass_mu[flow_idx], cells_per_biomass_sigma[flow_idx]);

    // Size inference model
    volume_mu ~ normal(0, 3);
    volume_sigma ~ normal(0, 0.5);
    volumes ~ normal(volume_mu[size_idx], volume_sigma[size_idx]);

    surface_area_mu ~ normal(0, 3);
    surface_area_sigma ~ std_normal();
    surface_areas ~ normal(surface_area_mu[size_idx], surface_area_sigma[size_idx]);
}

generated quantities {
    // Theory prediction
    vector<lower=0>[N_pred] pred_lam = (pred_phiRb - phiRb_intercept) ./ phiRb_slope;
    vector<lower=0>[N_pred] pred_prot = exp(log_prot_intercept + log_prot_slope .* pred_lam);
    vector<lower=0>[N_pred] pred_SAV = phi_mem_mu .* kappa ./ (2 .* ((pred_phiRb./0.4558) + 1 - phi_mem_mu - m_peri_mu./pred_prot));
    
    // Estimated SAV from phiRb
    vector<lower=0>[N_phiRb_lit] phiRb_SAV;
    for (i in 1:N_phiRb_lit) {
        phiRb_SAV[i] = SAV_intercept + SAV_slope * lit_phiRb_lam[i];
    }
}