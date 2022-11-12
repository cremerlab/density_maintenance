data {
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] m_cyto; // Cytoplasmic protein
    vector<lower=0>[N] m_peri; // Periplasmic protein protein
}

parameters {
    // Hyperparameters
    real<lower=0> m_cyto_mu;
    real<lower=0> m_peri_mu;
    real<lower=0> homosced_sigma;

}

model {
    m_cyto_mu ~ normal(0, 100);
    m_peri_mu ~ normal(0, 10);
    homosced_sigma ~ normal(0, 1);
    m_cyto ~ normal(m_cyto_mu, homosced_sigma);
    m_peri ~ normal(m_peri_mu, homosced_sigma);
}