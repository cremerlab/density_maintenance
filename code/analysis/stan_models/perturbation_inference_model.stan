data {
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_biomass;
    int<lower=1> N_cal;
    int<lower=1> N_brad;
    int<lower=1> N_flow;
    int<lower=1> N_growth;

    // Size measurements
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] volumes;
    vector<lower=0>[N_size] peri_volumes;

    // Total protein measurements
    vector<lower=0>[N_prot] total_protein;

}