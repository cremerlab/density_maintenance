data {
    // Dimensional information
    int<lower=1> J_size; // Number of unique size conditions
    int<lower=1> J_growth; // Number of unique growth conditons
    int<lower=1> J_brad; // Number of unique bradford condition
    int<lower=1> J_flow; // Number of unique flow cytometry conditions
    int

    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_brad;// Number of bradford measurements
    int<lower=1> N_growth; // Number of growth measurements
    int<lower=1> N_flow; // Number of unique flow cytometry measurements

    // ID vectors
    array[N_size] int<lower=1, upper=J_size> size_idx; 
    array[N_growth] int<lower=1, upper=J_growth> growth_idx;
    array[N_flow] int<lower=1, upper=J_flow> flow_idx;
    array[N_brad] int<lower=1, upper=J_brad> brad_idx;
}
