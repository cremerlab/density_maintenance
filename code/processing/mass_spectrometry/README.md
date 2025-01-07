#`processing/mass_spectrometry`

This directory contains all experimental data sets (images, growth rates, RNA/Protein,
and mass spectrometry) measurements conducted for the main analysis of this work. 

All folders beginning with `2024-MM-DD_r#` contain that day's experimental 
runs, with `r#` indicating the replicate conducted within that single day. 
The days  `2022-07-21_r1` and `2022-10-25_r1` contain the processing Python 
scripts for analyzing and standardizing the results from mass spectrometry of 
samples collected in the previous day.  

Not all folders contained here are used in the final analysis. The file 
`valid_experimental_metadata.csv` contains information of all of the data sets 
that were free of experimental errors or obvious outliers in measurements (e.g. 
clearly non-exponential growth curves or significantly slow/fast growth rates 
inconsistent with previous results).

Note that total RNA and total Protein measurements were collected for most mass
spectrometry samples as well. However, the sample storage protocol omitted an 
important washing step (before freezing cell pellets) which non-systematically 
skews the resulting quantitation. As it is not easy to systematically identify 
the effect, total RNA and total Protein from these samples were not used in the 
downstream analysis. 


