---
status: Accepted
reason: All looks good, save for the second replicate of glucose+acetate
---

# 2024-05-22 MS Sample Preparation for Wildtype in Varied Growth Media 
This experiment was run to gather samples for total proteome quantification of 
wildtype *E. coli* grown in fast growth media. As part of the sample preparation,
growth rates, total RNA, total Protein, and cell size measurements were acquired.


## Materials & Equipment
### Media

| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| glucose+acetate | N-C- | 10mM glucose + 30mM acetate|
| glucose | N-C-| 10mM glucose |
| sorbitol | N-C-| 10mM sorbitol |
| glycerol | N-C- | 20mM glycerol |
| acetate | N-C- | 30mM acetate | 


### Strains
|**Label**|**Strain Identifier**|**Description**|
|:--:|:--:|:--:|
|WT | GE046 | Wildtype NCM3722 *E. coli*|


## Results
* The second replicate of growth in glucose+acetate was dropped due to timing. I
had to leave before it would be ready for harvest, so I will repeat this run 
tomorrow. 

### Harvest & Residual Optical Densities
|**Label**| **Use**| **Harvest OD 25mm tube (cuvette)**| **Residual Cuvette OD**| **Net OD Units**|
|:--:|:--:|:--:|:--:|:--:|
|glucose 1A | Total Protein | 0.465 (0.223) | 0.039 | 0.276 |
|glucose 1B | Total RNA | 0.465 (0.223) | 0.039 | 0.276 |
|glucose 1M | Mass Spec| 0.465 (0.223) | 0.028 | 2.35 |
|glucose 2A | Total Protein | 0.409 (0.200) | 0.053 | 0.225 |
|glucose 2B | Total RNA | 0.409 (0.200) | 0.053 | 0.225 |
|glucose 2M | Mass Spec | 0.409 (0.200) | 0.026 | 2.08 |
|glucose+acetate 1A | Total Protein | 0.416 (0.204) | 0.040 | 0.248 |
|glucose+acetate 1B | Total RNA | 0.416 (0.204) | 0.069 | 0.202 | 
|glucose+acetate 1M | Mass Spec | 0.416 (0.204) | 0.034 | 2.04| 
|glucose+acetate 2A | Total Protein | 0.352 (0.172) | 0.020 | 0.228 |
|glucose+acetate 2B | Total RNA | 0.352 (0.172) | 0.025 | 0.228 |
|glucose+acetate 2M | Mass Spec | 0.352 (0.172) | 0.020 | 1.8 |
|glycerol 1A | Total Protein | 0.416 (0.204) | 0.009 | 0.292 |
|glycerol 1B | Total RNA | 0.416 (0.204) | 0.009 | 0.292 |
|glycerol 1M | Mass Spec | 0.416 (0.204) | 0.116 | 1.05 |
|glycerol 2A | Total Protein | 0.482 (0.237) | 0.009 | 0.342 |
|glycerol 2B | Total RNA| 0.482 (0.237) | 0.010 | 0.340 |
|glycerol 2M | Mass Spec | 0.482 (0.237) | 0.036 | 2.41 |
|sorbitol 1A | Total Protein | 0.391 (0.192) | 0.015 | 0.266 |
|sorbitol 1B | Total RNA | 0.391 (0.192) | 0.014 | 0.267 | 
|sorbitol 1M | Mass Spec | 0.391 (0.192) | 0.118 |0.888 | 
|sorbitol 2A | Total Protein | 0.413 (0.202) | 0.013 | 0.284 |
|sorbitol 2B | Total RNA | 0.413 (0.202) | 0.021 | 0.272 |
|sorbitol 2M | Mass Spec | 0.413 (0.202) | 0.125 | 0.924 |

### Growth Rate Estimation
|**Label** | **Growth Rate [inv. hr]** |
|:--:|:--:|
|glucose 1 | 0.86 |
|glucose 2 | 0.88 |
|glucose+acetate 1 | 0.85 |
|glucose+acetate 2 | 0.70 |
|glycerol 1| 0.67 |
|glycerol 2 | 0.63 |
|sorbitol 1 | 0.62 |
|sorbitol 2 | 0.61 |

![](viz/2024-05-22_r1_growth_curves.png)

### Cell Segmentation
|**Label**| **Average Width [µm]** | **Average Length [µm]** | **Average SA/V [µm^-1]** | **Average Volume [fL]**|
|:--:|:--:|:--:|:--:|:--:|
| glucose 1 | 0.66 | 2.63 | 6.67 | 0.835|
| glucose 2 | 0.70 | 2.37 | 6.42 | 0.818|
| glucose+acetate 1 | 0.66 | 2.53 | 6.70 | 0.800 |
| glucose+acetate 2 | 0.66 | 2.43 | 6.72 | 0.760 |
| glycerol 1 | 0.59 | 2.42 |  7.43 | 0.615 | 
| glycerol 2 | 0.60 | 2.27 | 7.36 | 0.588 | 
| sorbitol 1 | 0.605 | 2.32 | 7.35 | 0.617 |
| sorbitol 2 | 0.620 | 2.31 | 7.18 | 0.642 |

![](./viz/2024-05-22_r1_size_cdfs.png)

# Protocol
## Cell Husbandry & Growth Measurements
1. Precultures were inoculated by adding a single colony from a recently struck 
wildtype *E. coli* agar plate to 3 mL of growth medium in 16 mm glass test tubes.
For glucose and glucose+acetate sample, 10 mL of growth medium was used in a 25 mm 
glass test tube. These samples were incubated at 37° C in a waterbath with rapid (250 rpm) agitation.
2. Once precultures reached a cuvette-corrected OD$_{600nm}$ of ≈ 0.3, samples 
were diluted into 15 mL of prewarmed media in 25 mm diameter glass tubes to a final
OD of ≈ 0.02.
3. The cells were grown until a cuvette-corrected OD$_{600nm}$ of ≈ 0.4, with 
measurements of the optical density being taken by hand every 15 - 20 minutes. These 
measurements were later used to estimate the growth rate. 

## Sample Preparation For Mass Spectrometry
1. Once samples reached an OD$_{600nm}$ ≈ 0.4, 12 mL of the culture was transferred 
to 14 mL conical falcon tubes. 
2. The cells were pelleted at ≈3000 RPM for 10 minutes. 11 mL of the supernatant 
was removed and transferred to a clean 14 mL falcon tube. 
3. The pellet was gently resuspended in the remaining 1mL of supernatant and transferred
to a 1.7mL eppendorf tube.
4. The resuspended cell pellet was spun at 15000 RPM for 2 min and the supernatant 
was removed and mixed with the 11 mL of residual supernatant.  The optical density 
of this combined supernatant was measured using a qcm quartz cuvette.
5. The total number of OD units was calculated and written on the falcon tube, 
and was then transferred to the -80 °C freezer. 

## Single Cell Imaging 
1. During exponential growth, 1.5% agarose pads in N+C- buffer base were prepared. 
Briefly, 1 mL of the mixture was made molten in a 98° C heat block. 50µL of the
agar slurry was deposited onto a glass slide with a coverslip on each side, and 
then sandwiched with another glass coverslip. The pad was allowed to dry for 10 - 
20 minutes.
2. When the cultures reached OD ≈ 0.15 - 0.2, 2 µL of the culture was added 
to  the agar slab on each slide. The slides were then transferred to the 
microscope and imaged with a lamp power of 9.8V (3200K).
3. After approximately 10^2 cells were imaged (10 - 15 images), the images were 
transferred to cloud storage.

## Sample Preparation for Total RNA and Total Protein Quantification 
1. While samples were spinning for the mass spec preparation, 2 x 1.5 mL aliquots 
from the culture were transferred into clean 1.7 mL eppendorf tubes.  
2. The samples were spun at 15000 RPM for 2 min. 
3. The supernatant was removed and the residual optical density was measured 
in the quartz cuvette. 
4. The samples were then transferred to the -80° C freezer until total RNA 
and protein could be quantified. 

