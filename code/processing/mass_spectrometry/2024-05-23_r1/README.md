---
status: Accepted
reason: All looks good except for the second replicate of pRelA + 2 ng/mL dox.
---

# 2024-05-23 MS Sample Preparation for ppGpp Perturbations
This experiment was run to gather samples for total proteome quantification of 
*E. coli* grown in glucose media under various ppGpp perturbations. As part of the sample preparation,
growth rates, total RNA, total protein, and cell size measurements were acquired.


## Materials & Equipment
### Media

| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| glucose | N-C- | 10mM glucose |

### Strains
|**Label**|**Strain Identifier**|**Description**|
|:--:|:--:|:--:|
|meshI | GE463 | Wildtype NCM3722 *E. coli* with pMeshI|
|relA | GE462 | Wildtype NCM3722 *E. coli* with pRelA|

## Results
* The second replicate of pRelA + 2 ng/mL dox had a surprisingly slow growth 
rate, comparable to an induction level of 4 ng/mL as can be seen in the experiments 
from 2024-05-25. This replicate is dropped from this analysis.

### Harvest & Residual Optical Densities
|**Label**| **Use**| **Harvest OD 25mm tube (cuvette)**| **Residual Cuvette OD**| **Net OD Units**|
|:--:|:--:|:--:|:--:|:--:|
|meshI 0 µM IPTG 1A | Total Protein | 0.440 (0.216) | 0 | 0.324 |
|meshI 0 µM IPTG 1B | Total RNA | 0.440 (0.216) | 0 | 0.324 |
|meshI 0 µM IPTG 1M | Mass Spec | 0.440 (0.216) | 0.013 | 2.44 |
|meshI 0 µM IPTG 2A | Total Protein | 0.470 (0.230) | 0.003| 0.340|
|meshI 0 µM IPTG 2B | Total RNA | 0.470 (0.230) | 0| 0.345|
|meshI 0 µM IPTG 2M | Mass Spec | 0.470 (0.230) | 0.019| 2.53 |
|meshI 100 µM IPTG 1A | Total Protein | 0.467 (0.229) | 0.006 | 0.335 | 
|meshI 100 µM IPTG 1B | Total RNA | 0.467 (0.229) | 0.008 | 0.332 | 
|meshI 100 µM IPTG 1M | Mass Spec| 0.467 (0.229) | 0.004 | 2.70 | 
|meshI 100 µM IPTG 2A | Total Protein | 0.500 (0.245) | 0.006 | 0.359 | 
|meshI 100 µM IPTG 2B | Total RNA | 0.500 (0.245) | 0.006 | 0.359 | 
|meshI 100 µM IPTG 2M | Mass Spec | 0.500 (0.245) | 0.004 | 2.90 | 
|relA  0 ng/mL dox 1A | Total Protein | 0.404 (0.198) | 0.017| 0.272 |
|relA  0 ng/mL dox 1B | Total RNA | 0.404 (0.198) | 0.085| 0.170|
|relA  0 ng/mL dox 1M | Mass Spec | 0.404 (0.198) | 0.016| 2.18|
|relA  0 ng/mL dox 2A | Total Protein | 0.413 (0.202) | 0.034| 0.252|
|relA  0 ng/mL dox 2B | Total RNA | 0.413 (0.202) | 0.037| 0.248|
|relA  0 ng/mL dox 2M | Mass Spec | 0.413 (0.202) | 0.008| 2.33|
|relA  2 ng/mL dox 1A | Total Protein | 0.441 (0.216) | 0.033| 0.275|
|relA  2 ng/mL dox 1B | Total RNA | 0.441 (0.216) | 0.021| 0.293 |
|relA  2 ng/mL dox 1M | Mass Spec | 0.441 (0.216) | 0.014| 2.43 |
|relA  2 ng/mL dox 2A | Total Protein | 0.412 (0.202) | 0.016| 0.279|
|relA  2 ng/mL dox 2B | Total RNA | 0.412 (0.202) | 0.056| 0.219|
|relA  2 ng/mL dox 2M | Mass Spec  | 0.412 (0.202) | 0.016| 2.23 |

### Growth Rate Estimation
|**Label** | **Growth Rate [inv. hr]** |
|:--:|:--:|
| meshI 0µM IPTG 1 | 0.89 |
| meshI 0µM IPTG 2 | 0.82 | 
| meshI 100µM IPTG 1 | 0.72 |
| meshI 100µM IPTG 2 | 0.67 |
| relA 0 ng/mL dox 1 | 0.83 | 
| relA 0 ng/mL dox 2 | 0.77 |
| relA 2 ng/mL dox 1 | 0.60 |
| relA 2 ng/mL dox 2 | 0.48 |

![](viz/2024-05-23_r1_growth_curves.png)

### Cell Segmentation
|**Label**| **Average Width [µm]** | **Average Length [µm]** | **Average SA/V [µm^-1]** | **Average Volume [fL]**|
|:--:|:--:|:--:|:--:|:--:|
| meshI 0µM IPTG 1 | 0.66 | 2.66 | 6.72 | 0.84 |
| meshI 0µM IPTG 2 | 0.71 | 2.39 | 6.30 | 0.86 | 
| meshI 100µM IPTG 1 | 0.75 | 2.39 | 6.03 | 0.95 |
| meshI 100µM IPTG 2 | 0.75 | 2.51 | 6.01 | 1.00 |
| relA 0 ng/mL dox 1 | 0.69 | 2.40 | 6.51 | 0.81 | 
| relA 0 ng/mL dox 2 | 0.65 |  2.47 |  6.89 | 0.75 |
| relA 2 ng/mL dox 1 | 0.66 | 2.17 | 6.78 | 0.67 |
| relA 2 ng/mL dox 2 | 0.61 | 2.15 | 27.32 | 0.58 |

![](./viz/2024-05-23_r1_size_cdfs.png)

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

## Quantification of Total RNA


## Quantification of Total Protein