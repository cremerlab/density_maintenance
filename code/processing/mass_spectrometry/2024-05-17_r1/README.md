---
status: rejected
reason: Experiment in progress
---

# 2024-05-17 MS Sample Preparation for Wildtype in Glucose + CAA
This experiment was run to gather samples for total proteome quantification of 
wildtype *E. coli* grown in fast growth media. As part of the sample preparation,
growth rates, total RNA, total Protein, and cell size measurements were acquired.

## Materials & Equipment
### Media

| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| glucoseCAA| N-C- | 10mM glucose + 0.1% (w/v) casamino acids|


### Strains
|**Label**|**Strain Identifier**|**Description**|
|:--:|:--:|:--:|
|WT | GE046 | Wildtype NCM3722 *E. coli*|


## Results

### Residual Optical Densities
|**Label**| **Use**| **Residual OD**|
|:--:|:--:|:--:|
|glucose+CAA 1A | Total Protein | 0.051 |
|glucose+CAA 1B | Total RNA | 0.035 |
|gluose+CAA 1M| Mass Spec.| 0.037 |
|glucose+CAA 2A | Total Protein | 0.016 |
|glucose+CAA 2B | Total RNA | 0.017 |
|glucose+CAA 2M | Mass Spec.| 0.058 |

### Total OD For Mass Spec 
|**Label** | **Harvest OD in 25mm Tube (in cuvette)** | **Residual OD in cuvette (in 25mm tube)** | **Net OD Units** |
|:--:|:--:|:--:|:--:|
|glucose+CAA 1M | 0.405 (0.49 • 0.405 = 0.198) | 0.037 (0.037/0.49 = 0.075) | 12 • (0.198 - 0.037) = 1.932 |
|glucose+CAA 2M | 0.4 (0.49 • 0.4 = 0.196) | 0.058 (0.058/0.49 = 0.118) | 12 • (0.196 - 0.058) = 1.656 | 

### Growth Rate Estimation
|**Label** | **Growth Rate [inv. hr]** |
|:--:|:--:|
| glucose+CAA rep 1 | 1.13 |
| glucose+CAA rep 2 | 0.96 |

![](viz/2024-05-17_r1_growth_curves.png)

### Cell Segmentation
|**Label**| **Average Width [µm]** | **Average Length [µm]** | **Average SA/V [µm^-1]** | **Average Volume [fL]**|
|:--:|:--:|:--:|:--:|:--:|
|glucose+CAA 1 | 0.824 | 3.238 | 5.367 | 1.603 |
|glucose+CAA 2 | 0.872 | 3.074 | 5.113 | 1.666 |

![](./viz/2024-05-17_r1_size_cdfs.png)

# Protocol
## Cell Husbandry & Growth Measurements
1. Precultures were inoculated by adding a single colony from a recently struck 
wildtype *E. coli* agar plate to 3 mL of growth medium in 16 mm glass test tubes.
These samples were incubated at 37° C in a waterbath with rapid (250 rpm) agitation.
2. Once precultures reached a cuvette-corrected OD$_{600nm}$ of ≈ 0.3, samples 
were diluted into 15 mL of prewarmed media in 25 mm diameter glass tubes to a final
OD of ≈ 0.02.
3. The cells were grown until a cuvette-corrected OD$_{600nm}$ of ≈ 0.4, with 
measurements of the optical density being taken by hand every 15 - 20 minutes. These 
measurements were later used to estimate the growth rate. 

## Sample Preparation For Mass Spectrometry
1. Once samples reached an OD$_{600nm}$ ≈ 0.4, 12 mL of the culture was tranferred 
to 14 mL conical falcon tubes. 
2. The cells were pelleted at ≈3000 RPM for 10 minutes. 
3. The supernatant was removed using a serological pipette and the optical 
density was recorded using a quartz cuvette. 
4. The total number of OD units was calculated and written on the falcon tube, 
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