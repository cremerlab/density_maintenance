---
status: rejected
reason: Experiment in progress
---

# 2024-05-29 MS Sample Preparation for LacZ OE
This experiment was run to gather samples for total proteome quantification of 
*E. coli* grown in glucose medium overexpressing lacZ at a mild induction level.

## Materials & Equipment
### Media

| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| glucose | N-C- | 10mM glucose|


### Strains
|**Label**|**Strain Identifier**|**Description**|
|:--:|:--:|:--:|
|lacZ | GE246 | Wildtype NCM3722 *E. coli* with pZA31-tetR + pZE11-lacZ|

## Results
### Harvest & Residual Optical Densities
|**Label**| **Use**| **Harvest OD 25mm tube (cuvette)**| **Residual Cuvette OD**| **Net OD Units**|
|:--:|:--:|:--:|:--:|:--:|
| lacZ + 3 ng/mL ctc 1A | Total Protein | 0.415 (0.203) | 0.005 | 0.297 | 
| lacZ + 3 ng/mL ctc 1B | Total RNA | 0.415 (0.203) | 0.007 | 0.294 | 
| lacZ + 3 ng/mL ctc 1M | Mass Spec | 0.415 (0.203) | 0.009 | 2.33 | 



### Growth Rate Estimation
|**Label** | **Growth Rate [inv. hr]** |
|:--:|:--:|
| lacZ + 3 ng/mL ctc, rep 1 | 0.46 |



![](viz/2024-05-29_r1_growth_curves.png)

### Cell Segmentation
|**Label**| **Average Width [µm]** | **Average Length [µm]** | **Average SA/V [µm^-1]** | **Average Volume [fL]**|
|:--:|:--:|:--:|:--:|:--:|
| lacZ + 3 ng/mL ctc, rep 1| 0.78 | 2.9 | 5.75 | 1.29|

![](./viz/2024-05-29_r1_size_cdfs.png)

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