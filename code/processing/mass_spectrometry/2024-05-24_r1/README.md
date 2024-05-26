---
status: rejected
reason: Experiment in progress
---

# 2024-05-24 MS Sample Preparation for Fast Growth 
This experiment was run to gather samples for total proteome quantification of 
*E. coli* grown in fast growth media. As part of the sample preparation,
growth rates, total RNA, total protein, and cell size measurements were acquired.


## Materials & Equipment
### Media

| **Label** | **Buffer Base** | **Carbon Source & Concentration** |
|:--:|:--:|:--:|
| glucoseCAA | N-C- | 10mM glucose + 0.1% (w/v) CAA|
| LB | -- | Rich Undefined | 

### Strains
|**Label**|**Strain Identifier**|**Description**|
|:--:|:--:|:--:|
|wildtype | GE046 | Wildtype NCM3722 *E. coli*|

## Results

### Harvest & Residual Optical Densities
|**Label**| **Use**| **Harvest OD 25mm tube (cuvette)**| **Residual Cuvette OD**| **Net OD Units**|
|:--:|:--:|:--:|:--:|:--:|
|LB 1A | Total Protein |0.445 (0.218) | 0.018 | 0.300 |
|LB 1B | Total RNA |0.445 (0.218) | 0.021 | 0.296 |
|LB 1M | Mass Spec |0.445 (0.218) | 0.020 | 2.38 |
|LB 2A | Total Protein | 0.500 (0.245) | 0.014 | 0.347 |
|LB 2B | Total RNA | 0.500 (0.245) | 0.046 | 0.299 |
|LB 2M | Mass Spec | 0.500 (0.245) | 0.030 | 2.58 | 
|glucoseCAA 1A | Total Protein |0.412 (0.202) | 0.026 | 0.264 |
|glucoseCAA 1B | Total RNA |0.412 (0.202) | 0.043 | 0.239|
|glucoseCAA 1M | Mass Spec |0.412 (0.202) | 0.030 | 2.06 |
|glucoseCAA 2A | Total Protein | 0.453 (0.222) | 0.031 | 0.287 |
|glucoseCAA 2B | Total RNA | 0.453 (0.222) | 0.036 | 0.279 |
|glucoseCAA 2M | Mass Spec| 0.453 (0.222) | 0.023 | 2.39|

### Growth Rate Estimation
|**Label** | **Growth Rate [inv. hr]** |
|:--:|:--:|
| LB 1 | 2.27 | 
| LB 2 | 2.27 | 
| glucoseCAA 1 | 1.48 |
| glucoseCAA | 1.30|

![](viz/2024-05-24_r1_growth_curves.png)

### Cell Segmentation
|**Label**| **Average Width [µm]** | **Average Length [µm]** | **Average SA/V [µm^-1]** | **Average Volume [fL]**|
|:--:|:--:|:--:|:--:|:--:|
| LB 1 | 0.92 | 3.85  | 4.77 | 2.37 |
| LB 2 | 0.92 | 4.21 | 4.73 | 2.69 |
| glucoseCAA 1 | 0.87 | 3.18 | 5.17 | 1.77 | 
| glucoseCAA 2 | 0.86  | 3.34 | 5.17 | 1.83 |


![](./viz/2024-05-24_r1_size_cdfs.png)

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