# Coordination Between Cytoplasmic and Envelope Densities Shapes Cellular Geometry in *Escherichia coli*

[![DOI](https://zenodo.org/badge/438837626.svg)](https://zenodo.org/doi/10.5281/zenodo.10048570)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Welcome to the research repository for our recent preprint "Coordination Between Cytoplasmic and Envelope Densities Shapes Cellular Geometry in *Escherichia coli*!"
This repo houses all code and (almost) all data used for the findings reported in 
the work. Microscopy images used in the final analysis are large (≈ 6GB) and hosted separately on the
Stanford Research Data Repository accessible via doi: 10.25740/ws785mz0287. All experimental data files, including those that were not included in the final analysis, are also hosted on the Stanford Research Data Repository via doi: 10.25740/ty214py1791.

This repository is still in an active research phase and will be further developed 
over the course of the review process. Below is a brief summary of each 
subdirectory and what it houses. 

## `code`
It will not shock you to hear that this houses all code used in this work. This 
code is almost exclusively written in Python, though statistical models were 
defined in the Stan probabilistic programming language. Within this directory 
are several subdirectories as follows:

* `analysis`: This directory houses all code used in the statistical analyses 
presented in this work. It is further broken down into segments for the inference 
of the model parameters and inference of the data parameters, primarily from 
our own experiments. 

* `processing`: This directory houses scripts used in the processing and cleaning 
of literature data used in this work and in the processing of raw experimental 
measurements. 

* `visualization`: This directory houses Python scripts used in generating all 
of the figures in this work. 

## `data`
This is where all of the (small) datasets used in this work live. It is broken 
down in to groups of literature data and our own measurements.

## `figures`
This directory houses all of the figures used in the manuscript, usually separated 
by figure panel. As of now, this houses many different versions used in different 
stages of the manuscript.

## `protocols`
This directory houses markdown files with media recipes and detailed step-by-step
guides for the various experimental protocols performed in this work.

## `software`
This project required the development of an independent software package, primarily 
for the image processing. This module (named `size`) is packaged and documented
as a typical Python software package.

Installing this package is necessary to run effectively all of the code used in 
this work. To install locally via pip, you can run the following from the command 
line in the root directory:

```
pip install -e ./software/
```

## License
All creative work (e.g. prose, artwork) is licensed under a [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.
The software herein is licensed under a standard MIT license, which reads as
follows:

```
MIT License

Copyright (c) 2025 The Authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

