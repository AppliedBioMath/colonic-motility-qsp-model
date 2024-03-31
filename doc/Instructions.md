---
header-includes:
  - \usepackage[margin=1in]{geometry}
  - \usepackage[T1]{fontenc}
  - \usepackage[charter, 11pt]{mathdesign}
  - \usepackage{amsmath, graphicx}
  - \everymath{\displaystyle}
  - \makeatletter
  - \def\hrulefill{\leavevmode\leaders \hrule height \rulethickness \hfill\kern\z@}
  - \makeatletter
---

# Colinic motility QSP model MATLAB code

### 2023-05-12

\hrulefill

## Welcome!

This folder contains MATLAB code to simulate the colonic motility 
quantitative systems pharmacology  model developed by 
Applied BioMath, LLC and Takeda Pharmaceuticals. The model is descibed 
in this publication: https://pubmed.ncbi.nlm.nih.gov/31432345/

## License

Colonic motility QSP model code: 2023-05-12

MIT License

Copyright (c) 2023 Raibatak Das and Applied BioMath, LLC

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

## Installation

Ensure that you have a recent version of MATLAB. This code has been tested with Release 2018b on PC and Mac.

1. Copy/download this entire folder to a convenient location on your machine.
2. Start MATLAB. Navigate to this folder.
3. Run the script `runMeFirst`.  $\longleftarrow$ **Always run this at the start of a session**. This will update your MATLAB path and print a short description on the console

## Contents

There are several demo scripts included:

- **`runHealthy`**  
        Simulate colon physiology for healthy phenotype. This will produce Fig. 2A of GI motility manuscript.

- **`runSTC`**  
        Simulate colon physiology for STC phenotype. This will produce Fig. 2B of GI motility manuscript.

- **`doseResponseHealthy`**  
        Simulate dose response for healthy phenotype. This will produce Fig. 4C of GI motility manuscript.

- **`doseResponseSTC`**  
        Simulate dose response for STC phenotype. This will produce Fig. 4D of GI motility manuscript.
        
- **`GImodeldemo`**  
        Run all of the above demos. 

These scripts use parameters defined in the following files:

- **`parameter-table-healthy.csv`**  
        Parameter table for healthy phenotype

- **`parameter-table-stc.csv`**  
        Parameter table for STRC phenotype

- **`hapc-schedule.csv`**  
        Timings for HAPCs and defecation checks

The recommended way to change simulation parameters is to create a new copy of a parameter file and update the simulation scripts to use the new parameter file. The file `parameter-table-test.csv` is included for this purpose. To change parameters programmatically, see the examples in `runParameterChange.m`

