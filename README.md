# ARGOS-RAL
 Automatically Extracting Partial Differential Equations from Data

- All testing codes are in Tests/R_codes, and output data are in Tests/Outputs, ending with *.RData*.
- The *pde_solver_data* folder contains all codes and data for solving PDEs used in the paper. 
- The *Functions* folder contains all R functions used to build the candidate library and the recurrent adaptive lasso. 
- The *Figures* folder contains the results *.png* graphs used in the paper and ann *.R* file used for creating *.png* graphs. 

# Requirements
## Programming 
We write most scripts in R, but some functions are written in Python, so we use the R package `reticulate` to call Python functions in R. 
```
R>=4.1.x
Python>=3.9.x
Matlab>=2019
```

## R packages
```
library(pracma)
library(R.matlab)
library(glmnet)
library(dplyr)
library(signal)
library(orthopolynom)
library(reticulate)
library(signal)
library(Deriv)
library(Metrics)
library(purrr)
library(SMFilter)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
library(cowplot)
library(ggh4x)
library(stringr)
library(ggpubr)
library(scales)
library(gridExtra)
library(tibble)
library(dplyr)
library(plotly)
library(plot3D)
```

## Python Packages
```
numpy
scipy
scikit_learn
opencv_python
pandas
matplotlib
multiprocessing
```

When running test,s please check whether the working folder is correct. 
