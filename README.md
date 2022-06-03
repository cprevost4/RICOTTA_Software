# MULTI-FRAME SUPER-RESOLUTION MRI USING COUPLED LOW-RANK TUCKER APPROXIMATION

Copyright (c) 2022 Clémence Prévost, Freddy Odille <br>
Contact: ```clemence.prevost@univ-lille.fr```

This software reproduces the results from the following:
```
@unpublished{prevost:hal-03617754,
  TITLE = {{MULTI-FRAME SUPER-RESOLUTION MRI USING COUPLED LOW-RANK TUCKER APPROXIMATION}},
  AUTHOR = {Pr{\'e}vost, Cl{\'e}mence and ODILLE, F},
  URL = {https://hal.archives-ouvertes.fr/hal-03617754},
  NOTE = {working paper or preprint},
  YEAR = {2022},
  MONTH = Mar,
  PDF = {https://hal.archives-ouvertes.fr/hal-03617754/file/IRM_Tucker.pdf},
  HAL_ID = {hal-03617754},
  HAL_VERSION = {v1},
}
```

<br><br>
[Link to the project](https://github.com/cprevost4/RICOTTA_Software)

## Content

 - /demos : contains demo files that produce tables and figures
 - /metrics : contains the metrics used to assess the quality of the reconstruction
 - /src : contains helpful files to run the demos

## Minimal requirements

 In order to run the demo file ```demo.m```, you will need to:
 - Download and install Tensorlab 3.0: https://www.tensorlab.net
 - Download the Sharpness Index toolbox: http://www.mi.parisdescartes.fr/~moisan/sharpness/
 - Beltrami primal-dual solver: https://math.montana.edu/dzosso/code/
 
 Please quote the corresponding papers if you decide to use these codes.

 ## How it works
 
 ### Generate coupled tensor model
 
 In this software, we use the "MRI" dataset of MATLAB. The low-resolution observations are generated from the super-resolution image with manually-specified degradation matrices.

 ### Run algorithms
 
 In ```reconstruction.m```, we showcase the performance of three algorithms:
  - RICOTTA
  - Block-RICOTTA (applies RICOTTA to corresponding subblocks of the observations)
  - RICOTTA without regularization + 3D-Beltrami regularization (to highlight the importance of the Tikhonov       regularization in the algorithm

The metrics and computation time are then displayed in a table.
Slices of the reference and reconstructions are plotted in a figure.

## Available demos

They are available in the ```/demos``` folder.
The table below summarized what does what

| Name                       | Content                                                        |
|----------------------------|----------------------------------------------------------------|
| ```reconstruction.m```     | Evaluates performance of the algorithms                        |
| ```choice_ranks.m```       | plots R-SNR, CC and RMSE as a function of the ranks            |
| ```choice_regul.m```       | plots R-SNR, CC and RMSE as a function of the regul. parameter |
| ```choice_weights.m```     | plots R-SNR, CC and RMSE as a function of the weights lambda   |
