This folder contains the scripts used for the paper [Rodero, C., Strocchi, M., Lee, A. W., Rinaldi, C. A., Vigmond, E. J., Plank, G., ... & Niederer, S. A. (2022). Impact of anatomical reverse remodelling in the design of optimal quadripolar pacing leads: A computational study. _Computers in biology and medicine_, 140, 105073](https://www.sciencedirect.com/science/article/pii/S0010482521008672). They are divided according to the language.

-In the Matlab folder the scripts were mainly used to create AHA maps. There's probably by now a more updated version of it.
-The python scripts contain mainly the functions to run the EP simulations.
-The sh (shell) scripts were used mainly for the sensitivity analysis. To see a more detailed explaination of the steps to recreate it, check ```lifesaver.sh```
-The R scripts were the main tool to pre- and post- process the data. The steps can be found in ```multipolar_pipeline.R```. Each functions should have a brief explanation of what it does if it's not obvious. Some functions are obsolete, when a version parameter is needed use the most recent one (version 4).
