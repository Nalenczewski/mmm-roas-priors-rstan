# Geo-Hierarchical Media Mix Model with ROAS Priors in Rstan
This code is an implementation of the simulation in the Google paper [Media Mix Model Calibration with Bayesian Priors](https://research.google/pubs/media-mix-model-calibration-with-bayesian-priors/), which shows how to directly incorporate ROAS values obtained through experimentation into a Bayesian media mix model. For a full overview, check out the [Medium article](https://medium.com/@nicklenczewski/building-a-media-mix-model-with-roas-priors-in-rstan-a25321a87724) or [YouTube tutorial](https://youtu.be/hM1jHpvI2jU).

## Getting Started
1. Install [R](https://cran.rstudio.com/).
2. Install [Rstudio](https://posit.co/download/rstudio-desktop/).
3. Install the following packages:
```
install.packages(c("dplyr", "ggplot2", "tidyverse", "cowplot", "rstan"))
``` 
If you're having trouble installing rstan, reference its [repository](https://github.com/stan-dev/rstan/wiki/Rstan-Getting-Started).

4. Clone this repository
```
git clone https://github.com/Nalenczewski/mmm-roas-priors-rstan.git
```
## Contact
nalenczewski@gmail.com, Nick Lenczewski, Data Scientist
