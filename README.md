# "A unifying framework for flexible excess hazard modelling with applications in cancer epidemiology"

### Alessia Eletti, Giampiero Marra, Manuela Quaresma, Rosalba Radice, and [F. Javier Rubio](https://sites.google.com/site/fjavierrubio67/)

This repository contains two data examples of link-based additive net survival models based on the paper 

> Eletti, A., Marra, G., Quaresma, M., Radice, R., and Rubio, F.J. (2021+). A unifying framework for flexible excess hazard modelling with applications in cancer epidemiology, *under review*.

The models are fitted using the R package `GJRM`.

1. The first example analyses the Simulacrum dataset (https://simulacrum.healthdatainsight.org.uk/) and focuses on differences in lung cancer survival according to different deprivation categories.

The R code and outputs for this example are also available in RPubs at: 

[A unifying framework for flexible excess hazard modelling with applications in cancer epidemiology: The Simulacrum](https://rpubs.com/FJRubio/SimulacrumLung)


2. The second example analyses the `LeukSurv` dataset from the R package `spBayesSurv` and provides additional insight on how to model spatial components. Note, however, that we do not have access to life tables from those times based on the Townsend score so for the purpose of the example we used generic values of the population hazards based on more recent life tables.

The R code and outputs for this example are also available in RPubs at: 

[A unifying framework for flexible excess hazard modelling with applications in cancer epidemiology: Leukemia data](https://rpubs.com/FJRubio/GJRMLeukSurv)


## Disclaimer
The data provided here and the numerical results represent numerical illustrations of the methods proposed in "A unifying framework for flexible excess hazard modelling with applications in cancer epidemiology". These data and results should not be used to make epidemiological claims. 

"Data for this study/project/report used artificial data from the Simulacrum, a synthetic dataset developed by Health Data Insight CiC derived from anonymous cancer data provided by the National Cancer Registration and Analysis Service, which is part of Public Health England."
