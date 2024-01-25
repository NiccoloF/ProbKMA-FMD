# ProbKMA-FMD

R code implementing **ProbKMA** (probabilistic K-means with local alignment) for local clustering of functional data and functional motif discovery, proposed in the paper [Probabilistic K-means with local alignment for clustering and motif discovery in functional data](https://doi.org/10.1080/10618600.2022.2156522), by Marzia A. Cremona and Francesca Chiaromonte. 

## Getting Started
The source code can be cloned or downloaded directly from Github.
An R studio project file is provided to open the project in RStudio. --> da fare 

### Prerequisites

The package is linked against **OpenMP**, the **BLAS** and **LAPACK** libraries in Makevars.

To install the package locally, **Rcpp** and **RcppArmadillo** packages need to be installed.\
Using install_github() will install the dependencies (Rcpp and RcppArmadillo) automatically.

### Installing

The package can be installed directly from Github but **Devtools** is required.

If devtools is not installed use the following comand to install it
```
install.packages('devtools') 
```
and then install the package
```
library(devtools)
install_github('NiccoloF/ProbKMA-FMD',ref = "main_3")
```


Usually the clag compiler doesn't support OpenMP. For this reason is available a release without -fopenmp.
```
library(devtools)
install_github('zitale/fdakmapp@v2.0.2.noomp')
```

Othewise check: https://clang-omp.github.io/.

### Example

```
res<-kmap(x=aneurisk65$x, y=aneurisk65$y, n_clust=2)
kmap_show_results(res,FALSE,FALSE)
```

## Documentation

The R documentation can be found in the main directory in fdakmapp.pdf.
The C++ documentation can be found at https://niccolof.github.io/ProbKMA-FMD/index.html.

### Code

Starting from R functions **ProbKMA** and **find_candidate_motifs** we have developed a C++ implementation.
- `ProbKMA`: central class that handles probabilistic k-means with local alignment to find candidate motifs
- `find_candidate_motifs`: run multiple times `ProbKMA` function with different K, c and initializations, with the aim to find a set of candidate motifs

## Functional motif discovery example
Functional motif discovery on simulated data: 20 curves embedding 2 functional motifs of length 60, each with 12 occurrences. 
- `len200_sd0.1.RData`: simulated curves
- `len200_sd0.1_simulated_curves_with_motifs.pdf`: plot of curves with true motifs
- `FMD_simulated_data.r`: script to run the example
- `results`: functional motif discovery results

## Probabilistic local clustering examples

### Berkley growth curves
Probabilitstic local clustering of the Berkley Growth Study dataset, provided within the R package `fda` and consisting of the heights of 39 boys and 54 girls recorded from age 1 to 18.
- `growth_smoothed.RData`: smoothed growth curves
- `growth_curves.pdf`: plot of growth curves and growth velocities
- `probKMA_growth.r`: script to run the example
- `results`: probabilistic local clustering results

### Italian Covid-19 excess mortality curves
Probabilitic local clustering of Covid-19 excess mortality rate curves (daily difference between 2020 deaths and average deaths in the period 2015-2019) in the 20 regions of Italy. These curves were estimated in the period from February 16, 2020 and April 30, 2020 using the mortality data (due to all causes) from the Italian Institute of Statistics ISTAT. Raw mortality data are available on [ISTAT website](https://www.istat.it/it/files/2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record-4giugno.zip).
- `istat_mortality_rates_smoothed.Rdata`: smoothed excess mortality rate curves
- `istat_mortality_rates_curves.pdf`: plot of excess mortality rate curves
- `probKMA_mortality_regions.r`: script to run the example
- `results`: probabilistic local clustering results
