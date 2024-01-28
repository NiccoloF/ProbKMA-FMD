# ProbKMA-FMD

C++ implementation of **ProbKMA** (probabilistic K-means with local alignment) for local clustering of functional data and functional motif discovery, proposed in the paper [Probabilistic K-means with local alignment for clustering and motif discovery in functional data](https://doi.org/10.1080/10618600.2022.2156522), by Marzia A. Cremona and Francesca Chiaromonte. 

## Getting Started
The source code can be cloned or downloaded directly from Github.

### Prerequisites

The package is linked against **OpenMP**, the **BLAS** and **LAPACK** libraries in Makevars.

To install the package locally, **Rcpp** and **RcppArmadillo** packages need to be installed.\
Using ```install_github()``` will install the dependencies (**Rcpp** and **RcppArmadillo**) automatically.

### Installing

The package can be installed directly from Github but **Devtools** is required.

If devtools is not installed use the following comand to install it
```
install.packages('devtools') 
```
and then install the package
```
library(devtools)
install_github('NiccoloF/ProbKMA-FMD',ref = "main_3",subdir = 'ProbKMAcpp')
```
Several configurations are possible :
- **main_3** : more performant package with a full C++ design layout for **ProbKMA**.
- **main_2** : full implementation of the **ProbKMA-FMD** package with C++ parallelism(for illustrative purposes only).
- **main_1** : full implementation of the **ProbKMA-FMD** package with only R parallelism.

Usually the clag compiler doesn't support **OpenMP**. For this reason is available a release without -fopenmp.
```
library(devtools)
install_github('NiccoloF/ProbKMA-FMD',ref = 'v1.0.1.nopenmp',subdir = 'ProbKMAcpp')
```

Othewise check: https://clang-omp.github.io/.

### Example

```
params <- list(standardize=FALSE, K=2,c = c_min,c_max = c_min,iter_max = 1000,
               quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
               iter4elong = 2,tol4elong = 1e-3,max_elong = 0.5,
               trials_elong = 10,
               deltaJK_elong = 0.05,max_gap = 0.2,iter4clean = 1000,
               tol4clean = 1e-4,
               quantile4clean = 1/2,return_options = TRUE,
               m = 2,w = 1, alpha = 0.0,seed = seed,exe_print = TRUE, 
               set_seed= FALSE, n_threads = 7) # Inizialize parameters
data <- initialChecks(Y0,NULL,P0,S0,params,diss,seed) # Y0 and eventually Y1 are functional data 
params <- data$Parameters # Get the updated parameters for the problem
data <- data$FuncData # Get the updated functional data 
prok <- new(ProbKMA,data$Y,params,data$P0,data$S0,"L2") # Initialization of the main class that handles ProbKMA
prok$probKMA_run() # Run ProbKMA algorithm 
```

### Code

Starting from R functions **ProbKMA** and **find_candidate_motifs** we have developed a full design layout C++ implementation.
- `ProbKMA`: central class that handles probabilistic k-means with local alignment to find candidate motifs
- `find_candidate_motifs`: run multiple times `ProbKMA` function with different number of clusters, minimum motif's length and initializations, with the aim to find a set of candidate motifs.

## Functional motif discovery example
Within the **Test comparisons folder** is the **R** script **Comparisons_vectorial_data.R** in which an example is analyzed on the 
discovery of functional motifs on simulated data with 20 curves incorporating 2 functional motifs of length 60, each with 12 occurrences.
In particular, the following is found:
- `simulated200`: simulated curves already built into the package. 
- `len200_sd0.1_simulated_curves_with_motifs.pdf`: plot of curves with true motifs.
- `Comparisons_vectorial_data.R`: script to run the example. \
Finally, the differences in computational terms between the previous motif discovery algorithm written entirely in **R** and the new **C++** version can be analyzed.<br>
<br> Similar to the previous case, within the **R** script **Comparisons_matrices_data.R**, an example can be run using multivariate data. 
In particular it can be found: 
- `Y.Rdata`: functional data to be loaded.
- `Comparisons_matrices_data.R`: script to run the example in the multivariate case.\
In particular, in this last example, one can appreciate the significant computational differences between the two code versions. 



## Probabilistic local clustering example

### Italian Covid-19 excess mortality curves
Probabilitic local clustering of Covid-19 excess mortality rate curves (daily difference between 2020 deaths and average deaths in the period 2015-2019) in the 20 regions of Italy. These curves were estimated in the period from February 16, 2020 and April 30, 2020 using the mortality data (due to all causes) from the Italian Institute of Statistics ISTAT. Raw mortality data are available on [ISTAT website](https://www.istat.it/it/files/2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record-4giugno.zip).
- `istat_mortality_rates_smoothed.Rdata`: smoothed excess mortality rate curves
- `istat_mortality_rates_curves.pdf`: plot of excess mortality rate curves
- `probKMA_mortality_regions.r`: script to run the example
- `results`: probabilistic local clustering results

## Documentation

The C++ documentation can be found at [https://niccolof.github.io/ProbKMA-FMD/index.html](https://niccolof.github.io/ProbKMA-FMD/).
