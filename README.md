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
install_github('NiccoloF/ProbKMA-FMD/ProbKMAcpp',ref = "main_3")
```
Several configurations are possible :
- **main_3** : more performant package with a full C++ design layout for **ProbKMA**.
- **main_2** : full implementation of the **ProbKMA-FMD** package with C++ parallelism(for illustrative purposes only).
- **main_1** : full implementation of the **ProbKMA-FMD** package with only R parallelism.

Usually the clag compiler doesn't support **OpenMP**. For this reason is available a release without -fopenmp.
```
library(devtools)
install_github('NiccoloF/ProbKMA-FMD/ProbKMAcpp@v1.0.1.nopenmp')
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
               set_seed= FALSE, n_threads = 7)
data <- initialChecks(Y0,NULL,P0,S0,params,diss,seed)
params <- data$Parameters
data <- data$FuncData
prok <- new(ProbKMA,data$Y,params,data$P0,data$S0,"L2")
prok$probKMA_run()
```

## Documentation

The C++ documentation can be found at [https://niccolof.github.io/ProbKMA-FMD/index.html](https://niccolof.github.io/ProbKMA-FMD/).

### Code

Starting from R functions **ProbKMA** and **find_candidate_motifs** we have developed a full design layout C++ implementation.
- `ProbKMA`: central class that handles probabilistic k-means with local alignment to find candidate motifs
- `find_candidate_motifs`: run multiple times `ProbKMA` function with different number of clusters, minimum motif's length and initializations, with the aim to find a set of candidate motifs.

## Functional motif discovery example
Functional motif discovery on simulated data: 20 curves embedding 2 functional motifs of length 60, each with 12 occurrences. 
- `len200_sd0.1.RData`: simulated curves
- `len200_sd0.1_simulated_curves_with_motifs.pdf`: plot of curves with true motifs
- `FMD_simulated_data.r`: script to run the example
- `results`: functional motif discovery results

## Probabilistic local clustering example

### Berkley growth curves
Probabilitstic local clustering of the Berkley Growth Study dataset, provided within the R package `fda` and consisting of the heights of 39 boys and 54 girls recorded from age 1 to 18.
- `growth_smoothed.RData`: smoothed growth curves
- `growth_curves.pdf`: plot of growth curves and growth velocities
- `probKMA_growth.r`: script to run the example
- `results`: probabilistic local clustering results
