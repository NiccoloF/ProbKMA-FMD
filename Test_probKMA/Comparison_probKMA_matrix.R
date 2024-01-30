# load library
library(ProbKMA.package)

# set seed
seed = 123
set.seed(seed)

# load data
load(file = "../Test_probKMA/Y.RData")

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0
trials_elong = 200
c_max = 53
K = 2
c = 40 

my_output = ProbKMA.package::probKMA(Y$Y0,Y$Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                                     diss=diss,alpha=alpha,w=c(0.5,0.5),m=2,
                                     iter_max=20,stop_criterion='max',quantile=NULL,tol=1e-8,
                                     iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,
                                     deltaJk_elong=0.05,max_gap=max_gap,
                                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                                     return_options=TRUE,return_init=TRUE,worker_number=NULL)



# load previous R implementation source file
source(file = "../probKMA-FMD_functions.r")


# set seed again otherwise you get different P0 and S0
set.seed(seed)


true_output = probKMA(Y$Y0,Y$Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                      diss=diss,alpha=alpha,w=c(0.5,0.5),m=2,
                      iter_max=20,stop_criterion='max',quantile=NULL,tol=1e-8,
                      iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,deltaJk_elong=0.05,max_gap=max_gap,
                      iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                      return_options=TRUE,return_init=TRUE,worker_number=NULL)


# computational times comparison
system.time(ProbKMA.package::probKMA(Y$Y0,Y$Y1,standardize=FALSE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                                     diss=diss,alpha=alpha,w=c(0.5,0.5),m=2,
                                     iter_max=200,stop_criterion='max',quantile=NULL,tol=1e-8,
                                     iter4elong=2,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,deltaJk_elong=0.05,max_gap=max_gap,
                                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                                     return_options=TRUE,return_init=TRUE,worker_number=NULL))


###############################################################
############## find_candidate_motifs test part ################
###############################################################
rm(list = ls()) # clean the environment
load("../Test_probKMA/Y.RData") # load data
# set seed
seed = 123
set.seed(seed)

diss='d0_L2'
iter_max = 250
standardize = FALSE
alpha=0
max_gap = 0.2 # no gaps allowed
iter4elong = 5 # perform elongation
trials_elong = 100
c_max = 150 
### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(40, 50) # minimum motif lengths to try
n_init = 10 # number of random initializations to try
return_options = TRUE


# C++ version
system.time(ProbKMA.package::find_candidate_motifs(Y$Y0, NULL, K, c, n_init,
                                  name = '../Test_probKMA/results/our/matrix_data.1', names_var = 'x(t)',
                                  probKMA_options = list(c_max = c_max, standardize = standardize, iter_max = iter_max,
                                                         iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                         return_options = return_options, return_init = TRUE,
                                                         diss = diss, alpha = alpha),
                                  plot = FALSE, worker_number = NULL))







