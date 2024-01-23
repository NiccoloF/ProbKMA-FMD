# load library
library(ProbKMA.package)

# set seed
seed = 123
set.seed(seed)

# load data
data(sim_clusters)

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0
trials_elong = 201
c_max = 71
K = 3
c = 41

my_output = ProbKMA.package::probKMA(Y0,Y1,standardize=TRUE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                                     diss=diss,alpha=alpha,w=1,m=2,
                                     iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                                     iter4elong=1,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,
                                     deltaJk_elong=0.05,max_gap=max_gap,
                                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                                     return_options=TRUE,return_init=TRUE,worker_number=NULL)


# load previous R implementation source file
# warning: for comparing the results reload all the data and clean the envir
source(file = "../probKMA-FMD_functions.r")

true_output = probKMA(Y0,Y1,standardize=TRUE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                     diss=diss,alpha=alpha,w=1,m=2,
                     iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                     iter4elong=1,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,
                     deltaJk_elong=0.05,max_gap=max_gap,
                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                     return_options=TRUE,return_init=TRUE,worker_number=NULL)

# computational times comparison
system.time(ProbKMA.package::probKMA(Y0,Y1,standardize=TRUE,K=K,c=c,c_max=c_max,P0=NULL,S0=NULL,
                                     diss=diss,alpha=alpha,w=1,m=2,
                                     iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                                     iter4elong=1,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,
                                     deltaJk_elong=0.05,max_gap=max_gap,
                                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                                     return_options=TRUE,return_init=TRUE,worker_number=NULL))

