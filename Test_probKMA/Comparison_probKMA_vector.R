# load library
library(ProbKMA.package)

# set seed
seed = 123
set.seed(seed)

# load data
data(sim_motifs)

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0
trials_elong = 201
c_max = 71
K = 3
c = 61
P0 = NULL
S0 = NULL

my_output = ProbKMA.package::probKMA(Y0,Y1,standardize=TRUE,K=K,c=c,c_max=c_max,P0=P0,S0=S0,
                                     diss=diss,alpha=alpha,w=1,m=2,
                                     iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                                     iter4elong=1,tol4elong=1e-3,max_elong=0.5,trials_elong=trials_elong,
                                     deltaJk_elong=0.05,max_gap=max_gap,
                                     iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                                     return_options=TRUE,return_init=TRUE,worker_number=1)


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



#########################################################
##### computational times comparison ####################
#########################################################

# ProbKMA
rm(list = ls()) # clean the environment
# set seed
seed = 123
set.seed(seed)
# load data
data(sim_motifs)

diss = 'd0_d1_L2' 
P0= matrix() 
S0= matrix()
standardize=TRUE
K=2
c = 41
c_max = 71
iter_max = 1000
quantile = 0.25
stopCriterion = 'max'
tol = 1e-8
iter4elong = 1
tol4elong = 1e-3
max_elong = 0.5
trials_elong = 201 
deltaJK_elong = 0.05
max_gap = 0
iter4clean = 50
tol4clean = 1e-4
quantile4clean = 1/2
return_options = TRUE
m = 2
w = 1
alpha = 0.5


system.time(ProbKMA.package::probKMA(Y0,Y1,standardize=standardize,K=K,c=c,c_max=c_max,P0=P0,S0=S0,
                                     diss=diss,alpha=alpha,w=w,m=m,
                                     iter_max=iter_max,stop_criterion='max',quantile=quantile,tol=tol,
                                     iter4elong=iter4elong,tol4elong=tol4elong,
                                     max_elong=max_elong,trials_elong=trials_elong,
                                     deltaJk_elong=deltaJK_elong,max_gap=max_gap,
                                     iter4clean=iter4clean,tol4clean=tol4clean,quantile4clean=quantile4clean,
                                     return_options=return_options,return_init=TRUE,worker_number=NULL))

#############################################
####### Find candidate Motifs ###############
#############################################

rm(list = ls()) # clean the environment
# load data
data(sim_motifs)

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0 # no gaps allowed
iter4elong = 1 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 71 # maximum motif length 70
#run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(61, 51, 41) # minimum motif lengths to try
n_init = 10 # number of random initializations to try
iter_max = 1000
return_options = TRUE

system.time(ProbKMA.package::find_candidate_motifs(Y0,Y1, K, c, n_init,
                                                   name = '../Test_probKMA/results/our/len200_sd0.1', names_var = 'x(t)',
                                                   probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = iter_max,
                                                                           iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                           return_options = return_options, return_init = TRUE,
                                                                           diss = diss, alpha = alpha),
                                                   plot = FALSE,worker_number = NULL))





