# set the directory of the package, otherwise initialChecks does not work
setwd("../ProbKMAcpp")
devtools::load_all()

# seed for random initialization of P0 and S0
seed = 1
set.seed(seed)

# set type of distance
diss = 'd0_d1_L2' # try with d0_L2 d0_d1_L2 d1_L2

# null matrix for random initialization
P0= matrix() 
S0= matrix() 

load("../Test_comparisons/Y.RData")

params <- list(standardize=TRUE, K=2,c = 40,c_max = 53,iter_max = 2,
               quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
               iter4elong = 2,tol4elong = 10,max_elong = 0.5,
               trials_elong = 30,
               deltaJK_elong = 0.05,max_gap = 0,iter4clean = 2,
               tol4clean = 1,
               quantile4clean = 1/2,return_options = TRUE,
               m = 2,w = c(0.5,0.5),alpha = 0.5,seed = seed,exe_print = TRUE, set_seed= TRUE)

# checks the parameters
a <- ProbKMAcpp::initialChecks(Y$Y0,Y$Y1,P0,S0,params,diss,seed)
params <- a$Parameters
data <- a$FuncData

# initialize the class ProbKMA
prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"H1")

# run ProbKMA algo written in C++
output <- prok$probKMA_run()

# comparison with previous implementation
source(file ="../Test_comparisons/previous_ProbKMA.R") # @TODO: load using the library

true_output <- probKMA(Y0=Y$Y0,Y1=Y$Y1,standardize=params$standardize,K=params$K,c=params$c,c_max=params$c_max,
                       P0=data$P0,S0=data$S0,
                       diss=diss,alpha=params$alpha,w=params$w,m=params$m,iter_max=params$iter_max,
                       stop_criterion=params$stopCriterion,
                       quantile=params$quantile,tol=params$tol,iter4elong=params$iter4elong,
                       tol4elong=params$tol4elong,max_elong=params$max_elong,
                       trials_elong=params$trials_elong,deltaJk_elong=params$deltaJK_elong,
                       max_gap=params$max_gap,params$iter4clean,params$tol4clean,
                       params$quantile4clean,params$return_options,TRUE,NULL)


###############################################################
############## find_candidate_motifs test part ################
###############################################################

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0 # no gaps allowed
iter4elong = 20 # perform elongation
trials_elong = 30 
c_max = 53 
### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(40, 45, 50) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

find_candidate_motifs_results = ProbKMAcpp::find_candidate_motifs(Y$Y0, Y$Y1, K, c, n_init,
                                                                  name = '../Test_comparisons/results/our/matrix_data.1', names_var = 'x(t)',
                                                                  probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 100,
                                                                                         iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                                         return_options = TRUE, return_init = TRUE,
                                                                                         diss = diss, alpha = alpha),
                                                                  plot = FALSE,exe_print=TRUE)




true_find_candidate_motifs_results = find_candidate_motifs(Y$Y0, Y$Y1, K, c, n_init,
                                                           name = '../Test_comparisons/results/prof/matrix_data.1', names_var = 'x(t)',
                                                           probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 100,
                                                                                  iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                                  return_options = TRUE, return_init = TRUE,
                                                                                  diss = diss, alpha = alpha),
                                                           plot = TRUE, worker_number = NULL)




