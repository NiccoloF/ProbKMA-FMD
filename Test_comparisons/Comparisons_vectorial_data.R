library(ProbKMAcpp)

# seed for random initialization of P0 and S0
seed = 1
set.seed(seed)

# set type of distance
diss = 'd0_d1_L2' # try with d0_L2 d0_d1_L2 d1_L2

# null matrix for random initialization
P0= matrix() 
S0= matrix() 

params <- list(standardize=TRUE, K=2,c = 41,c_max = 71,iter_max = 1000, 
               quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
               iter4elong = 1,tol4elong = 1e-3,max_elong = 0.5, 
               trials_elong = 201, deltaJK_elong = 0.05,max_gap = 0,iter4clean = 50,
               tol4clean = 1e-4,
               quantile4clean = 1/2,return_options = TRUE,
               m = 2,w = 1,alpha = 0.5,seed = seed,exe_print = FALSE, #TRUE
               set_seed = TRUE, n_threads = 3) 

# check input data are correct
a <- initialChecks(simulated200$Y0,simulated200$Y1,P0,S0,params,diss,seed)

# take checked data and parameters  
params <- a$Parameters
data <- a$FuncData

# create an object of the class ProbKMA 
prok = new(ProbKMA,data$Y,params,data$P0,data$S0,"H1")

# run the probKMA algorithm 
output <- prok$probKMA_run()

# for the plot add other info
output <- c(output, list(Y0 = data$Y$Y0,Y1 = data$Y$Y1,
                         diss = diss,w = params$w,alpha = params$alpha))

# plot the results of probKMA
pdf(paste0('our_plot_vectorial','.pdf'),width=20,height=10)
probKMA_plot(output,cleaned = TRUE)
dev.off()

# comparison with previous implementation
source(file ="../Test_comparisons/previous_ProbKMA.R") # @TODO: load using the library

true_output <- probKMA(Y0=simulated200$Y0,Y1=simulated200$Y1,standardize=params$standardize,K=params$K,c=params$c,c_max=params$c_max,
                       P0=data$P0,S0=data$S0,
                       diss=diss,alpha=params$alpha,w=params$w,m=params$m,iter_max=params$iter_max,
                       stop_criterion=params$stopCriterion,
                       quantile=params$quantile,tol=params$tol,iter4elong=params$iter4elong,
                       tol4elong=params$tol4elong,max_elong=params$max_elong,
                       trials_elong=params$trials_elong,deltaJk_elong=params$deltaJK_elong,
                       max_gap=params$max_gap,params$iter4clean,params$tol4clean,
                       params$quantile4clean,params$return_options,TRUE,NULL)

# true plot
pdf(paste0('true_plot_vectorial','.pdf'),width=20,height=10)
probKMA_plot(true_output,cleaned = TRUE)
dev.off()

#############################################################
########### computational time comparison probKMA ###########
#############################################################

# reinitialize motifs for a new run
prok$reinit_motifs(params$c,ncol(simulated200$Y0[[1]]))

# computational time c++ imp.
system.time(prok$probKMA_run())

# computational time R imp.
system.time(probKMA(Y0=simulated200$Y0,Y1=simulated200$Y1,standardize=params$standardize,K=params$K,c=params$c,c_max=params$c_max,
                    P0=data$P0,S0=data$S0,
                    diss=diss,alpha=params$alpha,w=params$w,m=params$m,iter_max=params$iter_max,
                    stop_criterion=params$stopCriterion,
                    quantile=params$quantile,tol=params$tol,iter4elong=params$iter4elong,
                    tol4elong=params$tol4elong,max_elong=params$max_elong,
                    trials_elong=params$trials_elong,deltaJk_elong=params$deltaJK_elong,
                    max_gap=params$max_gap,params$iter4clean,params$tol4clean,
                    params$quantile4clean,params$return_options,TRUE,NULL))


###############################################################
############## find_candidate_motifs test part ################
###############################################################

library(ProbKMAcpp)
library(parallel)

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


find_candidate_motifs_results = ProbKMAcpp::find_candidate_motifs(simulated200$Y0,simulated200$Y1,K,c,n_init,
                                                                  name = '../Test_comparisons/results/our/len200_sd0.1', names_var = 'x(t)',
                                                                  probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                                         iter4elong = iter4elong, trials_elong = trials_elong, 
                                                                                         max_gap = max_gap,
                                                                                         return_options = TRUE, return_init = TRUE,
                                                                                         diss = diss, alpha = alpha),
                                                                  plot = TRUE,exe_print = FALSE,set_seed = FALSE)

# time c++:
system.time(ProbKMAcpp::find_candidate_motifs(simulated200$Y0, simulated200$Y1, K, c, n_init,
                                              name = '../Test_comparisons/results/our/len200_sd0.1', names_var = 'x(t)',
                                              probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                     iter4elong = iter4elong, trials_elong = trials_elong,
                                                                     max_gap = max_gap,
                                                                     return_options = TRUE, return_init = TRUE,
                                                                     diss = diss, alpha = alpha),
                                              plot = FALSE,exe_print=FALSE,set_seed = FALSE))



# comparison with previous implementation
source(file ="../Test_comparisons/previous_ProbKMA.R") 

true_find_candidate_motifs_results = find_candidate_motifs(simulated200$Y0, simulated200$Y1, K, c, n_init,
                                                           name = '../Test_comparisons/results/prof/len200_sd0.1', names_var = 'x(t)',
                                                           probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                                  iter4elong = iter4elong, trials_elong = trials_elong,
                                                                                  max_gap = max_gap,
                                                                                  return_options = TRUE, return_init = TRUE,
                                                                                  diss = diss, alpha = alpha),
                                                           plot = TRUE, worker_number = NULL)


# time R previous imp:
system.time(find_candidate_motifs(simulated200$Y0, simulated200$Y1, K, c, n_init,
                                  name = '../Test_comparisons/results/prof/len200_sd0.1', names_var = 'x(t)',
                                  probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                         iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                         return_options = TRUE, return_init = TRUE,
                                                         diss = diss, alpha = alpha),
                                  plot = FALSE, worker_number = NULL))



