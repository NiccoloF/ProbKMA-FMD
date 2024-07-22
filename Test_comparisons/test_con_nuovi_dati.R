library(ProbKMAcpp)

# seed for random initialization of P0 and S0
seed = 1
set.seed(1)

# set type of distance
diss = 'd0_d1_L2' # try with d0_L2 d0_d1_L2 d1_L2

# null matrix for random initialization
P0= matrix() 
S0= matrix() 

params <- list(standardize=FALSE, K=2,c = 40,c_max = 70,iter_max = 1000, 
               quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
               iter4elong = 1,tol4elong = 1e-3,max_elong = 0.5, 
               trials_elong = 201, deltaJK_elong = 0.05,max_gap = 0,iter4clean = 50,
               tol4clean = 1e-4,
               quantile4clean = 1/2,return_options = TRUE,
               m = 2,w = 1,alpha = 0.5,seed = seed,exe_print = TRUE, 
               set_seed = TRUE, n_threads = 7, transformed = FALSE) 

# check input data are correct
# per usare v_init_test load(v_init_tes.RData), da gestire il caso in cui v_init è null
#a <- initialChecks(Y0 = Y0, Y1 = Y1,P0 = P0,S0 = S0,params = params,diss = diss, v_init = v_init_test) #v_init_test
# caso in cui v_init è nullo
a <- initialChecks(Y0 = Y0, Y1 = Y1,P0 = P0,S0 = S0,params = params,diss = diss, v_init = NULL)

# take checked data and parameters  
params <- a$Parameters
data <- a$FuncData

# create an object of the class ProbKMA 
#prok = new(ProbKMA,data$Y,params,data$P0,data$S0,"H1",data$v_init) #data$v_init to be removed if testing without v_init_test
# caso in cui v_init è nullo
prok = new(ProbKMA,data$Y,params,data$P0,data$S0,"H1")

# run the probKMA algorithm 
output <- prok$probKMA_run()

sil <- prok$compute_silhouette(TRUE) # try TRUE


# for the plot add other info
output <- c(output, list(Y0 = data$Y$Y0,Y1 = data$Y$Y1,
                         diss = diss,w = params$w,alpha = params$alpha))

# plot the results of probKMA
pdf(paste0('our_plot_vectorial','.pdf'),width=20,height=10)
probKMA_plot(output)
dev.off()

# comparison with previous implementation
source(file ="../Test_comparisons/previous_ProbKMA.R")

true_output <- probKMA(Y0=Y0,Y1=Y1,standardize=params$standardize,
                       transformed = FALSE,K=params$K,c=params$c[1],c_max=params$c_max,
                       P0=data$P0,S0=data$S0,
                       diss=diss,alpha=params$alpha,w=params$w,m=params$m,v_init = NULL, # v_init_test
                       iter_max=params$iter_max,
                       stop_criterion=params$stopCriterion,
                       quantile=params$quantile,tol=params$tol,iter4elong=params$iter4elong,
                       tol4elong=params$tol4elong,max_elong=params$max_elong,
                       trials_elong=params$trials_elong,deltaJk_elong=params$deltaJK_elong,
                       max_gap=params$max_gap,params$iter4clean,params$tol4clean,
                       params$quantile4clean,params$return_options,return_init=TRUE,
                       worker_number=NULL)


true_silhouette <- probKMA_silhouette(true_output,align = TRUE, plot = FALSE) # try TRUE

# true plot
pdf(paste0('true_plot_vectorial','.pdf'),width=20,height=10)
probKMA_plot(true_output)
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


# test wrap 
# lanciare fino a riga 14
arguments =   list(Y0 = simulated200$Y0, 
                   Y1 = simulated200$Y1,
                   P0 = matrix(), 
                   S0 = matrix(),
                   standardize=params$standardize,
                   c_max = params$c_max,
                   iter_max = params$iter_max,
                   iter4elong = params$iter4elong,
                   trials_elong = params$trials_elong,
                   return_options = params$return_options,
                   alpha = params$alpha,
                   max_gap = params$max_gap,
                   diss = diss,
                   quantile = params$quantile, 
                   stopCriterion = params$stopCriterion, 
                   tol = params$tol, 
                   tol4elong = params$tol4elong, 
                   max_elong = params$max_elong, 
                   deltaJK_elong = params$deltaJK_elong, 
                   iter4clean = params$iter4clean, 
                   tol4clean = params$tol4clean,
                   quantile4clean = params$quantile4clean, 
                   m = params$m,
                   w = params$w, 
                   seed = seed, 
                   K = params$K, 
                   c = params$c[1],
                   exe_print = TRUE,
                   set_seed = TRUE,
                   n_threads = 7,
                   transformed = FALSE, # TRUE or FALSE
                   v_init = v_init_test, # v_init_test or NULL
                   silhouette = TRUE) # TRUE or FALSE

results = do.call(probKMA_wrap,arguments)
probKMA_silhouette_plot(probKMA_results = results[[1]],silhouette_results = results[[2]])

# my plot
pdf(paste0('my_plot_vectorial','.pdf'),width=20,height=10)
probKMA_plot(results[[1]])
dev.off()



# test find candidate
library(ProbKMAcpp)
library(parallel)
diss = 'd0_d1_L2'
alpha = 0.9
max_gap = 0 # no gaps allowed
iter4elong =1 #5-10 # perform elongation
trials_elong =150 # try all possible elongations 60-70
c_max = 150 # maximum motif length 60-70
initializations=10
n_init = initializations # number of random initializations to try
K = c(2, 3) # number of clusters to try
c=c(40,50,60)
load("../Test_comparisons/V_init.RData")

find_candidate_motifs_results = ProbKMAcpp::find_candidate_motifs(simulated200$Y0,simulated200$Y1,K,c,n_init,
                                                                  name = '../Test_comparisons/results/our/len200_sd0.1', names_var = 'x(t)',
                                                                  probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                                         iter4elong = iter4elong, trials_elong = trials_elong, 
                                                                                         max_gap = max_gap,
                                                                                         return_options = TRUE, return_init = TRUE,
                                                                                         diss = diss, alpha = alpha),
                                                                  plot = TRUE,exe_print = FALSE,set_seed = FALSE, 
                                                                  V_init = V_init, transformed = TRUE)


