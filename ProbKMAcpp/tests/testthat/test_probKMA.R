library(ProbKMAcpp)

context("ProbKMA")

test_that("Class ProbKMA works", {
  # seed for random initialization of P0 and S0
  seed = 1
  set.seed(seed)
  
  # set type of distance
  diss = 'd0_d1_L2' # try with d0_L2 d0_d1_L2 d1_L2
  
  # null matrix for random initialization
  P0= matrix() 
  S0= matrix() 
  
  params <- list(standardize=TRUE, K=2,c = 61,c_max = 71,iter_max = 1000,
                 quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
                 iter4elong = 1,tol4elong = 1e-3,max_elong = 0.5,
                 trials_elong = 201, deltaJK_elong = 0.05,max_gap = 0,iter4clean = 50,
                 tol4clean = 1e-4,
                 quantile4clean = 1/2,return_options = TRUE,
                 m = 2,w = 1,alpha = 0.5,seed = seed,
                 exe_print = FALSE,set_seed = TRUE, n_threads = 1)
  
  # check input data are correct
  a <- initialChecks(simulated200$Y0,simulated200$Y1,P0,S0,params,diss,seed)
  
  # take checked data and parameters  
  params <- a$Parameters
  data <- a$FuncData
  
  # create an object of the class ProbKMA 
  prok = new(ProbKMA,data$Y,params,data$P0,data$S0,"H1")
  
  # run the probKMA algorithm 
  output <- prok$probKMA_run()
  
  # check length output
  testthat::expect_equal(length(output), 32) 
  
  # compare output with prof
  testthat::expect_true( abs(sum(trueOutputSim200$BC_dist_iter - output$BC_dist_iter) ) < 1e-15)

})

