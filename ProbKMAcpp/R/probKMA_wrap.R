#' @title probKMA_wrap
#'
#' @description R wrap for function probKMA used to parallelize find_candidate_motifs
#' 
#' @param Y0 ...
#' 
#' @return outputs of probKMA
#'
#' @author Riccardo Lazzarini Niccolò Feresini
#' @export
probKMA_wrap <- function(Y0 = NULL,Y1 = NULL,P0 = matrix(),S0 = matrix(),
                         standardize= FALSE,c_max = Inf,iter_max = 1000,
                         iter4elong = 10,trials_elong = 10,return_options = TRUE,
                         alpha = 0,max_gap = 0.2,quantile = 0.25, stopCriterion = 'max', 
                         tol = 1e-8, tol4elong = 1e-3, max_elong = 0.5, deltaJK_elong = 0.05, 
                         iter4clean = 50, tol4clean = 1e-4,m = 2,w = 1, seed = 1, 
                         K = 2, c = 40, quantile4clean = 1/K, exe_print = FALSE,
                         set_seed = FALSE,n_threads = 7,diss = 'd0_2'){
  
  params = list(standardize=standardize,c_max = c_max,iter_max = iter_max,
                iter4elong = iter4elong,trials_elong = trials_elong,
                return_options = return_options, alpha = alpha,
                max_gap = max_gap,quantile = quantile, 
                stopCriterion = stopCriterion, tol = tol, 
                tol4elong = tol4elong, max_elong = max_elong, 
                deltaJK_elong = deltaJK_elong, iter4clean = iter4clean, 
                tol4clean = tol4clean,quantile4clean = quantile4clean, 
                m = m, w = w, seed = seed, K = K, c = c, exe_print = exe_print,
                set_seed = set_seed,n_threads = n_threads) 
  
  checked_data <- initialChecks(Y0,Y1,P0,S0,params,diss,seed)
  
  params <- checked_data$Parameters
  
  data <- checked_data$FuncData
  
  if ( alpha == 0 ||  alpha == 1)
  {
    string  = "L2"
  } 
  else 
  {
    string = "H1"
  }
  
  prok = new(ProbKMA,data$Y,params,data$P0,data$S0,string)
  
  rm(params)
  
  probKMA_results_1 = list(Y0 = data$Y$Y0,Y1 = data$Y$Y1,
                           diss = diss, w = w, alpha = alpha)
  
  rm(data)
  
  probKMA_results_2 = prok$probKMA_run() 
  
  rm(prok)
  
  return(c(probKMA_results_1,probKMA_results_2))
}
