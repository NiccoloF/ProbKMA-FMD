library(ProbKMAcpp)
# load the package to be sure everything works, sometime without loading initialChecks doesn't work
devtools::load_all()

# data for find_candidate_motifs
diss = 'd0_d1_L2' # type of distance
alpha = 0.5 # sobolev-type distance
max_gap = 0 # no gaps allowed
iter4elong = 1 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 71 # maximum motif length 70
### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(61, 51, 41) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

# if you want to really tests find_candidate_motifs remove all the elements in the folder results
find_candidate_motifs_results = ProbKMAcpp::find_candidate_motifs(simulated200$Y0, simulated200$Y1, 
                                                                  K, c, n_init,
                                                                  name = './results/len200_sd0.1', names_var = 'x(t)',
                                                                  probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                                         iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                                         return_options = TRUE, return_init = TRUE,
                                                                                         diss = diss, alpha = alpha),
                                                                  plot = FALSE,set_seed = TRUE)


# tests
expect_equal(length(find_candidate_motifs_results),6)

expect_true(sum(colSums(find_candidate_motifs_results$silhouette_average_sd[[1]][[1]] 
                        - true_find_candidate_motifs_results$silhouette_average_sd[[1]][[1]])) == 0)
