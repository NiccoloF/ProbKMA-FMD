#ifndef COMPUTE_JK_HPP
#define COMPUTE_JK_HPP
#include "RcppArmadillo.h"
#include "find_min_diss.h"
#include <ranges>
#include <algorithm>
#include <limits>
#include <omp.h>

double compute_Jk_rcpp(const arma::field<arma::mat> & v,
                       const arma::ivec & s_k,
                       const arma::vec & p_k,
                       const arma::field<arma::mat> & Y,
                       double alpha,
                       const arma::vec & w,
                       int m,
                       bool use0,
                       bool use1, 
                       double c_k, // in reality is an int
                       arma::vec keep_k); // in reality is an uvec

#endif //COMPUTE_JK_HPP