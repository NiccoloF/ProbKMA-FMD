#ifndef ELONGATE_MOTIFS_HPP
#define ELONGATE_MOTIFS_HPP
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>
#include <variant>
#include <omp.h>
#include "utilities.h"
#include "compute_jk.h"

void elongation_rcpp(arma::field<arma::mat> & V_new, 
                     std::vector<arma::uvec> & V_dom,  
                     arma::imat & S_k, 
                     const arma::vec & p_k, 
                     const arma::ivec& len_elong_k, 
                     const arma::uvec& keep_k,  
                     double c, 
                     bool use0,
                     bool use1,
                     const arma::vec& w, 
                     double alpha, double max_gap,  
                     const arma::field<arma::mat> Y, 
                     int m, 
                     double deltaJk_elong,
                     const unsigned int index);

void elongate_motifs(Rcpp::List & V_new,
                     Rcpp::List & V_dom,
                     Rcpp::List & S_k,
                     const Rcpp::List & P_k,
                     const Rcpp::List & Y,
                     const arma::vec & w, 
                     int m, 
                     bool use0,
                     bool use1,
                     double alpha,
                     const arma::ivec & c,
                     const arma::ivec & c_max, // is a vector to be understood
                     double max_elong, 
                     double deltaJk_elong,
                     int trials_elong,
                     const arma::mat & D,
                     unsigned int K,
                     double max_gap);


#endif //ELONGATE_MOTIFS_HPP