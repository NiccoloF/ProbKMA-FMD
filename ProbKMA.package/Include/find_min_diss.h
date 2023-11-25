#ifndef FIND_MIN_DISS_H
#define FIND_MIN_DISS_H
#include "RcppArmadillo.h"
#include <ranges>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <limits>
#include "utilities.h"


Rcpp::List find_shift_warp_min(const Rcpp::List & Y, 
                               const Rcpp::List & V_new,
                               const arma::vec & w,
                               const arma::ivec & c_k,
                               unsigned int K,
                               unsigned int d,
                               double max_gap,
                               double alpha,
                               bool use0,
                               bool use1);

arma::vec find_diss_aligned_rcpp(const arma::mat &y0,
                                 const arma::mat &y1,
                                 const arma::mat &v0,
                                 const arma::mat &v1, 
                                 const arma::vec & w, 
                                 double alpha,
                                 bool aligned,
                                 unsigned int d,
                                 bool use0,
                                 bool use1);

arma::vec find_diss_rcpp(const arma::mat &y0,
                         const arma::mat &y1,
                         const arma::mat &v0,
                         const arma::mat &v1,
                         const arma::vec & w, 
                         double alpha, unsigned int c_k,
                         unsigned int d,bool use0,bool use1);

#endif // FIND_MIN_DISS_H