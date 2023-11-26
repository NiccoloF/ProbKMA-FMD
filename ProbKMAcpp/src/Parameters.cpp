#include "Parameters.hpp"

 Parameters::Parameters(const Rcpp::List& params)
  {
    _K = params[0];
    _c = Rcpp::as<arma::ivec>(params[1]);
    _c_max = Rcpp::as<arma::ivec>(params[2]);
    _iter_max = params[3];
    _quantile = params[4];
    _stopCriterion = Rcpp::as<std::string>(params[5]);
    _tol = params[6];
    _iter4elong = params[7];
    _tol4elong = params[8];
    _max_elong = params[9];
    _trials_elong = params[10];
    _deltaJK_elong = params[11];
    _max_gap = params[12];
    _iter4clean = params[13];
    _tol4clean = params[14];
    _quantile4clean = params[15];
    _parallel = params[16];
  }


