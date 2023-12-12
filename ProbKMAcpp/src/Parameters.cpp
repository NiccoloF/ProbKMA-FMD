#include "Parameters.hpp"

 Parameters::Parameters(const Rcpp::List& params)
  {
    _standardize = params["standardize"];
    _K = params["K"];
    _c = Rcpp::as<arma::ivec>(params["c"]);
    _c_max = Rcpp::as<arma::ivec>(params["c_max"]);
    _iter_max = params["iter_max"];
    _quantile = params["quantile"];
    _prob = params["prob"];
    _stopCriterion = Rcpp::as<std::string>(params["stopCriterion"]);
    _m = params["m"];
    _w = Rcpp::as<arma::vec>(params["w"]);
    _alpha = params["alpha"];
    _tol = params["tol"];
    _iter4elong = params["iter4elong"];
    _tol4elong = params["tol4elong"];
    _max_elong = params["max_elong"];
    _trials_elong = params["trials_elong"];
    _deltaJK_elong = params["deltaJK_elong"];
    _max_gap = params["max_gap"];
    _iter4clean = params["iter4clean"];
    _tol4clean = params["tol4clean"];
    _quantile4clean = params["quantile4clean"];
    _return_options = params["return_options"];
  }


