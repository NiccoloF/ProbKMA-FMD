#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <string>
#include "RcppArmadillo.h"


struct Parameters
{ 
    /*
     *bool standardize;
     *unsigned int K,
     *const arma::ivec& c,
     *const arma::ivec& c_max,
     *unsigned int iter_max,
     *double quantile,
     *const std::string& stopCriterion,
     *double tol,
     *unsigned int iter4elong,
     *double tol4elong,
     *double max_elong,
     *unsigned int trials_elong,
     *double deltaJK_elong,
     *double max_gap,
     *unsigned int iter4clean,
     *double tol4clean,
     *double quantile4clean,
     *long long _seed;
     *bool return_options;
     */
    
    Parameters() = delete;
    Parameters(const Parameters&) = default;
    Parameters(const Rcpp::List& params);
    
    bool _standardize;
    unsigned int _K;
    arma::ivec _c;
    arma::ivec _c_max;
    unsigned int _iter_max;
    double _quantile;
    std::string _stopCriterion;
    double _m;
    arma::vec _w;
    double _alpha;
    double _tol;
    unsigned int _iter4elong;
    double _tol4elong;
    double _max_elong;
    unsigned int _trials_elong;
    double _deltaJK_elong;
    double _max_gap;
    unsigned int _iter4clean;
    double _tol4clean;
    double _quantile4clean;
    long long _seed;
    bool _return_options;
};

#endif // PARAMETERS_HPP