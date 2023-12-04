#ifndef PROBKMA_HPP
#define PROBKMA_HPP

#include "RcppArmadillo.h"
#include "Parameters.hpp"
#include "Dissimilarity.hpp"
#include "Motif.hpp"
#include <vector>
#include <ranges>
#include <algorithm>
#include <memory>
#include <Rcpp.h>

// Forward declaration
class _probKMAImp;

class ProbKMA
{
 public:
    using matrix = arma::mat;
    using imatrix = arma::imat;
    using umatrix = arma::umat;
    using vector = arma::vec;
    using ivector = arma::ivec;
    using uvector = arma::uvec;
  
    // Y: a list containing two list -> Y0 and Y1
    // V: a list containing two list -> V0 and V1
    ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,
            const Rcpp::List& parameters,
            const matrix& P0,const imatrix& S0,
            const std::string& diss);
   
    virtual ~ProbKMA() = default;
    
    // run probKMA's algorithm
    Rcpp::List probKMA_run() const;

    void set_parameters(const Rcpp::List& parameters);
    
 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<_probKMAImp> _probKMA;
};

  Rcpp::List initialChecks(const Rcpp::List& Y0,const Rcpp::List& Y1,
                           const Rcpp::NumericMatrix& P0,
                           const Rcpp::NumericMatrix& S0,
                           const Rcpp::List& params,
                           const Rcpp::String diss,
                           const double alpha,
                           const Rcpp::NumericVector& w);


#endif // PROBKMA_HPP




