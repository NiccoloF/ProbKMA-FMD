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
    ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,const Rcpp::List& parameters,
            const arma::mat& P0,const arma::mat& S0);
    
    virtual ~ProbKMA() = default;
    
    // run probKMA's algorithm
    Rcpp::List probKMA_run(const SEXP& dissimilarity,const SEXP& motif) const;

    void set_parameters(const Rcpp::List& parameters);
    
 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<_probKMAImp> _probKMA;
};

#endif // PROBKMA_HPP




