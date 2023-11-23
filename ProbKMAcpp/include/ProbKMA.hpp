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

// Forward declaration
class _probKMAImp;

class ProbKMA
{
 public:
    ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
            const arma::mat& S0,SEXP parameters,
            SEXP dissimilarity,SEXP motif);
    
    ~ProbKMA() = default;
    
    // run probKMA's algorithm
    Rcpp::List probKMA_run() const;

    void set_parameters(SEXP parameters);
    void set_distance(SEXP pnewDistance);
    void set_motif(SEXP pnewMotif);
    
 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<_probKMAImp> _probKMA;
};

#endif // PROBKMA_HPP




