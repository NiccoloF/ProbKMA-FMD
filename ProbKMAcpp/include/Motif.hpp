#ifndef MOTIF_HPP
#define MOTIF_HPP
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <numeric>
#include <ranges>
#include <algorithm>
#include <variant>
#include "Utilities.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class MotifPure
{
  public:
    using indexField = std::pair<arma::field<arma::mat>,arma::sword>;
    
    MotifPure() = default;
    
    virtual std::variant<indexField,arma::field<arma::mat>>
      compute_motif(const arma::uvec& v_dom,
                    const arma::ivec& s_k,
                    const arma::vec& p_k,
                    const arma::field<arma::mat>& Y,
                    double m) const = 0;
                          
    virtual ~MotifPure() = default;
    
  protected:
    
    arma::mat compute_v_new(const arma::field<arma::mat> & Y_inters_k,
                            const arma::umat & Y_inters_supp,
                            const arma::uvec & v_dom,
                            arma::uword v_len,
                            const arma::vec & p_k,
                            arma::uword d,
                            arma::uword m) const;
  
};

#endif // MOTIF_HPP

