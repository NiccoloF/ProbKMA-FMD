#ifndef MOTIF_HPP
#define MOTIF_HPP
#include "RcppArmadillo.h"
#include <Rcpp.h>

class MotifBase
{
  public:
    
    MotifBase() = default;
    
    virtual arma::field<arma::mat> compute(const arma::uvec& v_dom,
                                          const arma::vec& s_k,
                                          const arma::vec& p_k,
                                          const arma::field<arma::mat>& Y,
                                          double m) const = 0;
  
  virtual ~MotifBase() = default;
};

#endif // MOTIF_HPP

