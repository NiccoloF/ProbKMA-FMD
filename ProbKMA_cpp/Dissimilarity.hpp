#ifndef DISSIMILARITY_HPP
#define DISSIMILARITY_HPP
#include "RcppArmadillo.h"

// Abstract class for dissimilarities
class Dissimilarity
{
public:
  
  Dissimilarity() = default;
  
  // compute dissimilarity 
  // C'è da controllare che una subview di field possa essere convertita dirretamente in field
  virtual double compute(const arma::field<arma::mat>& Y_i,
                         const arma::field<arma::mat>& V_i) const = 0; // to be declared as const

private:

  double distance(const arma::mat& y,
                  const arma::mat& v) const = 0;
};


#endif // DISSIMILARITY_HPP