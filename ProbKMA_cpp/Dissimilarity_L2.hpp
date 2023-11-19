#ifndef DISSIMILARITY_L2_HPP
#define DISSIMILARITY_L2_HPP

#include "Dissimilarity.hpp"

class L2: public Dissimilarity
{
    L2();

    double distance(const arma::mat& v,
                    const arma::mat& y) const;
    
    virtual double compute(const arma::field<arma::mat>& Y_i,
                           const arma::field<arma::mat>& V_i) const override;
}



#endif // DISSIMILARITY_L2_HPP