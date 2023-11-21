#ifndef DISSIMILARITY_H1_HPP
#define DISSIMILARITY_H1_HPP

#include "Dissimilarity.hpp"

class H1: public Dissimilarity
{
    public:

    H1();

    virtual double compute(const arma::field<arma::mat>& Y_i,
                           const arma::field<arma::mat>& V_i) const override;

    private:

    virtual double distance(const arma::mat& y,
                            const arma::mat& v) const override ;

    protected:
    
    alpha;
}



#endif // DISSIMILARITY_H1_HPP