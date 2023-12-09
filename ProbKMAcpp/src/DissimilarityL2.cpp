#include "Dissimilarity.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

L2::L2(const KMA::vector& w):SobolDiss(w) {};


double L2::computeDissimilarity(const KMA::Mfield& Y_i,
                                const KMA::Mfield& V_i) const
{
    return this->distance(Y_i(0,0),V_i(0,0));
}

KMA::vector L2::find_diss(const KMA::Mfield Y,
                          const KMA::Mfield V,
                          const KMA::vector& w, 
                          double alpha, unsigned int c_k) const
{
  return find_diss_helper<false>(Y,V,w,alpha,c_k);
}

