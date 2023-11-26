#include "Motif_H1.hpp"

arma::field<arma::mat> Motif_H1::compute(const arma::uvec& v_dom,
                                         const arma::vec& s_k,
                                         const arma::vec& p_k,
                                         const arma::field<arma::mat>& Y,
                                         double m) const
{
  return arma::field<arma::mat>();
}


RCPP_MODULE(MotifH1Module) {
  Rcpp::class_<Motif_H1>("Motif_H1")
  .constructor()
  .method("compute", &Motif_H1::compute);
}