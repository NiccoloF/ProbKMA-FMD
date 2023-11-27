#include "Motif_L2.hpp"

 arma::field<arma::mat> Motif_L2::compute(const arma::uvec& v_dom,
                                          const arma::vec& s_k,
                                          const arma::vec& p_k,
                                          const arma::field<arma::mat>& Y,
                                          double m) const
 {
   return arma::field<arma::mat>();
 }

  double Motif_L2::test() const
  {
    return x;
  }


RCPP_MODULE(MotifL2Module) {
  Rcpp::class_<Motif_L2>("Motif_L2")
  .constructor()
  .constructor<double>()
  .method("compute", &Motif_L2::compute)
  .field("x",&Motif_L2::x)
  .method("test",&Motif_L2::test);
}