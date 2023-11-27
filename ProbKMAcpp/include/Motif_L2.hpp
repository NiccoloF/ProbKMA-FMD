#ifndef MOTIF_L2_HPP
#define MOTIF_L2_HPP
#include "Motif.hpp"

class Motif_L2 : public MotifPure
{
public:
  
  Motif_L2() = default;
  Motif_L2(double x_):x(x_){};
  
  virtual arma::field<arma::mat> compute_motif(const arma::uvec& v_dom,
                                               const arma::vec& s_k,
                                               const arma::vec& p_k,
                                               const arma::field<arma::mat>& Y,
                                               double m) const override;
  virtual double test() const override;
  
  virtual ~Motif_L2() = default;
  double x;
  
};

#endif // MOTIF_L2_HPP