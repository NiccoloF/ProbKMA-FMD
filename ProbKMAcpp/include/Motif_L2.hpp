#ifndef MOTIF_L2_HPP
#define MOTIF_L2_HPP
#include "Motif.hpp"

class Motif_L2 : public MotifBase
{
public:
  
  Motif_L2() = default;
  
  virtual arma::field<arma::mat> compute(const arma::uvec& v_dom,
                                         const arma::vec& s_k,
                                         const arma::vec& p_k,
                                         const arma::field<arma::mat>& Y,
                                         double m) const override;
  
  virtual ~Motif_L2() = default;
  
};

#endif // MOTIF_L2_HPP