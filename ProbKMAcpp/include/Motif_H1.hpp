#ifndef MOTIF_H1_HPP
#define MOTIF_H1_HPP
#include "Motif.hpp"

class Motif_H1 : public MotifPure
{
public:
  
  Motif_H1() = default;
  
  virtual ~Motif_H1() = default;
  
  virtual std::variant<indexField,arma::field<arma::mat>>
    compute_motif(const arma::uvec& v_dom,
                  const arma::ivec& s_k,
                  const arma::vec& p_k,
                  const arma::field<arma::mat>& Y,
                  double m) const override;
  
};

#endif // MOTIF_H1_HPP