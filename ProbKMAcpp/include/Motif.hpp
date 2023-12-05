#ifndef MOTIF_HPP
#define MOTIF_HPP
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include <Rcpp.h>
#include <numeric>
#include <ranges>
#include <algorithm>
#include <variant>
#include "Utilities.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class MotifPure
{
  public:
    using indexField = std::pair<KMA::Mfield,arma::sword>;
    
    MotifPure() = default;
  
    virtual std::variant<indexField,KMA::Mfield>
      compute_motif(const arma::urowvec& v_dom,
                    const KMA::ivector& s_k,
                    const KMA::vector& p_k,
                    const KMA::Mfield& Y,
                    double m) const = 0;
                          
    virtual ~MotifPure() = default;
    
  protected:
    
    KMA::matrix compute_v_new(const KMA::Mfield& Y_inters_k,
                            const KMA::umatrix& Y_inters_supp,
                            const arma::urowvec & v_dom,
                            arma::uword v_len,
                            const KMA::vector & p_k,
                            arma::uword d,
                            arma::uword m) const;
  
};


class Motif_L2 final: public MotifPure
{
public:
  
  Motif_L2() = default;
  
  virtual std::variant<indexField,KMA::Mfield>
    compute_motif(const arma::urowvec& v_dom,
                  const KMA::ivector& s_k,
                  const KMA::vector& p_k,
                  const KMA::Mfield& Y,
                  double m) const override;
  
  virtual ~Motif_L2() = default;
  
};

class Motif_H1 final: public MotifPure
{
public:
  
  Motif_H1() = default;
  
  virtual ~Motif_H1() = default;
  
  virtual std::variant<indexField,KMA::Mfield>
    compute_motif(const arma::urowvec& v_dom,
                  const KMA::ivector& s_k,
                  const KMA::vector& p_k,
                  const KMA::Mfield& Y,
                  double m) const override;

};

#endif // MOTIF_HPP

