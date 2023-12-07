#ifndef DISSIMILARITY_HPP
#define DISSIMILARITY_HPP
#include "TypeTraits.hpp"
#include "RcppArmadillo.h"
#include <Rcpp.h>

// Abstract class for dissimilarities
class Dissimilarity
{
public:
  
  Dissimilarity() = default;
  
  virtual ~Dissimilarity() = default;
  
  // compute dissimilarity 
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const = 0; 
  
protected:

  virtual double distance(const KMA::matrix& y,
                          const KMA::matrix& v) const = 0;
};

class SobolDiss : public Dissimilarity
{
public:
  
    SobolDiss(const KMA::vector& w);
    
protected:
  
    virtual double distance(const KMA::matrix& y,
                            const KMA::matrix& v) const override;
  
    KMA::vector _w;
    
};


class L2 final: public SobolDiss
{
public:
  
  L2(const KMA::vector& w);
  virtual ~L2() = default;
  
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;
  
};

class H1 final: public SobolDiss
{
public:
  
  H1(const KMA::vector& w,double alpha);
  virtual ~H1() = default;
  
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;
  
  double _alpha;

};




#endif // DISSIMILARITY_HPP