#ifndef DISSIMILARITY_HPP
#define DISSIMILARITY_HPP
#include "TypeTraits.hpp"
#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <Utilities.hpp>
#include <ranges>
#include "Parameters.hpp"

// Abstract class for dissimilarities
class Dissimilarity
{
public:
  
  Dissimilarity() = default;
  
  virtual ~Dissimilarity() = default;
  
  // compute dissimilarity 
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const = 0; 
  
  // Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w, 
                                double alpha, unsigned int c_k) const = 0;

  virtual void set_parameters(const Parameters & newParameters) = 0;
  
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
  
    template<bool use1>
    KMA::vector find_diss_helper(const KMA::Mfield Y,
                                 const KMA::Mfield V,
                                 const KMA::vector& w, 
                                 double alpha, unsigned int c_k) const;
  
    KMA::vector _w;
    
};

#include "Dissimilarity.ipp"

class L2 final: public SobolDiss
{
public:
  
  L2(const KMA::vector& w);
  virtual ~L2() = default;
  
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;
  
  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w, 
                                double alpha, unsigned int c_k) const override;

   void set_parameters(const Parameters & newParameters) override;
  
};

class H1 final: public SobolDiss
{
public:
  
  H1(const KMA::vector& w,double alpha);
  virtual ~H1() = default;
  
  virtual double computeDissimilarity(const KMA::Mfield& Y_i,
                                      const KMA::Mfield& V_i) const override;
  
  virtual KMA::vector find_diss(const KMA::Mfield Y,
                                const KMA::Mfield V,
                                const KMA::vector& w, 
                                double alpha, unsigned int c_k) const override;
  
  void set_parameters(const Parameters & newParameters) override;
  
  double _alpha;

};




#endif // DISSIMILARITY_HPP