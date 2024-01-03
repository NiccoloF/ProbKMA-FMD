#ifndef __PERFORMANCE_INDEX__
#define __PERFORMANCE_INDEX__
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"
#include "Utilities.hpp"
#include "Dissimilarity.hpp"
#include <memory>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class PerformanceIndexAB
{
public:
  
  virtual double compute_Jk(const KMA::Mfield& v,
                            const KMA::ivector& s_k,
                            const KMA::vector& p_k,
                            const KMA::Mfield& Y,
                            const KMA::vector& w,
                            int m,
                            double c_k, // actually is an int
                            const KMA::uvector & keep_k, // actually is an KMA::uvector
                            const std::shared_ptr<Dissimilarity>& diss) const = 0; 

  protected:
  
  template<bool use1>
  void shiftCurveHandle(KMA::Mfield& Y_inters_k,
                        const KMA::Mfield& Y,
                        const KMA::ivector& s_k,
                        const arma::urowvec& v_dom) const
  {
    std::size_t Y_size = Y.n_rows;
    std::size_t v_len = v_dom.size();
    KMA::uvector indeces_dom = arma::find(v_dom==0);
    arma::uword index_row;
    
    for (unsigned int i = 0; i < Y_size; ++i){
      int s_k_i = s_k[i];
      arma::ivec index = arma::regspace<arma::ivec>(1, v_len - std::max(0, 1-s_k_i))+std::max(1,s_k_i)-1;
      unsigned int index_size = index.size();
      
      arma::mat new_y0(index_size + std::max(0, 1-s_k_i), Y(0,0).n_cols);
      const int y_len = Y(i,0).n_rows;
      new_y0.fill(arma::datum::nan);
      for(unsigned int j = 0; j < index_size; ++j) {
        if (index[j]  <= y_len){
          index_row = std::max(0, 1-s_k_i) + j;
          new_y0.row(index_row) =  Y(i,0).row(index[j] - 1);
        }
      }
      new_y0.shed_rows(indeces_dom);
      Y_inters_k(i,0) = new_y0;
      
      if constexpr(use1){
        arma::mat new_y1(index_size + std::max(0, 1-s_k_i), Y(0,1).n_cols);
        const int y_len = Y(i,1).n_rows;
        new_y1.fill(arma::datum::nan);
        for(unsigned int j = 0; j < index_size; ++j) {
          if (index[j] <= y_len){
            index_row = std::max(0, 1-s_k_i) + j;
            new_y1.row(index_row) =  Y(i,1).row(index[j] - 1);
          }
        }
        new_y1.shed_rows(indeces_dom);
        Y_inters_k(i,1) = new_y1;
      }
    } 
  }

};

class PerformanceSobol: public PerformanceIndexAB
{
protected:
  
  template<bool use1>
  double compute_Jk_helper(const KMA::Mfield& v,const KMA::ivector& s_k,
                           const KMA::vector& p_k,const KMA::Mfield& Y,
                           const KMA::vector& w,int m,double c_k, 
                           const KMA::uvector& keep_k,
                           const std::shared_ptr<Dissimilarity>& diss) const;
};


class PerformanceL2 final: public PerformanceSobol
{
  
  double compute_Jk(const KMA::Mfield& v,const KMA::ivector& s_k,
                    const KMA::vector& p_k,const KMA::Mfield& Y,
                    const KMA::vector& w,int m,double c_k, 
                    const KMA::uvector& keep_k,
                    const std::shared_ptr<Dissimilarity>& diss) const override;
};


class PerformanceH1 final: public PerformanceSobol
{
  
  double compute_Jk(const KMA::Mfield& v,const KMA::ivector& s_k,
                    const KMA::vector& p_k,const KMA::Mfield& Y,
                    const KMA::vector& w,int m,double c_k,
                    const KMA::uvector& keep_k,
                    const std::shared_ptr<Dissimilarity>& diss) const override;
  
};

#include "PerformanceIndex.ipp"

#endif //__PERFORMANCE_INDEX__