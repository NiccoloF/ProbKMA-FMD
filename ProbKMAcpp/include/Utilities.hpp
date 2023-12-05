#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__
#include "RcppArmadillo.h"
#include "TypeTraits.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

namespace util
{

  template <typename MatType, bool Use>
  auto selectDomainHelper(const MatType& mat_v, const KMA::uvector& dom) {
    if constexpr (Use) {
      return mat_v.rows(dom);
    } else {
      return MatType();  // Return an empty matrix if not used
    }
  }
  
  // restituire direttamente Mfield
  template <typename... MatTypes, bool... Uses>
  KMA::Mfield selectDomain(const KMA::uvector& dom, const MatTypes&... mat_v) {
    static_assert(sizeof...(MatTypes) == sizeof...(Uses), "Error: Mismatch in the number of matrices and use flags.");
    
    KMA::Mfield result(sizeof...(MatTypes));
    
    auto processMatrix = [&dom](const auto& mat_v, bool use) {
      return selectDomainHelper<decltype(mat_v), use>(mat_v, dom);
    };
    
    size_t index = 0;
    ((result(index++) = processMatrix(mat_v, Uses)), ...);
    
    return result;
  }
  
  // returns a rowvector 
  template <typename MatType>
  arma::urowvec findDomain(const MatType& v) {
    arma::urowvec result(v.n_rows, arma::fill::zeros);
    for (arma::uword i = 0; i < v.n_rows; ++i) {
      const KMA::uvector& finite_row = arma::find_finite(v.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return result;
  }

}
  


#endif // __UTILITIES_HPP__