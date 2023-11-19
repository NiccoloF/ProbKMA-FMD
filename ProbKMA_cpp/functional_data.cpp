#include "RcppArmadillo.h"
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


// Abstract class for functional data
class FuncData
{
public:
  
  FuncData() = default; 
  //@TODO: avoid conversion from uword to FuncData
  explicit FuncData(arma::uword N_): Y0(std::vector<arma::mat>(N_)), 
                            Y1(std::vector<arma::mat>(N_)), 
                            lengths(arma::uvec(N_,arma::fill::zeros)),
                            d(0),
                            N(N_) {} 
    
  
  // method used for computing the domain of each curves 
  virtual std::vector<arma::uvec> find_domain_curves() = 0 const; // to be declared as const method
  
  // method used for finding for all the curves the intervals contained in the support
  // it returns a vector such that its element are a list of vector of 3 elements: starting point, end point, length
  virtual std::vector<std::list<arma::uvec>> find_intervals_domains() = 0; // to be declared as const method
  
protected:
  // @TODO: make this class template for the type of datum in mat, i.e. using arma::Mat<T>
  std::vector<arma::mat> Y0;
  // Vector of the derivates of the curves
  std::vector<arma::mat> Y1;
  // Vector of the lengths of the curves
  arma::uvec lengths;
  // Dimensionality of the curves
  arma::uword d;
  // Number of curves
  arma::uword N; 
};


class FuncData_L2d0: public FuncData //only curves 
{
public:
  
  FuncData_L2d0() = default; 
  FuncData_L2d0(const Rcpp::List& Y0_,
                arma::uword N_): FuncData(N_){
   
    for(arma::uword i=0; i < N; ++i){
      Y0[i] = as<arma::mat>(Y0_[i]);
      lengths[i] = Y0[i].n_rows;
    }
    d = Y0[0].n_cols;
  }
  
  virtual std::vector<arma::uvec> find_domain_curves() const 
  {
    std::vector<arma::uvec> domains(N);
    auto zip_views = std::views::zip(domains,Y0,lengths);
    for(auto& tuple : zip_views)
    {
      auto& domains_i = std::get<0>(tuple);
      const auto& Y0_i = std::get<1>(tuple);
      const auto& length_i = std::get<2>(tuple);
      domains_i = arma::uvec{length_i, arma::fill::zeros};
      std::views::iota(0,length_i) 
      | std::views::filter([&Y0_i](auto j){return arma::is_finite(Y0_i.row(j));})
      | std::views::transform([&domains_i](auto j){domains_i[j] = 1;});
    }
    return domains;
  }
  
  std::vector<std::list<arma::uvec>> find_intervals_domains()
  {
    std::vector<std::list<arma::uvec>> Y_intervals(N);
    auto domains = find_domain_curves();
    arma::uword start = 0;
    arma::uword length = 0;
    for (arma::uword i=0; i < N; ++i){
       std::list<arma::uvec> y_intervals;
       arma::uvec y_supp_i(3,arma::fill::zeros);
       arma::uvec indeces_domain = arma::find(domains[i]==1);
       if (indeces_domain(0) == 0 && indeces_domain(1) != 1) 
       { 
         y_intervals.emplace_back(0,1,1);
       }
       for(arma::uword j = 1; j < indeces_domain.size(); ++j){
         if (indeces_domain(j) == indeces_domain(j-1) + 1)
           length++;
         else {
           start = j;
           length = 0;
         }
       }
      y_intervals.emplace_back(start, start + length - 1, length);
    }
    return Y_intervals;  
  }
};


class FuncData_L2d1: public FuncData //only derivatives 
{
public:
  
  FuncData_L2d1() = default; 
  explicit FuncData_L2d1(const Rcpp::List & Y1_,
                arma::uword N_): FuncData(N_) {
    
    for(unsigned int i=0; i < N; ++i){
      Y1[i] = as<arma::mat>(Y1_[i]);
      lengths[i] = Y1[i].n_rows;
    }
    d = Y1[0].n_cols;
  }
  
  virtual std::vector<arma::uvec> find_domain_curves()
  {
    std::vector<arma::uvec> domains(N);
    auto zip_views = std::views::zip(domains,Y1,lengths);
    for(auto& tuple : zip_views)
    {
      auto&& domains_i = std::get<0>(tuple);
      const auto& Y0_i = std::get<1>(tuple);
      const auto& length_i = std::get<2>(tuple);
      domains_i = arma::uvec{length_i, arma::fill::zeros};
      std::views::iota(0,length_i) 
      | std::views::filter([&Y1_i](auto j){return arma::is_finite(Y1_i.row(j));})
      | std::views::transform([&domains_i](auto j){domains_i[j] = 1;});
    }
    return domains;
};


class FuncDataH1: public FuncData //both curves and derivatives also in the zipped version
{
public:
  
  FuncDataH1() = default; 
  FuncDataH1(const Rcpp::List & Y0_,
              const Rcpp::List & Y1_,
              arma::uword N_): FuncData(N_);
 {
    for(auto& tuple: std::views::zip(Y0_,Y0,Y1_,Y1,lengths))
    {
      std::get<1>(tuple) = std::get<0>(tuple);
      std::get<3>(tuple) = std::get<2>(tuple);
      std::get<4>(tuple) = std::get<1>(tuple).n_rows;
    }
    
    d = Y0_[0].n_cols;
  }
  
  virtual std::vector<arma::uvec> find_domain_curves()
  {
    std::vector<arma::uvec> domains(N);
    auto zip_views = std::views::zip(domains,Y0,lengths);
    for(auto& tuple : zip_views)
    {
      auto&& domains_i = std::get<0>(tuple);
      const auto& Y0_i = std::get<1>(tuple);
      const auto& length_i = std::get<2>(tuple);
      domains_i = arma::uvec{length_i, arma::fill::zeros};
      std::views::iota(0,length_i) 
      | std::views::filter([&Y0_i](auto j){return arma::is_finite(Y0_i.row(j));})
      | std::views::transform([&domains_i](auto j){domains_i[j] = 1;});
    }
    return domains;
  }

  // returns a zip views of the objects 
  template <typename A, typename B>
  decltype(auto) zip(const std::vector<A> &a,
                     const std::vector<B> &b)
  { 
    return std::views::zip(a, b);
  }

};


// [[Rcpp::export]]
void test_functional_class(const Rcpp::List & Y0,
                           const Rcpp::List & Y1,
                           arma::uword N){
  FuncDataH1 test
  return;
}

