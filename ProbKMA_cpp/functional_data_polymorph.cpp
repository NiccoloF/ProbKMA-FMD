#include "RcppArmadillo.h"
using namespace Rcpp;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>
#include <array>
#include <memory>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

class Abstract_Curve // curva astratta
{
public:
  
  Abstract_Curve() = default; 
  
  Abstract_Curve(arma::uword len_,
                 arma::uword d_): len(len_), d(d_) {}
  
  // find the domain of the curve that is an uvec such that = 1 iff the row(i) belongs to the support of curve
  virtual arma::uvec domain() const = 0;
  
  // select the portion of the curve that belongs to the support
  virtual void select_domain() = 0;
  
  // @TODO : compute the distance between two curves
  virtual double distance(const std::shared_ptr<Abstract_Curve>& v, const arma::vec& w, double alpha = 0) const = 0;
  
  // find all the intervals in the support of the curve
  std::vector<std::array<arma::uword,3>>  find_intervals() const{
    std::vector<std::array<arma::uword,3>> y_intervals;
    arma::uword start  = 1;
    arma::uword length = 0;
    arma::uvec indeces_domain = find(domain() ==1); 
    if (indeces_domain(0) == 0 && indeces_domain(1) != 1) // in the case only the first observation is not NA then we have NA
    { 
      y_intervals.emplace_back(std::array<arma::uword,3>{0,1,1});
      start = 2;
    }
    for(arma::uword j = start; j < indeces_domain.size(); ++j){
      if (indeces_domain(j) == indeces_domain(j-1) + 1){
        length++;
        if (j == indeces_domain.size() - 1)  // for the last interval if we are ending with not a NA
          y_intervals.emplace_back(std::array<arma::uword,3>{1 + indeces_domain(j) - length, 1 + indeces_domain(j), length + 1});
      } 
      else 
      {
        y_intervals.emplace_back(std::array<arma::uword,3>{1 + indeces_domain(j-1) - length, 1 + indeces_domain(j-1), 1 + length});
        length = 0;
      }
    }
    return y_intervals;
  }
  
  virtual const arma::mat & get_curve() const  = 0;
  
  virtual const arma::mat & get_derivative() const  = 0;
  
  virtual std::shared_ptr<Abstract_Curve> shift_curve(arma::sword s, arma::uword v_len) const = 0;
  
protected:
  
  arma::uword len;
  arma::uword d;
  
};

class Curve : public Abstract_Curve
{
public:
  Curve() = default;
  
  Curve(const arma::mat& y): Abstract_Curve(y.n_rows,y.n_cols), y(y) {}
  
  // domain based on y0
  arma::uvec domain() const override{
    arma::uvec result(len,arma::fill::zeros);
    for (arma::uword i=0; i < len; ++i){
      const arma::uvec & finite_row = find_finite(y.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return result;
  }
  
  void select_domain() override{
    arma::uvec indeces = find(domain() == 1);
    y = y.rows(indeces); 
  }
  
  const arma::mat & get_curve() const override {
    return y;
  }
  
  const arma::mat & get_derivative() const override { // if Curve is used for storing only the derivative
    return y;
  }
  
  double distance(const std::shared_ptr<Abstract_Curve>& v, //v is a Curve in reality 
                  const arma::vec& w,
                  double alpha = 0) const override {
    arma::mat diff = arma::square(y - v->get_curve()); //(y-v)^2
                    
    diff.replace(arma::datum::nan,0);
                    
    const arma::rowvec & col_sum = arma::sum(diff,0); //colSums
                    
    unsigned int num_rows = 0;
    for(arma::uword i = 0; i < len; ++i)
      if(is_finite(y.row(i))){num_rows += 1;}
                      
    arma::urowvec length_dom(d,arma::fill::value(num_rows)); //length of the domain
                      
    return sum((col_sum/length_dom)%w.t())/d; 
  }
  
  std::shared_ptr<Abstract_Curve> shift_curve(arma::sword s, arma::uword v_len) const override {
    arma::ivec index = std::max(1,s) - 1 + arma::regspace<arma::ivec>(1, v_len - std::max(0,1-s));
    const int index_size = index.size();
    arma::mat shift_y(std::max(0,1-s) + index_size,d);
    shift_y.fill(arma::datum::nan);
    auto filtered_j = std::views::iota(0,index_size)
      | std::views::filter([this,&index](int j){return (index[j] <= this->len);});
    for(int j : filtered_j) 
      shift_y.row(std::max(0,1-s) + j) =  y.row(index(j) - 1); 
    return std::make_shared<Curve>(Curve(shift_y));
    
  }
  
protected:
  arma::mat y;
  
};


class Curve_d0d1 : public Abstract_Curve
{
public:
  Curve_d0d1(const arma::mat& y0,
             const arma::mat& y1): Abstract_Curve(y0.n_rows,y0.n_cols) 
  {
    y = std::make_pair<Curve,Curve>(Curve(y0),Curve(y1));
  }
  
  arma::uvec domain() const override {
    return y.first.domain();
  }
  
  void select_domain() override {
    y.first.select_domain();
    y.second.select_domain();
  }
  
  double distance(const std::shared_ptr<Abstract_Curve>& v, //v is a Curve_d0_d1
                  const arma::vec& w,
                  double alpha = 0) const override {
    return (1-alpha)*y.first.distance(v,w) + alpha*y.second.distance(std::make_shared<Curve>(v->get_derivative()),w);
  }
  
  const arma::mat & get_curve() const override {
    return y.first.get_curve();
  }
  
  const arma::mat & get_derivative() const override {
    return y.second.get_curve();
  }
  
  std::shared_ptr<Abstract_Curve> shift_curve(arma::sword s, arma::uword v_len) const override {
    std::shared_ptr<Abstract_Curve> y0_shift = y.first.shift_curve(s,v_len);
    std::shared_ptr<Abstract_Curve> y1_shift = y.second.shift_curve(s,v_len);
    return std::make_shared<Curve_d0d1>(Curve_d0d1(y0_shift->get_curve(),y1_shift->get_curve()));
  }
  
protected:
  
  std::pair<Curve,Curve> y;
  
};

// Abstract class for functional data
class FuncData
{
public:
  
  FuncData() = default; 
  //@TODO: avoid conversion from uword to FuncDataset
  FuncData(arma::uword N_): Y(std::vector<std::shared_ptr<Abstract_Curve>>(N_)), N(N_) {}
  
  // this constructor is the one that is called from R
  FuncData(const Rcpp::List & Y0, //I am assuming that Y0 is always passed for initialize the correct N
           const Rcpp::List & Y1,
           bool use0, bool use1): FuncData(Y0.size()) {
    if (use0 && !use1) {
      for(arma::uword i=0; i < N; ++i){
        Y[i] = std::make_shared<Curve>(Rcpp::as<arma::mat>(Y0[i]));
      }
    } else if (use1 && !use0){
      for(arma::uword i=0; i < N; ++i){
        Y[i] = std::make_shared<Curve>(Rcpp::as<arma::mat>(Y1[i]));
      }
    } else {
      for(arma::uword i=0; i < N; ++i){
        Y[i] = std::make_shared<Curve_d0d1>(Rcpp::as<arma::mat>(Y0[i]),
                                            Rcpp::as<arma::mat>(Y1[i]));
      }
    }
  }
    
  // method used for computing the domain of each curves 
  std::vector<arma::uvec> find_domains_curves() const {
    std::vector<arma::uvec> domains(N);
    for(arma::uword i=0; i < N; ++i)
      domains[i] = Y[i]->domain();
    return domains;
  }
  
  // method used for finding for all the curves the intervals contained in the support
  // it returns a vector whose elements are a list of vector of 3 elements: starting point, end point, length
  std::vector<std::vector<std::array<arma::uword,3>>> find_intervals_domains() const {
    std::vector<std::vector<std::array<arma::uword,3>>> Y_intervals(N);
    for(arma::uword i=0; i < N; ++i)
      Y_intervals[i] = Y[i]->find_intervals();
    return Y_intervals;
  }
  
  // shift curves shitfs all the curves in the dataset according to their specific shift 
  FuncData shift_curves(const arma::ivec & S_k, const arma::uword v_len){
    FuncData result(N);
    auto & shifted_curves = result.Y; // we should add a getter
    for(arma::uword i=0; i < N; ++i)
      shifted_curves[i] = Y[i]->shift_curve(S_k[i], v_len);
    return result;
  }
  
protected:
  // all the curves are stored in Y as a shared_ptr in order to use polymorph
  std::vector<std::shared_ptr<Abstract_Curve>> Y;
  // number of curves
  arma::uword N; 
  
};

// [[Rcpp::export]]
Rcpp::List test_functional_class(const Rcpp::List & Y0,
                                 const Rcpp::List & Y1,
                                 const Rcpp::List & V0,
                                 const Rcpp::List & V1,
                                 const arma::vec w,
                                 double alpha){
  std::shared_ptr<Abstract_Curve> datum = std::make_shared<Curve_d0d1>(Y0[0],Y1[0]);
  std::shared_ptr<Abstract_Curve> v = std::make_shared<Curve_d0d1>(V0[0],V1[0]);
  FuncData dataset(Y0,Y1,true,true);
  return List::create(dataset.find_intervals_domains(), datum->distance(v,w,alpha), datum->shift_curve(29,40)->get_curve(),datum->shift_curve(29,40)->get_derivative());
}

