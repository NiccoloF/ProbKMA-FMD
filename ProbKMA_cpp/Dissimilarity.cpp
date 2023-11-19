#include "RcppArmadillo.h"
using namespace Rcpp;
#include <numeric>
#include <vector>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <typeinfo>



// L2 distance (change name)
class L2: public Dissimilarity
{
public:
  L2():Dissimilarity() {};  
  L2(const arma::vec& w_):Dissimilarity(w_) {};  // constructor taking arma::vec w_ as input
  
  // sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y)
  double distance(const arma::mat& v,
                  const arma::mat& y)
                  {
                    arma::mat diff = arma::square(y - v); //(y-v)^2
                    
                    diff.replace(arma::datum::nan,0);
                    
                    const arma::rowvec & col_sum = arma::sum(diff,0); //colSums
                    
                    unsigned int n_rows = 0;
                    for(arma::uword i = 0; i < y.n_rows; ++i)
                      if(is_finite(y.row(i))){n_rows += 1;}
                      
                    arma::urowvec length_dom(y.n_cols,arma::fill::value(n_rows)); //length of the domain
                    
                    return sum((col_sum/length_dom)%w.t())/y.n_cols;
                  };
      
  virtual double compute(const std::pair<arma::mat,arma::mat> &v,
                         const std::pair<arma::mat,arma::mat> &y) override 
                         {
                          return this->distance(v.first,y.first);
                         }
  
};

class H1 final: public L2 
{
public:
  
  H1():L2(), alpha(0.0) {};  
  H1(const arma::vec& w_, double alpha_):L2(w_), alpha(alpha_) {};  // constructor taking arma::vec w_ and double alpha_ as input
  
  // Override the compute method
  virtual double compute(const std::pair<arma::mat,arma::mat> &v,
                         const std::pair<arma::mat,arma::mat> &y) override
                         {
                          return alpha == 1? distance(v.second,y.second) : 
                                 (1-alpha)*distance(v.first,y.first) + 
                                 alpha*distance(v.second,y.second);
                         }
  
protected:
  
  double alpha;  
};


// [[Rcpp::export]]
double test_dissimilarities_class(const List & v,
                                  const List & y,
                                  const arma::vec &w,
                                  double alpha){
  arma::mat y0 = as<arma::mat>(y[0]);
  arma::mat y1 = as<arma::mat>(y[1]);
  arma::mat v0 = as<arma::mat>(v[0]);
  arma::mat v1 = as<arma::mat>(v[1]);
  std::pair<arma::mat,arma::mat> y_pair(std::make_pair<arma::mat,arma::mat>(std::move(y0),std::move(y1)));
  std::pair<arma::mat,arma::mat> v_pair(std::make_pair<arma::mat,arma::mat>(std::move(v0),std::move(v1)));
  H1 dist(w,alpha);
  return dist.compute(v_pair,y_pair);
}