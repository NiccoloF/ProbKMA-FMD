#include "Dissimilarity_H1.hpp"

H1::H1(const arma::vec& w,double alpha):Dissimilarity(),_w(w),_alpha(alpha){};

double H1::distance(const arma::mat& y,
                    const arma::mat& v) const
{
    arma::mat diff = arma::square(y - v); //(y-v)^2
    
    diff.replace(arma::datum::nan,0);
    
    const arma::rowvec & col_sum = arma::sum(diff,0); //colSums
    
    unsigned int n_rows = 0;
    for(arma::uword i = 0; i < y.n_rows; ++i)
        if(is_finite(y.row(i))){n_rows += 1;}
        
    arma::urowvec length_dom(y.n_cols,arma::fill::value(n_rows)); //length of the domain
    
    return sum((col_sum/length_dom)%_w.t())/y.n_cols;
};

double H1::computeDissimilarity(const arma::field<arma::mat>& Y_i,
                                const arma::field<arma::mat>& V_i) const
{
    return (1-_alpha) * this -> distance(Y_i(0,0),V_i(0,0)) + 
            _alpha * this -> distance(Y_i(0,1),V_i(0,1));
}


RCPP_MODULE(H1Module) {
  Rcpp::class_<H1>("H1")
  .constructor<arma::vec,double>()
  .field("w",&H1::_w)
  .field("alpha",&H1::_alpha);
}

