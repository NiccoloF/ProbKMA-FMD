#include "Dissimilarity_L2.hpp"

L2::L2(const arma::vec& w): Dissimilarity(),_w(w){};

// sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y)
double L2::distance(const arma::mat& y,
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



double L2::compute(const arma::field<arma::mat>& Y_i,
                   const arma::field<arma::mat>& V_i) const
{
    return this->distance(Y_i(0,0),V_i(0,0));
}

RCPP_MODULE(L2Module) {
  Rcpp::class_<L2>("L2")
  .constructor<arma::vec>()
  .method("compute", &L2::compute)
  .field("w",&L2::_w);
}

