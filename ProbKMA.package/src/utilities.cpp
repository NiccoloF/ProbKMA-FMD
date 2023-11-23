#include "utilities.h"

std::vector<arma::ivec> util::repeat_elements(const arma::imat& A,const arma::ivec & times) {
  arma::uword times_size = times.n_elem;
  std::vector<arma::ivec> result((times_size*(times_size+1))/2 - 1);
  std::size_t i = 0;
  for(arma::uword j = 0;j < times_size;++j)
  {
    const arma::imat& B = repmat(A.col(j),1,times[j]);
    B.each_col([&result,&i](const arma::ivec& v){result[i++] = v;});
  }
  return result;
}

Rcpp::IntegerVector util::myseq(int first, int last) {
  Rcpp::IntegerVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

arma::field<arma::mat> util::conv_to_field(const Rcpp::List & Y){
  arma::field<arma::mat> result(Y.size(),2);
  for(unsigned int i=0; i < Y.size(); ++i){
    const Rcpp::List & y = Y[i];
    result(i,0) = Rcpp::as<arma::mat>(y[0]);
    result(i,1) = Rcpp::as<arma::mat>(y[1]);
  }
  return result;
}
