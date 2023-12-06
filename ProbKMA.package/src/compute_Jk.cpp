#include "compute_jk.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

//[[Rcpp::export(.compute_Jk_rcpp)]]
double compute_Jk_rcpp(const arma::field<arma::mat> & v,
                       const arma::ivec & s_k,
                       const arma::vec & p_k,
                       const arma::field<arma::mat> & Y,
                       double alpha,
                       const arma::vec & w,
                       int m,
                       bool use0,
                       bool use1, 
                       double c_k, // in reality is an int
                       arma::vec keep_k) // in reality is an uvec
{
  // domain of the centroid
  const auto & v_dom = domain_rcpp(v(0,0),v(0,1),use0);
  
  // length of the domain
  const unsigned int v_len = v_dom.size();
  
  // select the part of the domain of the centroid
  const auto & v_new = select_domain_rcpp(v(0,0),v(0,1), v_dom, use0, use1);
  
  // dimensionality of the curves
  const unsigned int d = use0 ? Y(0,0).n_cols : Y(0,1).n_cols;
  
  const unsigned int Y_size = Y.n_rows;
  
  arma::field<arma::mat> Y_inters_k(Y_size,2);
  
  arma::uword index_row;
  
  arma::uvec indeces_dom = arma::find(v_dom==0);
  
  for (unsigned int i = 0; i < Y_size; ++i){
    int s_k_i = s_k[i];
    arma::ivec index = arma::regspace<arma::ivec>(1, v_len - std::max(0, 1-s_k_i))+std::max(1,s_k_i)-1;
    unsigned int index_size = index.size();
    if (use0){
      arma::mat new_y0(index_size + std::max(0, 1-s_k_i), d);
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
    }
    if (use1){
      arma::mat new_y1(index_size + std::max(0, 1-s_k_i), d);
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
  
  if(std::isfinite(c_k) && arma::is_finite(keep_k))
  { 
    arma::uvec keep_k_= arma::conv_to<arma::uvec>::from(keep_k);
    arma::ivec supp_inters_length(arma::accu(keep_k_));

    unsigned int k = 0;
    for (arma::uword i = 0; i < Y_size; ++i){
      if (keep_k_(i)){
        supp_inters_length(k++)= arma::accu(domain_rcpp(Y_inters_k(i,0),Y_inters_k(i,1),use0));
      }
    }

    arma::uvec check_lengths = supp_inters_length < static_cast<int>(c_k);
    if (arma::accu(check_lengths) > 0){
      return NA_REAL;
    }
  }

  arma::vec dist(Y_size);
  for (arma::uword i = 0; i < Y_size; ++i){
    dist(i) = diss_d0_d1_L2_rcpp(Y_inters_k(i,0),Y_inters_k(i,1), v_new(0,0),v_new(0,1), w, alpha);
  }
  arma::vec result = dist % pow(p_k,m);
  return arma::accu(result.elem(arma::find_finite(result)));
}
