#include "find_min_diss.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


arma::uvec domain_rcpp(const arma::mat & v0,
                       const arma::mat & v1,
                       bool use0){
  if (use0){
    arma::uvec result(v0.n_rows,arma::fill::zeros);
    for (arma::uword i=0; i < v0.n_rows; ++i){
      const arma::uvec & finite_row = find_finite(v0.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return result;
  } else {
    arma::uvec result(v1.n_rows,arma::fill::zeros);
    for (arma::uword i=0; i < v1.n_rows; ++i){
      const arma::uvec & finite_row = find_finite(v1.row(i));
      if(finite_row.n_elem) 
        result(i) = 1;
    }
    return result;
  }
}

arma::urowvec domain_rcpp_base(const arma::mat & v)
{
    arma::urowvec result(v.n_rows,arma::fill::zeros);
    for (arma::uword i=0; i < v.n_rows; ++i){
      const arma::uvec & finite_row = arma::find_finite(v.row(i));
      if(finite_row.n_elem)
        result(i) = 1;
    }
    return result;
}

arma::field<arma::mat> select_domain_rcpp(const arma::mat & v0,
                                          const arma::mat & v1,
                                          const arma::uvec &v_dom,
                                          bool use0,
                                          bool use1){ 
  arma::uvec dom = arma::find(v_dom==1);
  arma::field<arma::mat> v(1,2);
  if(use0) {
    v(0,0) = v0.rows(dom);
  }
  if(use1) {
    v(0,1) = v1.rows(dom);
  }
  return v;
}

double distance(const arma::mat& v,
                const arma::mat& y,
                const arma::vec& w)
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

double diss_d0_d1_L2_rcpp(const arma::mat & y0,
                          const arma::mat & y1,
                          const arma::mat & v0,
                          const arma::mat & v1,
                          const arma::vec& w,
                          double alpha)
{
  if (alpha == 0)
    return distance(y0,v0,w);
  if (alpha == 1)
    return distance(y1,v1,w);

  return (1 - alpha)*distance(y0,v0,w) + alpha*distance(y1,v1,w);
}

// [[Rcpp::export()]]
arma::vec find_diss_rcpp(const arma::mat &y0,
                         const arma::mat &y1,
                         const arma::mat &v0,
                         const arma::mat &v1,
                         const arma::vec & w, 
                         double alpha, unsigned int c_k,
                         unsigned int d,bool use0,bool use1)
{
  // Convert domain and select_domain
  arma::uvec v_dom = domain_rcpp(v0,v1,use0);
  const auto & v_new = select_domain_rcpp(v0,v1,v_dom,use0,use1);
  int v_len = v_dom.size();
  int y_len = y0.n_rows;
  arma::ivec s_rep = arma::regspace<arma::ivec>(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
  const unsigned int s_rep_size = s_rep.size();
  arma::field<arma::mat> y_rep(s_rep_size,2);
  
  const int index_size = v_len;
  auto index_range = std::views::iota(0,index_size);
  arma::uvec indeces_dom = arma::find(v_dom==0);
  for (unsigned int i = 0; i < s_rep_size; ++i) {
    arma::ivec index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
    auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));}); // verificare il v_dom
    if (use0) {
      y_rep(i,0).resize(index_size, d);
      y_rep(i,0).fill(arma::datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&y_rep,&y0,&index,&i](int j){y_rep(i,0).row(j) = y0.row(index[j] - 1);});
      y_rep(i,0).shed_rows(indeces_dom);
    }
    if (use1) {
      y_rep(i,1).resize(index_size, d);
      y_rep(i,1).fill(arma::datum::nan);
      std::for_each(j_true.begin(),j_true.end(),[&y_rep,&y1,&index,&i](int j){y_rep(i,1).row(j) = y1.row(index[j] - 1);});
      y_rep(i,1).shed_rows(indeces_dom);
    }
  }
  
  arma::ivec length_inter(s_rep_size);
  arma::uvec non_na_indices;
  for (unsigned int i = 0; i < s_rep_size; ++i){
    if (use0) {
      non_na_indices = arma::find_nan(y_rep(i,0).col(0));
      length_inter(i) = y_rep(i,0).col(0).n_elem - non_na_indices.n_elem;
    }
    else {
      non_na_indices = arma::find_nan(y_rep(i,1).col(0));
      length_inter(i) = y_rep(i,1).col(0).n_elem - non_na_indices.n_elem;  
    }
  }
  
  arma::uvec valid = length_inter >= c_k;
  if (arma::accu(valid) == 0) {        
    valid.elem(arma::find(length_inter == arma::max(length_inter)))+= 1;  
  }
  
  double min_d = std::numeric_limits<double>::max();
  int min_s = 0;
  
  for (unsigned int i = 0; i < s_rep_size; i++) {
    if (valid(i)) {
      double dist = diss_d0_d1_L2_rcpp(y_rep(i,0), y_rep(i,1),v_new(0,0),v_new(0,1), w, alpha);
      if (dist < min_d){
        min_d = dist;
        min_s = s_rep[i];
      } 
    }
  }

  return arma::vec({static_cast<double>(min_s), min_d});  
}

arma::vec find_diss_aligned_rcpp(const arma::mat &y0,
                                 const arma::mat &y1,
                                 const arma::mat &v0,
                                 const arma::mat &v1, 
                                 const arma::vec & w, 
                                 double alpha,
                                 bool aligned,
                                 unsigned int d,
                                 bool use0,
                                 bool use1)
  {
    arma::uvec v_dom = domain_rcpp(v0,v1,use0);
    auto v_new = select_domain_rcpp(v0,v1,v_dom,use0,use1);
    int v_len = v_dom.size();
    int y_len = y0.n_rows;
    arma::ivec s_rep;
    if (aligned){
      s_rep = 1;
    } else {
      s_rep = arma::regspace<arma::ivec>(1, y_len - v_len + 1);
    }
    const unsigned int s_rep_size = s_rep.size();
    arma::field<arma::mat> y_rep(s_rep_size,2);
    const int index_size = v_len;
    auto index_range = std::views::iota(0,index_size);
    arma::uvec indeces_dom = arma::find(v_dom==0);
    for (unsigned int i = 0; i < s_rep_size; ++i) {
      arma::ivec index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
      auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
      if (use0) {
        y_rep(i,0).resize(index_size, d);
        y_rep(i,0).fill(arma::datum::nan);
        std::for_each(j_true.begin(),j_true.end(),[&y_rep,&y0,&index,&i](int j){y_rep(i,0).row(j) = y0.row(index[j] - 1);});
        y_rep(i,0).shed_rows(indeces_dom);
      }
      if (use1) {
        y_rep(i,1).resize(index_size, d);
        y_rep(i,1).fill(arma::datum::nan);
        std::for_each(j_true.begin(),j_true.end(),[&y_rep,&y1,&index,&i](int j){y_rep(i,1).row(j) = y1.row(index[j] - 1);});
        y_rep(i,1).shed_rows(indeces_dom);
      }
    }
    
    double min_d = std::numeric_limits<double>::max();
    int min_s = 0;
    
    for (unsigned int i = 0; i < s_rep_size; i++) {
      double dist = diss_d0_d1_L2_rcpp(y_rep(i,0),y_rep(i,1), v_new(0,0),v_new(0,1), w, alpha);
      if (dist < min_d){
        min_d = dist;
        min_s = s_rep(i);
      }
    }
    return arma::vec({static_cast<double>(min_s),min_d});
  }

// [[Rcpp::export(.find_shift_warp_min)]]
Rcpp::List find_shift_warp_min(const Rcpp::List & Y, 
                               const Rcpp::List & V_new,
                               const arma::vec & w,
                               const arma::ivec & c_k,
                               unsigned int K,
                               unsigned int d,
                               double max_gap,
                               double alpha,
                               bool use0,
                               bool use1) {
  
  
  const unsigned int Y_size = Y.size();
  const unsigned int V_new_size = V_new.size();
  arma::vec sd(2);
  arma::imat S_new(Y_size,V_new_size);
  arma::mat  D_new(Y_size,V_new_size);
  arma::field<arma::mat> Y_(util::conv_to_field(Y,use0,use1));
  arma::field<arma::mat> V_new_(util::conv_to_field(V_new,use0,use1));

  #ifdef _OPENMP
      #pragma omp parallel for collapse(2) firstprivate(sd)
  #endif
  for (unsigned int i = 0; i < V_new_size; ++i)
    for (unsigned int j = 0; j < Y_size; ++j){ 
      sd = find_diss_rcpp(Y_(j,0),Y_(j,1),V_new_(i,0),V_new_(i,1),w,alpha,c_k(i),d,use0,use1); 
      S_new(j,i) = sd(0);
      D_new(j,i) = sd(1);
    }
    
  return Rcpp::List::create(S_new,D_new);
}

