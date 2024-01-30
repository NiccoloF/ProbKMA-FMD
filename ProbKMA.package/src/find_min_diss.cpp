#include "find_min_diss.h"


// [[Rcpp::export(.find_diss)]]
arma::vec find_diss(const Rcpp::List &y,const Rcpp::List &v,
                    const arma::vec & w,
                    double alpha, unsigned int c_k,
                    unsigned int d,bool use0,bool use1,
                    const Rcpp::Function & domain,
                    const Rcpp::Function & select_domain,
                    const Rcpp::Function & diss_d0_d1_L2)
  {
    // Convert domain and select_domain
    Rcpp::LogicalVector v_dom = Rcpp::as<Rcpp::LogicalVector>(domain(v,use0));
    Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
    int v_len = v_dom.size();
    int y_len = use0? Rcpp::as<arma::mat>(y[0]).n_rows :  Rcpp::as<arma::mat>(y[1]).n_rows;
    arma::ivec s_rep = arma::regspace<arma::ivec>(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
    arma::uword s_rep_size = s_rep.size();
    Rcpp::List y_rep(s_rep_size);

    // Convert y[0] and y[1] to mat objects
    arma::mat temp_y0, temp_y1;
    if (use0) {
      temp_y0 = Rcpp::as<arma::mat>(y[0]);
    }
    if (use1) {
      temp_y1 = Rcpp::as<arma::mat>(y[1]);
    }

    arma::mat new_y0(v_len, d);
    arma::mat new_y1(v_len, d);
    arma::ivec index;
    Rcpp::List y_rep_i;
    arma::uvec neg_index;
    arma::uvec filtered_j;
    for (unsigned int i = 0; i < s_rep_size; ++i) {
      index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
      filtered_j = arma::find((index > 0) && (index <= y_len));
      neg_index = arma::find(index <= 0);
      y_rep_i = Rcpp::List::create(Rcpp::Named("y0") = R_NilValue, Rcpp::Named("y1") = R_NilValue);
      if (use0) {
        new_y0.fill(arma::datum::nan);
        new_y0.rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = temp_y0.rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
        y_rep_i["y0"] = new_y0;
      }
      if (use1) {
        new_y1.fill(arma::datum::nan);
        new_y1.rows(neg_index.n_elem,neg_index.n_elem + filtered_j.n_elem - 1) = temp_y1.rows(index(*(filtered_j.cbegin())) - 1, index(*(filtered_j.cend() - 1)) - 1);
        y_rep_i["y1"] = new_y1;
      }
      y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
      y_rep[i] = y_rep_i;
    }

    arma::ivec length_inter(s_rep_size);

    arma::mat temp_y;
    for (unsigned int i = 0; i < s_rep_size; ++i){
      const Rcpp::List & y_rep_i = y_rep[i];
      if (use0)
        temp_y = Rcpp::as<arma::mat>(y_rep_i["y0"]);
      else
        temp_y = Rcpp::as<arma::mat>(y_rep_i["y1"]);
      arma::uvec non_na_indices = find_nan(temp_y.col(0));
      length_inter(i) = temp_y.col(0).n_elem - non_na_indices.n_elem;
    }

    arma::uvec valid = length_inter >= c_k;
    if (arma::accu(valid) == 0) {
      valid.elem(arma::find(length_inter == arma::max(length_inter))).fill(1);
    }

    double min_d = std::numeric_limits<double>::max();
    int min_s = 0;

    for (unsigned int i = 0; i < s_rep_size; i++) {
      if(valid(i)){
        double dist = Rcpp::as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
        if (dist < min_d)
        {
          min_d = dist;
          min_s = s_rep[i];
        }
      }
    }

    return arma::vec({static_cast<double>(min_s), min_d});
  }

  // [[Rcpp::export(.find_diss_aligned_rcpp)]]
  arma::vec find_diss_aligned_rcpp(const Rcpp::List &y,
                                   const Rcpp::List &v,
                                   const arma::vec & w,
                                   double alpha,
                                   bool aligned,
                                   unsigned int d,
                                   bool use0,
                                   bool use1,
                                   const Rcpp::Function & domain,
                                   const Rcpp::Function & select_domain,
                                   const Rcpp::Function & diss_d0_d1_L2)
  {
    Rcpp::LogicalVector v_dom = Rcpp::as<Rcpp::LogicalVector>(domain(v,use0));
    Rcpp::List v_new = select_domain(v, v_dom, use0, use1);
    int v_len = v_dom.size();
    int y_len = use0? Rcpp::as<arma::mat>(y[0]).n_rows :  Rcpp::as<arma::mat>(y[1]).n_rows;
    arma::ivec s_rep;
    if (aligned){
      s_rep = 1;
    } else {
      s_rep = arma::regspace<arma::ivec>(1, y_len - v_len + 1);
    }
    arma::uword s_rep_size = s_rep.size();
    Rcpp::List y_rep(s_rep_size);

    // Convert y[0] and y[1] to mat objects
    const int index_size = v_len;
    arma::mat temp_y0, temp_y1;
    if (use0) {
      temp_y0 = Rcpp::as<arma::mat>(y[0]);
    }
    if (use1) {
      temp_y1 = Rcpp::as<arma::mat>(y[1]);
    }

    auto index_range = std::views::iota(0,index_size);
    arma::ivec index;
    for (arma::uword i = 0; i < s_rep_size; ++i) {
      index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
      Rcpp::List y_rep_i = Rcpp::List::create(Rcpp::Named("y0") = R_NilValue, Rcpp::Named("y1") = R_NilValue);
      auto j_true = index_range
      | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len));});
      if (use0) {
        arma::mat new_y0(index_size, d);
        new_y0.fill(arma::datum::nan);
        std::for_each(j_true.begin(),j_true.end(),[&new_y0,&temp_y0,&index](int j){new_y0.row(j) = temp_y0.row(index[j] - 1);});
        y_rep_i["y0"] = new_y0;
      }
      if (use1) {
        arma::mat new_y1(index_size, d);
        new_y1.fill(arma::datum::nan);
        std::for_each(j_true.begin(),j_true.end(),[&new_y1,&temp_y1,&index](int j){new_y1.row(j) = temp_y1.row(index[j] - 1);});
        y_rep_i["y1"] = new_y1;
      }
      y_rep_i = select_domain(y_rep_i, v_dom, use0, use1);
      y_rep[i] = y_rep_i;
    }

    double min_d = std::numeric_limits<double>::max();
    int min_s = 0;

    for (arma::uword i = 0; i < s_rep_size; i++) {
      double dist = Rcpp::as<double>(diss_d0_d1_L2(y_rep[i], v_new, w, alpha));
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
                               bool use1,
                               const Rcpp::Function & domain,
                               const Rcpp::Function & select_domain,
                               const Rcpp::Function & diss_d0_d1_L2) {


  const unsigned int Y_size = Y.size();
  const unsigned int V_new_size = V_new.size();
  arma::vec sd(2);
  arma::imat S_new(Y_size,V_new_size);
  arma::mat  D_new(Y_size,V_new_size);
  for (unsigned int i = 0; i < V_new_size; ++i)
    for (unsigned int j = 0; j < Y_size; ++j)
    {
      sd = find_diss(Y[j],V_new[i],w,alpha,c_k(i),d,use0,use1,domain,select_domain,diss_d0_d1_L2);
      S_new(j,i) = sd(0);
      D_new(j,i) = sd(1);
    }
  return Rcpp::List::create(S_new,D_new);
}
