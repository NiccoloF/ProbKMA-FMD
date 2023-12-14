#include "Dissimilarity.hpp"

template<bool use1>
KMA::vector SobolDiss::find_diss_helper(const KMA::Mfield Y,
                                        const KMA::Mfield V,
                                        const KMA::vector& w, 
                                        double alpha, unsigned int c_k) const
    {
      // Convert domain and select_domain
      unsigned int d = Y(0,0).n_cols;
      arma::urowvec v_dom = util::findDomain(V(0,0));
      const KMA::Mfield& v_new = util::selectDomain<use1,arma::urowvec>(v_dom,V);
      int v_len = v_dom.size();
      int y_len = Y(0,0).n_rows;
      arma::ivec s_rep = arma::regspace<arma::ivec>(1 - (v_len - c_k), y_len - v_len + 1 + (v_len - c_k));
      const unsigned int s_rep_size = s_rep.size();
      KMA::Mfield y_rep(s_rep_size,Y.n_cols);
      
      const int index_size = v_len;
      auto index_range = std::views::iota(0,index_size);
      arma::uvec indeces_dom = arma::find(v_dom==0);
      for (unsigned int i = 0; i < s_rep_size; ++i) {
        arma::ivec index = s_rep(i) - 1 + arma::regspace<arma::ivec>(1,v_len);
        auto j_true = index_range
        | std::views::filter([&index,&y_len](int j){return((index[j] > 0) && (index[j] <= y_len)) ;}); 
        
        y_rep(i,0).set_size(index_size, d);
        y_rep(i,0).fill(arma::datum::nan);
        std::for_each(j_true.begin(),j_true.end(),[&y_rep,&Y,&index,&i](int j){y_rep(i,0).row(j) = Y(0,0).row(index[j] - 1);});
        y_rep(i,0).shed_rows(indeces_dom); 
        
        if constexpr(use1) {
          y_rep(i,1).set_size(index_size, d);
          y_rep(i,1).fill(arma::datum::nan);
          std::for_each(j_true.begin(),j_true.end(),[&y_rep,&Y,&index,&i](int j){y_rep(i,1).row(j) = Y(0,1).row(index[j] - 1);});
          y_rep(i,1).shed_rows(indeces_dom);
        }
      }
      
      KMA::ivector length_inter(s_rep_size);
      KMA::uvector non_na_indices;
      for (unsigned int i = 0; i < s_rep_size; ++i){
        non_na_indices = arma::find_nan(y_rep(i,0).col(0));
        length_inter(i) = y_rep(i,0).col(0).n_elem - non_na_indices.n_elem;
      }
      arma::uvec valid = length_inter >= c_k;
      if (arma::accu(valid) == 0) {        
        valid.elem(arma::find(length_inter == arma::max(length_inter))).fill(1);
      }
      
      double min_d = std::numeric_limits<double>::max();
      int min_s = 0;
    
      for (unsigned int i = 0; i < s_rep_size; i++) {
        if (valid(i)) {
        const double dist = computeDissimilarity(y_rep.row(i),v_new);
        Rcpp::Rcout<<"dist="<<dist<<std::endl;
          if (dist < min_d){
            min_d = dist;
            min_s = s_rep[i];
          } 
        }
      }
      
      return arma::vec({static_cast<double>(min_s), min_d}); 
  }