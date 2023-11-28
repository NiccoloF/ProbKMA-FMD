#include "Motif_L2.hpp"

std::pair<arma::field<arma::mat>,arma::sword>
  Motif_L2::compute_motif(const arma::uvec& v_dom,
                          const arma::ivec& s_k,
                          const arma::vec& p_k,
                          const arma::field<arma::mat>& Y,
                          double m) const
 {
    if(arma::accu(p_k) == 0) {
      
      Rcpp::stop("Motif with no members! Degenerate cluster!");
    }
    
    arma::uword v_len = v_dom.n_elem; //length of the current domain
    
    arma::uword d = Y(0,0).n_cols; //number of cols of the curves
    
    arma::uvec p_k_pos = arma::find(p_k > 0); //we will consider only the curves assigned to motif k
    
    arma::field<arma::mat> Y_inters_k(p_k_pos.n_elem,2); //list of shifted curves
    
    arma::umat Y_inters_supp(p_k_pos.n_elem,v_dom.n_elem); // for each shifted curve contains an uvec with the domain
    
    for (arma::uword i = 0; i < p_k_pos.n_elem; ++i){
      
      const int y_len =  Y(p_k_pos(i),0).n_rows; //length of the curve
      
      arma::ivec index = std::max(1,s_k(p_k_pos(i))) - 1 + arma::regspace<arma::ivec>(1,v_len - std::max(0,1 - s_k(p_k_pos(i))));
      
      const int index_size = index.size();
      
      auto filtered_j = std::views::iota(0,index_size) // filtering of the indexes
        | std::views::filter([&y_len,&index](int j){return (index[j] <= y_len);});
      
      arma::mat y0(v_len,d); 
      
      y0.fill(arma::datum::nan);
      
      for(int j : filtered_j) // se questi sono consecutivi c'� modo migliore di agire
        y0.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),0).row(index(j) - 1);
      
      y0.shed_rows(arma::find(v_dom==0));
      
      Y_inters_supp.row(i) = util::findDomain<arma::mat>(y0);
      
      y0.replace(arma::datum::nan,0);
      
      Y_inters_k(i,0) = y0;
    }
    
    arma::field<arma::mat> v_new(1,2); 
    
    v_new(0,0) = compute_v_new(Y_inters_k.col(0),
                               Y_inters_supp,v_dom,
                               v_len,p_k, d,m);
                      
    arma::uvec v_new_domain = arma::find(util::findDomain<arma::mat>(v_new(0,0)) == 1);
    arma::sword index_min = v_new_domain.min();
    arma::sword index_max = v_new_domain.max();
    
    v_new(0,0) = v_new(0,0).rows(index_min,index_max);
    
    if (index_min > 1) {
      return std::make_pair(v_new,index_min - 1);
    }
    return std::make_pair(v_new,arma::datum::nan);
    // slight different to the R case in both case we return two elements in this case, see in elongate_motifs
 }


RCPP_MODULE(MotifL2Module) {
  Rcpp::class_<Motif_L2>("Motif_L2")
  .constructor();
  //.method("compute_motif", &Motif_L2::compute_motif);
}
    
    
    