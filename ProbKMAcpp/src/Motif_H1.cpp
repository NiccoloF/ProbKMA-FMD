#include "Motif.hpp"

std::variant<MotifPure::indexField,KMA::Mfield>
  Motif_H1::compute_motif(const arma::urowvec& v_dom,
                          const KMA::ivector& s_k,
                          const KMA::vector& p_k,
                          const KMA::Mfield& Y,
                          double m) const
{
  if(arma::accu(p_k) == 0) {
    Rcpp::stop("Motif with no members! Degenerate cluster!");
  }
  
  arma::uword v_len = v_dom.n_elem; //length of the current domain
  
  arma::uword d = Y(0,0).n_cols; //number of cols of the curves
  
  KMA::uvector p_k_pos = arma::find(p_k > 0); //we will consider only the curves assigned to motif k
  
  KMA::Mfield Y_inters_k(p_k_pos.n_elem,2); //list of shifted curves
 
  KMA::umatrix Y_inters_supp(p_k_pos.n_elem,v_dom.n_elem); // for each shifted curve contains an uvec with the domain
  
  for (arma::uword i = 0; i < p_k_pos.n_elem; ++i){
    
    const int y_len = Y(p_k_pos(i),0).n_rows; //length of the curve
    
    KMA::ivector index = std::max(1,s_k(p_k_pos(i))) - 1 + arma::regspace<KMA::ivector>(1,v_len - std::max(0,1 - s_k(p_k_pos(i))));
    
    const int index_size = index.size();
    
    auto filtered_j = std::views::iota(0,index_size) // filtering of the indexes
      | std::views::filter([&y_len,&index](int j){return (index[j] <= y_len);});
  
    KMA::matrix y0(v_len,d);

    y0.fill(arma::datum::nan);

    for(int j : filtered_j) // se questi sono consecutivi c'è modo migliore di agire
      y0.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),0).row(index(j) - 1);
    
    y0.shed_rows(arma::find(v_dom==0));

    Y_inters_supp.row(i) = util::findDomain<KMA::matrix>(y0);
    
    y0.replace(arma::datum::nan,0);
    
    Y_inters_k(i,0) = y0;
    
    KMA::matrix y1(v_len,d); 
    
    y1.fill(arma::datum::nan);
  
    for(int j : filtered_j) // se questi sono consecutivi c'è modo migliore di agire
      y1.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),1).row(index(j) - 1); 
    
    y1.shed_rows(arma::find(v_dom==0));
    Y_inters_supp.row(i) = util::findDomain<KMA::matrix>(y1);
    
    y1.replace(arma::datum::nan,0);
    
    Y_inters_k(i,1) = y1;
    
  }
  
  KMA::Mfield v_new(1,2); 
  v_new(0,0) = compute_v_new(Y_inters_k.col(0),
                             Y_inters_supp,v_dom,
                             v_len,p_k,d,m);
  
  v_new(0,1) = compute_v_new(Y_inters_k.col(1),
                             Y_inters_supp,v_dom,
                             v_len,p_k,d,m);
  
  KMA::uvector v_new_domain = arma::find(util::findDomain<KMA::matrix>(v_new(0,0)) == 1);
  arma::sword index_min = v_new_domain.min();
  arma::sword index_max = v_new_domain.max();
    
  v_new(0,0) = v_new(0,0).rows(index_min,index_max);
 
  v_new(0,1) = v_new(0,1).rows(index_min,index_max);
 
  if (index_min > 1) {
    return std::make_pair(v_new,index_min - 1);
  }
  return v_new;
  // slight different to the R case in both case we return two elements in this case, see in elongate_motifs
}  

RCPP_MODULE(MotifH1Module) {
  Rcpp::class_<Motif_H1>("Motif_H1")
  .constructor();
  //.method("compute_motif", &Motif_H1::compute_motif);
}