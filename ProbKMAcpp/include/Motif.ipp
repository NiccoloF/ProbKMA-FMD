#ifndef __MOTIF_IPP__
#define __MOTIF_IPP__
#include "Motif.hpp"

template<bool use1>
std::variant<MotifPure::indexField,KMA::Mfield> 
MotifSobol::compute_motif_helper(const arma::urowvec& v_dom,
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
  
  KMA::Mfield Y_inters_k(p_k_pos.n_elem,1); //list of shifted curves
  // If H1
  if constexpr(use1)
  {
    Y_inters_k.set_size(p_k_pos.n_elem,2);
  }
  
  KMA::umatrix Y_inters_supp(p_k_pos.n_elem,v_dom.n_elem); // for each shifted curve contains an uvec with the domain
  
  for (arma::uword i = 0; i < p_k_pos.n_elem; ++i){
    
    const int y_len = Y(p_k_pos(i),0).n_rows;//length of the curve
    
    KMA::ivector index = std::max(1,s_k(p_k_pos(i))) - 1 + arma::regspace<KMA::ivector>(1,v_len - std::max(0,1 - s_k(p_k_pos(i))));
    
    const int index_size = index.size();
    
    auto filtered_j = std::views::iota(0,index_size) // filtering of the indexes
      | std::views::filter([&y_len,&index](int j){return (index[j] <= y_len);});
    
    KMA::matrix y0(v_len,d); 
    
    y0.fill(arma::datum::nan);
    
    for(int j : filtered_j) // se questi sono consecutivi c'� modo migliore di agire
      y0.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),0).row(index(j) - 1);
      
    y0.shed_rows(arma::find(v_dom==0));
    
    Y_inters_supp.row(i) = util::findDomain(y0);
      
    y0.replace(arma::datum::nan,0);
    
    Y_inters_k(i,0) = y0;
    
    if constexpr(use1) {
      KMA::matrix y1(v_len,d); 
      
      y1.fill(arma::datum::nan);
      
      for(int j : filtered_j) // se questi sono consecutivi c'è modo migliore di agire
        y1.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),1).row(index(j) - 1); 
        
      y1.shed_rows(arma::find(v_dom==0));
      
      Y_inters_supp.row(i) = util::findDomain(y1);
      
      y1.replace(arma::datum::nan,0);
      
      Y_inters_k(i,1) = y1;
    }
  }
  
  KMA::Mfield v_new(1,1);
  
  if constexpr (use1) {
    v_new.set_size(1,2); // resize
    v_new(0,1) = compute_v_new(Y_inters_k.col(1),
                               Y_inters_supp,
                               v_dom,v_len,p_k,d,m);
  }
  
  v_new(0,0) = compute_v_new(Y_inters_k.col(0),
                               Y_inters_supp,
                               v_dom,v_len,p_k,d,m);

  
  KMA::uvector v_new_domain = arma::find(util::findDomain<KMA::matrix>(v_new(0,0)) == 1);
  arma::sword index_min = v_new_domain.min();
  arma::sword index_max = v_new_domain.max();

  v_new(0,0) = v_new(0,0).rows(index_min,index_max);
  
  if constexpr (use1){
      v_new(0,1) = v_new(0,1).rows(index_min,index_max);
  }
  if (index_min > 1) {
    return std::make_pair(v_new,index_min - 1);
  }
  return std::make_pair(v_new,arma::datum::nan);
  // slight different to the R case in both case we return two elements in this case, see in elongate_motifs                     
}

template<bool use1>
void MotifSobol::elongation(KMA::Mfield& V_new, 
                            std::vector<KMA::uvector> & V_dom,  
                            KMA::imatrix & S_k, 
                            const arma::vec & p_k, 
                            const arma::ivec& len_elong_k, 
                            const arma::uvec& keep_k,
                            double c,
                            const KMA::Mfield Y, 
                            const unsigned int index,
                            const Parameters& param,
                            const std::shared_ptr<PerformanceIndexAB>& perf,
                            const std::shared_ptr<Dissimilarity>& diss) const
{
  if(len_elong_k.empty()) return;
  
  const KMA::Mfield & v_new_k = V_new.row(index);
  const arma::uvec & v_dom_k = V_dom[index];
  const arma::ivec & s_k = S_k.col(index);
  
  // new vec with zero at the top
  arma::ivec len_elong_k_zero(len_elong_k.size() + 1, arma::fill::zeros);
  std::copy(len_elong_k.begin(), len_elong_k.end(), len_elong_k_zero.begin() + 1);
  
  // create a matrix whose column_i contains the vector s_k - len_elong_k_zero[i]
  arma::uword len_elong_k_zero_size = len_elong_k_zero.size();
  arma::imat s_k_elong_left_right_temp(s_k.n_elem, len_elong_k_zero_size);
  
  //s_k_elong_left_right_temp.each_col() += s_k -  len_elong_k_zero(i);
  for (arma::uword i=0; i < len_elong_k_zero_size;++i) {
    s_k_elong_left_right_temp.col(i) = s_k - len_elong_k_zero(i);
  }
  
  // create a sequence of integer from len_elong_k_zero.size() to 1
  arma::ivec reversedSequence = arma::regspace<arma::ivec>(len_elong_k_zero_size,-1,1);
  reversedSequence(0) -= 1;
  
  // repeat each col of s_k_elong_left_right a number of times specified by reversedSequence 
  std::vector<arma::ivec> s_k_elong_left_right = util::repeat_elements(s_k_elong_left_right_temp, reversedSequence);
  
  std::vector<arma::ivec> len_elong_k_right_list(len_elong_k_zero_size);
  const int max_len_elong_k = len_elong_k.back();
  unsigned int v_dom_elong_size = 0;
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i) {
    len_elong_k_right_list[i] = len_elong_k_zero.elem(find(len_elong_k_zero <=  max_len_elong_k - len_elong_k_zero(i)));
    v_dom_elong_size += len_elong_k_right_list[i].size();
  }
  
  //  v_dom_elong_left_right will be a vector of arma::uvec containing all the elongated version of v_dom_k
  std::vector<arma::urowvec> v_dom_elong_left_right(v_dom_elong_size);  
  const int v_dom_k_len = v_dom_k.n_elem;
  unsigned int k = 0;
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i){
    const arma::ivec & leng_elong_right_vector = len_elong_k_right_list[i];
    const unsigned int leng_elong_left = len_elong_k_zero(i);
    for (unsigned int j = 0; j <  leng_elong_right_vector.size(); ++j) {
      arma::urowvec temp(leng_elong_left + v_dom_k_len + leng_elong_right_vector(j), arma::fill::ones);
      temp.rows(leng_elong_left, leng_elong_left + v_dom_k.n_elem - 1) = v_dom_k;
      v_dom_elong_left_right[k++] = temp;
    }
  }
  
  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  arma::uvec not_start_with_NA(v_elong_left_right_size,arma::fill::zeros);
  KMA::Mfield v_elong_left_right(v_elong_left_right_size,1);
  if constexpr(use1)
  {
    v_elong_left_right.set_size(v_elong_left_right_size,2);
  }
  for (int i = 0; i < v_elong_left_right_size; i++) {
    const auto & pair_motif_shift = compute_motif_helper<use1>(v_dom_elong_left_right[i+1], s_k_elong_left_right[i], p_k, Y, param._m);
    if (std::holds_alternative<KMA::Mfield>(pair_motif_shift)){
      v_elong_left_right.row(i) = *(std::get_if<KMA::Mfield>(&pair_motif_shift)); 
      not_start_with_NA(i) = 1;
    }
  }
  // filter centroid and shifts that are not in NA positions
  auto not_NA_index = std::views::iota(0,v_elong_left_right_size) 
    | std::views::filter([&not_start_with_NA](int index_j){return(not_start_with_NA(index_j));});
  
  arma::field<arma::mat> filtered_v_elong(arma::accu(not_start_with_NA),1);
  if constexpr(use1)
  {
    v_elong_left_right.set_size(arma::accu(not_start_with_NA),2);
  }
  
  int i = 0;
  for(const auto & index : not_NA_index){
    filtered_v_elong.row(i++) = v_elong_left_right.row(index);
  }
  
  auto filtered_s_k = not_NA_index 
  | std::views::transform([&s_k_elong_left_right](int j){return s_k_elong_left_right[j];});
  
  v_elong_left_right = filtered_v_elong;
  
  s_k_elong_left_right = std::vector<arma::ivec>(filtered_s_k.begin(),filtered_s_k.end());
  const unsigned int s_k_elong_left_right_size = s_k_elong_left_right.size();
  
  // compute performance index before elongation
  double Jk_before = perf->compute_Jk(v_new_k,s_k,p_k,Y,param._w,param._m,
                                      arma::datum::nan, KMA::vector(arma::datum::nan), 
                                      diss);

  // compute performance indexes for all possible elongations
  arma::vec c_k_after(s_k_elong_left_right_size);
  arma::vec Jk_after(s_k_elong_left_right_size);
  
  for (arma::uword i = 0; i < s_k_elong_left_right_size; i++) {
   
    int c_i =  std::max(floor(util::findDomain(v_elong_left_right(i,0)).n_elem*(1 - param._max_gap)),c); 
    
    c_k_after[i] = c_i;
    Jk_after[i] = perf->compute_Jk(v_elong_left_right.row(i),s_k_elong_left_right[i],
                                   p_k,Y,param._w,param._m,static_cast<double>(c_i),
                                   arma::conv_to<KMA::vector>::from(keep_k),diss);
  }
  
  // find the best elongation in terms of perf. index
  arma::vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  arma::uword best_elong = arma::index_min(diff_perc);
  
  // check that the min really exists
  bool elongate = false;
  if (best_elong < diff_perc.size())
    elongate = diff_perc(best_elong) < param._deltaJK_elong;
  
  // evaluate if elongate or not
  if(elongate) {
    V_new.row(index) =  v_elong_left_right.row(best_elong);
    V_dom[index] =  v_dom_elong_left_right[best_elong + 1];
    S_k.col(index) = s_k_elong_left_right[best_elong];
  } else {
    return;
  }               
}


#endif 