#include "elongate_motifs.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

arma::mat compute_v_new_rcpp(const arma::field<arma::mat> & Y_inters_k,
                             const arma::umat & Y_inters_supp,
                             const arma::uvec & v_dom,
                             arma::uword v_len,
                             const arma::vec & p_k,
                             arma::uword d,
                             arma::uword m)
{
  arma::mat v_new(v_len,d,arma::fill::zeros);
  
  if (Y_inters_k.n_rows == 1){
    v_new.rows(arma::find(v_dom==1)) = Y_inters_k(0);
    v_new.rows(arma::find(v_dom==0)).fill(arma::datum::nan);
    return v_new;
  } 
  
  arma::vec coeff_k = arma::pow(p_k.elem(arma::find(p_k > 0)), m) / arma::sum(Y_inters_supp, 1);
  arma::vec coeff_x = arma::sum(arma::conv_to<arma::mat>::from(Y_inters_supp) % arma::repmat(coeff_k,1,Y_inters_supp.n_cols), 0).t();
  coeff_x.elem(arma::find(arma::sum(Y_inters_supp, 1) == 0)).fill(arma::datum::nan);
  
  for (arma::uword i = 0; i < Y_inters_k.n_rows; ++i){
    v_new.rows(arma::find(v_dom==1)) += (Y_inters_k(i)*coeff_k(i)) / arma::repmat(coeff_x,1,Y_inters_k(i).n_cols);
  }
  
  v_new.rows(arma::find(v_dom==0)).fill(arma::datum::nan);
  return v_new;
 
}
  
std::variant<std::pair<arma::field<arma::mat>,arma::sword>, arma::field<arma::mat>> compute_motif_rcpp(const arma::uvec & v_dom,
                                                                                                      const arma::ivec & s_k,
                                                                                                      const arma::vec & p_k,
                                                                                                      const arma::field<arma::mat> & Y,
                                                                                                      arma::uword m,
                                                                                                      bool use0,
                                                                                                      bool use1)
{
  if(arma::accu(p_k) == 0) {
    Rcpp::stop("Motif with no members! Degenerate cluster!");
  }
  arma::uword v_len = v_dom.n_elem; //length of the current domain
  arma::uword d = use0? Y(0,0).n_cols : Y(0,1).n_cols; //number of cols of the curves
  arma::uvec p_k_pos = arma::find(p_k > 0); //we will consider only the curves assigned to motif k
  arma::field<arma::mat> Y_inters_k(p_k_pos.n_elem,2); //list of shifted curves
  arma::umat Y_inters_supp(p_k_pos.n_elem,v_dom.n_elem); // for each shifted curve contains an uvec with the domain
  arma::mat y0(v_len,d); 
  arma::mat y1(v_len,d); 
  for (arma::uword i = 0; i < p_k_pos.n_elem; ++i){
    const int y_len = use0 ? Y(p_k_pos(i),0).n_rows : Y(p_k_pos(i),1).n_rows; //length of the curve
    arma::ivec index = std::max(1,s_k(p_k_pos(i))) - 1 + arma::regspace<arma::ivec>(1,v_len - std::max(0,1 - s_k(p_k_pos(i))));
    const int index_size = index.size();
    auto filtered_j = std::views::iota(0,index_size) // filtering of the indexes
      | std::views::filter([&y_len,&index](int j){return (index[j] <= y_len);});
    if (use0) {
      y0.fill(arma::datum::nan);
      for(int j : filtered_j) 
        y0.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),0).row(index(j) - 1);
      y0.shed_rows(arma::find(v_dom==0));
      Y_inters_supp.row(i) = domain_rcpp_base(y0);
      y0.replace(arma::datum::nan,0);
      Y_inters_k(i,0) = y0;
    }
    if (use1) {
      y1.fill(arma::datum::nan);
      for(int j : filtered_j) 
        y1.row(std::max(0,1 - s_k(p_k_pos(i))) + j) =  Y(p_k_pos(i),1).row(index(j) - 1); 
      y1.shed_rows(arma::find(v_dom==0));
      Y_inters_supp.row(i) = domain_rcpp_base(y1);
      y1.replace(arma::datum::nan,0);
      Y_inters_k(i,1) = y1;
    }
  }
  arma::field<arma::mat> v_new(1,2); 
  if (use0) {
    v_new(0,0) = compute_v_new_rcpp(Y_inters_k.col(0),
                                    Y_inters_supp,
                                    v_dom,
                                    v_len,
                                    p_k,
                                    d,
                                    m);
  }
  if (use1) {
    v_new(0,1) = compute_v_new_rcpp(Y_inters_k.col(1),
                                    Y_inters_supp,
                                    v_dom,
                                    v_len,
                                    p_k,
                                    d,
                                    m);
  }
  arma::uvec v_new_domain = find(domain_rcpp(v_new(0,0),v_new(0,1),use0) == 1);
  arma::sword index_min = v_new_domain.min();
  arma::sword index_max = v_new_domain.max();
  if (use0){
      v_new(0,0) = v_new(0,0).rows(index_min,index_max);
  }
  if (use1){
      v_new(0,1) = v_new(0,1).rows(index_min,index_max);
  }
  if (index_min > 1) {
    return std::make_pair(v_new,index_min - 1);
  }
  return v_new;
}

void elongation_rcpp(arma::field<arma::mat> & V_new, 
                     std::vector<arma::uvec> & V_dom,  
                     arma::imat & S_k, 
                     const arma::vec & p_k, 
                     const arma::ivec& len_elong_k, 
                     const arma::uvec& keep_k,  
                     double c, 
                     bool use0,
                     bool use1,
                     const arma::vec& w, 
                     double alpha, double max_gap,  
                     const arma::field<arma::mat> Y, 
                     int m, 
                     double deltaJk_elong,
                     const unsigned int index) 
{
  if(len_elong_k.empty()) return;

  const arma::field<arma::mat> & v_new_k = V_new.row(index);
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
  std::vector<arma::uvec> v_dom_elong_left_right(v_dom_elong_size);  
  const int v_dom_k_len = v_dom_k.n_elem;
  unsigned int k = 0;
  for (unsigned int i = 0; i < len_elong_k_zero_size; ++i){
    const arma::ivec & leng_elong_right_vector = len_elong_k_right_list[i];
    const unsigned int leng_elong_left = len_elong_k_zero(i);
    for (unsigned int j = 0; j <  leng_elong_right_vector.size(); ++j) {
      arma::uvec temp(leng_elong_left + v_dom_k_len + leng_elong_right_vector(j), arma::fill::ones);
      temp.rows(leng_elong_left, leng_elong_left + v_dom_k.n_elem - 1) = v_dom_k;
      v_dom_elong_left_right[k++] = temp;
    }
  }
  
  // create the list containing all the possible v_dom_k elongated using compute_motif
  const int v_elong_left_right_size = s_k_elong_left_right.size();
  arma::uvec not_start_with_NA(v_elong_left_right_size,arma::fill::zeros);
  arma::field<arma::mat> v_elong_left_right(v_elong_left_right_size,2); 
  for (int i = 0; i < v_elong_left_right_size; i++) {
    const auto & pair_motif_shift = compute_motif_rcpp(v_dom_elong_left_right[i+1], s_k_elong_left_right[i], p_k, Y, m, use0, use1);
    if (std::holds_alternative<arma::field<arma::mat>>(pair_motif_shift)){
      v_elong_left_right.row(i) = *(std::get_if<arma::field<arma::mat>>(&pair_motif_shift)); // errore
      not_start_with_NA(i) = 1;
    }
  }
  // filter centroid and shifts that are not in NA positions
  auto not_NA_index = std::views::iota(0,v_elong_left_right_size) 
    | std::views::filter([&not_start_with_NA](int index_j){return(not_start_with_NA(index_j));});
  
  arma::field<arma::mat> filtered_v_elong(arma::accu(not_start_with_NA),2);
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
  double Jk_before = compute_Jk_rcpp(v_new_k, s_k, p_k, Y, alpha, w, m, use0, use1, arma::datum::nan, arma::vec(arma::datum::nan));
  
  // compute performance indexes for all possible elongations
  arma::vec c_k_after(s_k_elong_left_right_size);
  arma::vec Jk_after(s_k_elong_left_right_size);
  
  for (arma::uword i = 0; i < s_k_elong_left_right_size; i++) {
    int c_i =  std::max(floor(domain_rcpp(v_elong_left_right(i,0),v_elong_left_right(i,1), use0).n_elem*(1 - max_gap)),c); 
    c_k_after[i] = c_i;
    Jk_after[i] = compute_Jk_rcpp(v_elong_left_right.row(i), s_k_elong_left_right[i], p_k, Y, alpha, w, m,use0,use1, static_cast<double>(c_i), arma::conv_to<arma::vec>::from(keep_k));
  }
  
  // find the best elongation in terms of perf. index
  arma::vec diff_perc = ((Jk_after-Jk_before)/Jk_before);
  arma::uword best_elong = arma::index_min(diff_perc);
  
  // check that the min really exists
  bool elongate = false;
  if (best_elong < diff_perc.size())
    elongate = diff_perc(best_elong) < deltaJk_elong;
  
  // evaluate if elongate or not
  if(elongate) {
    V_new.row(index) =  v_elong_left_right.row(best_elong);
    V_dom[index] =  v_dom_elong_left_right[best_elong + 1];
    S_k.col(index) = s_k_elong_left_right[best_elong];
  } else {
    return;
  }
}


//[[Rcpp::export(.elongate_motifs)]]
void elongate_motifs(Rcpp::List & V_new,
                     Rcpp::List & V_dom,
                     Rcpp::List & S_k,
                     const Rcpp::List & P_k,
                     const Rcpp::List & Y,
                     const arma::vec & w, 
                     int m, 
                     bool use0,
                     bool use1,
                     double alpha,
                     const arma::ivec & c,
                     const arma::ivec & c_max, // is a vector to be understood
                     double max_elong, 
                     double deltaJk_elong,
                     int trials_elong,
                     const arma::mat & D,
                     unsigned int K,
                     double max_gap)
 {
   // cast to imat for S_k, for each col we have a vector of shift
   arma::imat S_k_(Y.size(),K,arma::fill::zeros);
   for(unsigned int k=0; k < K; ++k){
    S_k_.col(k) = Rcpp::as<arma::ivec>(S_k[k]);
   }
   
   // cast to mat for P_k
   arma::mat P_k_(Y.size(),K,arma::fill::zeros);
   for(unsigned int k=0; k < K; ++k){
    P_k_.col(k) = Rcpp::as<arma::vec>(P_k[k]);
   }
   
   // cast to arma::field for Y
   auto Y_(util::conv_to_field(Y,use0,use1));
   
   // cast to arma::field for V_new
   auto V_new_(util::conv_to_field(V_new,use0,use1));
   
   // cast to std::vector<arma::uvec> for V_dom
   unsigned int V_dom_size = V_dom.size(); 
   std::vector<arma::uvec> V_dom_(V_dom_size);
   for (unsigned int i=0; i < V_dom_size; ++i){
    V_dom_[i] = Rcpp::as<arma::uvec>(V_dom[i]);
   }
   
   arma::ivec len_dom(V_dom_size);
   arma::uvec bool_gaps(V_dom_size);
   for(unsigned int i = 0; i < V_dom_size; ++i){
      len_dom(i) = V_dom_[i].size();
      bool_gaps(i) = (arma::accu(V_dom_[i]) < len_dom(i));
   }
   
   arma::uvec with_gaps = arma::find(bool_gaps == 1);
   unsigned int with_gaps_size = with_gaps.size();
   
   // @TODO: check this part for domains with NaN
   if (with_gaps_size != 0){
   
     arma::field<arma::mat> V_filled(with_gaps_size,2);
     
     std::vector<arma::uvec> V_dom_filled(with_gaps_size);
     
     arma::vec Jk_before(with_gaps_size);
     
     arma::vec Jk_after(with_gaps_size);
     
     // fill the domains of the motifs with gaps and recompute the motifs with the filled domains
     // and compute the perf.indexes before and after the filling
     for (unsigned int i = 0; i < with_gaps_size; ++i){
     
       V_dom_filled[i] = arma::uvec(V_dom_[with_gaps(i)].n_elem,arma::fill::ones);
       const auto & variant_motif = compute_motif_rcpp(V_dom_filled[i],  // errore
                                                        S_k_.col(with_gaps(i)),
                                                        P_k_.col(with_gaps(i)),
                                                        Y_,
                                                        m,
                                                        use0,
                                                        use1);
       V_filled.row(i) = *(std::get_if<arma::field<arma::mat>>(&variant_motif));
     }
     
     for (unsigned int i=0; i < with_gaps_size; ++i){
      Jk_before(i) = compute_Jk_rcpp(V_new_.row(with_gaps(i)), // va adattato
                                      S_k_.col(with_gaps(i)),
                                      P_k_.col(with_gaps(i)),
                                      Y_,
                                      alpha,
                                      w,
                                      m,
                                      use0,
                                      use1,
                                      arma::datum::nan,
                                      arma::vec(arma::datum::nan));
      
      Jk_after(i) = compute_Jk_rcpp(V_filled.row(i), // va adattato
                                    S_k_.col(with_gaps(i)),
                                    P_k_.col(with_gaps(i)),
                                    Y_,
                                    alpha,
                                    w,
                                    m,
                                    use0,
                                    use1,
                                    arma::datum::nan,
                                    arma::vec(arma::datum::nan));
     }
     // if filling the domain improves the perf. index over a certain threshold replace the domain and the motifs with the filled one
    const arma::uvec& fill = arma::find((Jk_after-Jk_before)/Jk_before < deltaJk_elong);
     for(unsigned int i=0; i < fill.size(); ++i){
       V_dom_[with_gaps(i)] = V_dom_filled[i];
       V_new_.row(with_gaps(i)) = V_filled.row(i);
     }

   }
   
   std::vector<arma::ivec> len_elong(V_dom_size);
   for (unsigned int i = 0; i < V_dom_size; ++i){
    const int len_max_elong_i = std::min<int>(std::floor(len_dom[i]*max_elong),c_max[i] - len_dom[i]);
     if (len_max_elong_i == 0){
       len_elong[i] = arma::ivec{};
     } 
     else{
       len_elong[i] =  (len_max_elong_i <= trials_elong) ? 
       arma::regspace<arma::ivec>(1, len_max_elong_i): 
       round(arma::linspace<arma::ivec>(1, len_max_elong_i, trials_elong));
    }
   }
   
   // vector of probabilities for the quantile function , to be checked this part
   arma::vec prob(1,arma::fill::value(0.25));
   // compute the quantile of the distance matrix
   arma::vec quantile = arma::quantile(vectorise(D), prob);
   // keep will be a matrix whose value i,j will be D(i,j) < quantile(0)
   arma::umat keep = D < quantile(0);
   // col-wise sum of the matrix keep
   const arma::uvec& col_sum_keep = (sum(keep, 0)).t();
   // vector of bool = true iff col_sum_keep[i]==0
   const arma::uvec& col_sum_keep_zero = (col_sum_keep==0);
   // empty_k stores the indexes of the col that have col_sum_keep = 0 
   const arma::uvec& empty_k = find(col_sum_keep_zero);
   
   for (auto k : empty_k){
   const unsigned int min_col_k_D = index_min(D.col(k));
   keep(min_col_k_D, k) = true;
   } 
   
   for (unsigned int i = 0; i < V_dom_size; ++i){
      elongation_rcpp(V_new_, 
                      V_dom_, 
                      S_k_,
                      P_k_.col(i),
                      len_elong[i],
                      keep.col(i),
                      c[i],
                      use0,
                      use1,
                      w, 
                      alpha,
                      max_gap,
                      Y_,
                      m,
                      deltaJk_elong,
                      i);  
   }
   // una volta modificato i V_new_, V_dom_,S_k_ devo ritrasformarli in V_new e in V_dom e in S_k
   for(unsigned int i=0; i < V_new.size(); ++i){
     V_new[i] = Rcpp::List::create(V_new_(i,0),V_new_(i,1));
     V_dom[i] = V_dom_[i];
     S_k[i] = S_k_.col(i);
   } 
 
 }
