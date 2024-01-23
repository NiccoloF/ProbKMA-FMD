#' @title find_shift_warp_min
#'
#' @description Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
#'
#' @param y list of two elements y0=y(x), y1=y'(x), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x), matrices with d columns.
#' @param alpha weight coefficient between d0.L2 and d1.L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param c_k minimum length of supp(y_shifted) and supp(v) intersection.
#' @return Shift warping and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
find_shift_warp_min <- function(y,v,w,c_k,K,d,max_gap,alpha,use0,use1,domain,select_domain,diss_d0_d1_L2){
  
  out<-.find_shift_warp_min(y,v,w,c_k,K,d,max_gap,alpha,use0,use1,domain,select_domain,diss_d0_d1_L2)
}