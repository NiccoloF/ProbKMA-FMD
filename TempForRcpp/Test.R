set.seed(1)
standardize = TRUE
diss = 'd0_d1_L2' # try with d1_L2 d0_d1_L2 d0_L2
alpha = 0.5 #                 0.0    0.5     0.0
max_gap = 0 
trials_elong = 200
c_max = 53
K = 2
c = 40 
standardize=FALSE
c=c
c_max=c_max
N = 43
P0= matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
S0= matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
diss=diss
alpha=alpha
w=c(0.5,0.5)
m=2
iter_max=1
stop_criterion='max'
quantile=0.25
tol=1e-8
iter4elong=10
tol4elong=100
max_elong=0.5
trials_elong=10
deltaJk_elong=0.05
max_gap=max_gap
iter4clean=50
tol4clean=1e-4
quantile4clean=1/K
return_options=TRUE
prob <- 0.5
#return_init=TRUE


load("../TempForRcpp/Y_data.Rdata")


params <- list(standardize=standardize, K=K,c = c,c_max = c_max,iter_max = iter_max,
               quantile = quantile,stopCriterion = stop_criterion,tol = tol,
               iter4elong = iter4elong,tol4elong = tol4elong,max_elong = max_elong,
               trials_elong = trials_elong,
               deltaJK_elong = deltaJk_elong,max_gap = max_gap,iter4clean = iter4clean,
               tol4clean = tol4clean,
               quantile4clean = quantile4clean,return_options = return_options,
               m = m,w = w,alpha = alpha,prob = prob)

#library(ProbKMAcpp)
Y0_f <- function(Y_i)
{
  return(Y_i$y0)
}
Y1_f <- function(Y_i)
{
  return(Y_i$y1)
}

library(ProbKMAcpp)
Y0 <- lapply(Y,Y0_f)
Y1 <- lapply(Y,Y1_f)

a <- ProbKMAcpp::initialChecks(Y0,Y1,P0,S0,params,diss,alpha,w)
params <- a$Parameters
data <- a$FuncData

prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"H1")
b <- prok$probKMA_run()

.mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}
.diss_d0_d1_L2 <- function(y,v,w,alpha){
  # Dissimilarity index for multidimensional curves (dimension=d).
  # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
  # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
  
  .diss_L2 <- function(y,v,w){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # L2 distance with normalization on common support.
    # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    
    sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
  }
  
  if(alpha==0){
    return(.diss_L2(y[[1]],v[[1]],w))
  }else if(alpha==1){
    return(.diss_L2(y[[2]],v[[2]],w))
  }else{
    return((1-alpha)*.diss_L2(y[[1]],v[[1]],w)+alpha*.diss_L2(y[[2]],v[[2]],w))
  }
}
.domain <- function(v,use0){
  if(use0){
    rowSums(!is.na(v[[1]]))!=0
  }else{
    rowSums(!is.na(v[[2]]))!=0
  }
}
.select_domain <- function(v,v_dom,use0,use1){
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}

.find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1){
  # Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
  # Return shift and dissimilarity.
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  
  v_dom=.domain(v,use0)
  v=.select_domain(v,v_dom,use0,use1)
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1]
  s_rep=(1-(v_len-c_k)):(y_len-v_len+1+(v_len-c_k))
  index = seq_len(v_len)
  y_rep=lapply(s_rep,
               function(i){
                 index=i-1+seq_len(v_len)
                 y_rep_i=list(y0=NULL,y1=NULL)
                 if(use0){
                   y_rep_i$y0=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[1]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 if(use1){
                   y_rep_i$y1=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[2]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 y_rep_i=.select_domain(y_rep_i,v_dom,use0,use1)
                 return(y_rep_i)
               })
  
  length_inter=unlist(lapply(y_rep,
                             function(y_rep_i){
                               if(use0)
                                 return(sum((!is.na(y_rep_i$y0[,1]))))
                               return(sum((!is.na(y_rep_i$y1[,1]))))
                             }))
  valid=length_inter>=c_k
  if(sum(valid)==0){
    valid[length_inter==max(length_inter)]=TRUE
  }
  s_rep=s_rep[valid]
  y_rep=y_rep[valid]
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha)))
  return(c(s_rep[which.min(d_rep)],min(d_rep)))
}


.find_diss <- function(y,v,alpha,w,aligned,d,use0,use1){
  # Find dissimilarity between multidimensional curves (dimension=d), without alignment unless their lengths are different.
  # Return shift and dissimilarity.
  # To be used by probKMA_silhouette fucntion.
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # aligned: if TRUE, curves are already aligned. If FALSE, the shortest curve is aligned inside the longest.
  
  v_dom=.domain(v,use0)
  v=.select_domain(v,v_dom,use0,use1)
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1]
  if(aligned){
    s_rep=1
  }else{
    s_rep=1:(y_len-v_len+1)
  }
  y_rep=lapply(s_rep,
               function(i){
                 index=i-1+seq_len(v_len)
                 y_rep_i=list(y0=NULL,y1=NULL)
                 if(use0){
                   y_rep_i$y0=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[1]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 if(use1){
                   y_rep_i$y1=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[2]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 y_rep_i=.select_domain(y_rep_i,v_dom,use0,use1)
                 return(y_rep_i)
               })
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha)))
  return(c(s_rep[which.min(d_rep)],min(d_rep)))
}


.compute_motif <- function(v_dom,s_k,p_k,Y,m,use0,use1){
  # Compute the new motif v_new.
  # v_dom: TRUE for x in supp(v).
  # s_k: shift vector for motif k.
  # p_k: membership vector for motif k.
  # Y: list of N lists of two elements, Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
  
  .domain <- function(v,use0){
    if(use0){
      rowSums(!is.na(v[[1]]))!=0
    }else{
      rowSums(!is.na(v[[2]]))!=0
    }
  }
  .select_domain <- function(v,v_dom,use0,use1){
    if(use0)
      v[[1]]=as.matrix(v[[1]][v_dom,])
    if(use1)
      v[[2]]=as.matrix(v[[2]][v_dom,])
    return(v)
  }
  .compute_v_new <- function(Y_inters_k,Y_inters_supp,v_dom,v_len,p_k,d,m){
    v_new=matrix(NA,nrow=v_len,ncol=d)
    if(length(Y_inters_k)==1){
      v_new[v_dom,]=Y_inters_k[[1]]
      return(v_new)
    }
    Y_inters_supp=Reduce(rbind,Y_inters_supp)
    coeff_k=p_k[p_k>0]^m/(rowSums(Y_inters_supp)) # NB: divide for the length of the interval, not for the squared length!
    coeff_x=colSums(Y_inters_supp*coeff_k)
    coeff_x[colSums(Y_inters_supp)==0]=NA
    v_new[v_dom,]=Reduce('+',mapply('*',Y_inters_k,coeff_k,SIMPLIFY=FALSE))/coeff_x
    return(v_new)
  }
  
  if(sum(p_k)==0){
    stop('Motif with no members! Degenerate cluster!')
  }
  v_len=length(v_dom)
  d=unlist(lapply(Y[[1]],ncol))[1] # dimension of curves
  Y_inters_k=mapply(function(y,s_k_i,d,use0,use1){
    y_len=unlist(lapply(y,nrow))[1]
    y_inters_k=list(y0=NULL,y1=NULL)
    index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
    if(use0)
      y_inters_k$y0=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                          matrix(y$y0[index[index<=y_len],],ncol=d),
                          matrix(NA,nrow=sum(index>y_len),ncol=d))
    if(use1)
      y_inters_k$y1=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                          matrix(y$y1[index[index<=y_len],],ncol=d),
                          matrix(NA,nrow=sum(index>y_len),ncol=d))
    return(.select_domain(y_inters_k,v_dom,use0,use1))
  },Y[p_k>0],s_k[p_k>0],MoreArgs=list(d,use0,use1),SIMPLIFY=FALSE)
  Y_inters_supp=lapply(Y_inters_k,.domain,use0)
  Y_inters_k=mapply(function(y_inters_k,y_inters_supp,use0,use1){
    if(use0)
      y_inters_k$y0[!y_inters_supp]=0
    if(use1)
      y_inters_k$y1[!y_inters_supp]=0
    return(y_inters_k)},
    Y_inters_k,Y_inters_supp,SIMPLIFY=FALSE,MoreArgs=list(use0,use1))
  v_new=list(v0=NULL,v1=NULL)
  if(use0)
    v_new$v0=.compute_v_new(lapply(Y_inters_k,function(y_inters_k) y_inters_k$y0),
                            Y_inters_supp,v_dom,v_len,p_k,d,m)
  if(use1)
    v_new$v1=.compute_v_new(lapply(Y_inters_k,function(y_inters_k) y_inters_k$y1),
                            Y_inters_supp,v_dom,v_len,p_k,d,m)
  # check if there are NA at the corners (it can happen after cleaning), and remove it
  range_v_new=range(which(.domain(v_new,use0)))
  v_dom_new=rep(FALSE,v_len)
  v_dom_new[(range_v_new[1]):(range_v_new[2])]=TRUE
  v_new=.select_domain(v_new,v_dom_new,use0,use1)
  if(range_v_new[1]>1){
    return(c(v_new,list(shift=range_v_new[1]-1)))
  }
  return(v_new)
}


.compute_Jk <- function(v,s_k,p_k,Y,alpha,w,m,use0,use1,c_k=NULL,keep_k=NULL){
  # Compute the objective function J for the motif k.
  # v: list of two elements, v0=v(x), v1=v'(x), matrices with d columns.
  # s_k: shift vector for motif k.
  # p_k: membership vector for motif k.
  # Y: list of N lists of two elements, Y0=y_i(x), Y1=y'_i(x), matrices with d columns, for d-dimensional curves.
  # alpha: weight coefficient between d0_L2 and d1_L2.
  # m>1: weighting exponent in least-squares functional.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  # keep_k: check c_k only when keep=TRUE for y_shifted.
  
  .diss_d0_d1_L2 <- function(y,v,w,alpha){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
    # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
    # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
    
    .diss_L2 <- function(y,v,w){
      # Dissimilarity index for multidimensional curves (dimension=d).
      # L2 distance with normalization on common support.
      # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
      # w: weights for the dissimilarity index in the different dimensions (w>0).
      
      sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y) # NB: divide for the length of the interval, not for the squared length!
    }
    
    if(alpha==0){
      return(.diss_L2(y[[1]],v[[1]],w))
    }else if(alpha==1){
      return(.diss_L2(y[[2]],v[[2]],w))
    }else{
      return((1-alpha)*.diss_L2(y[[1]],v[[1]],w)+alpha*.diss_L2(y[[2]],v[[2]],w))
    }
  }
  .domain <- function(v,use0){
    if(use0){
      rowSums(!is.na(v[[1]]))!=0
    }else{
      rowSums(!is.na(v[[2]]))!=0
    }
  }
  .select_domain <- function(v,v_dom,use0,use1){
    if(use0)
      v[[1]]=as.matrix(v[[1]][v_dom,])
    if(use1)
      v[[2]]=as.matrix(v[[2]][v_dom,])
    return(v)
  }
  
  v_dom=.domain(v,use0)
  v_len=length(v_dom)
  v=.select_domain(v,v_dom,use0,use1)
  d=unlist(lapply(Y[[1]],ncol))[1]
  Y_inters_k=mapply(function(y,s_k_i,d){
    y_len=unlist(lapply(y,nrow))[1]
    y_inters_k=list(y0=NULL,y1=NULL)
    index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
    if(use0)
      y_inters_k$y0=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                          matrix(y$y0[index[index<=y_len],],ncol=d),
                          matrix(NA,nrow=sum(index>y_len),ncol=d))
    if(use1)
      y_inters_k$y1=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                          matrix(y$y1[index[index<=y_len],],ncol=d),
                          matrix(NA,nrow=sum(index>y_len),ncol=d))
    return(.select_domain(y_inters_k,v_dom,use0,use1))
  },Y,s_k,MoreArgs=list(d),SIMPLIFY=FALSE)
  if(!is.null(keep_k)){
    supp_inters_length=unlist(lapply(Y_inters_k[keep_k],function(y_inters_k) sum(.domain(y_inters_k,use0))))
    if(TRUE %in% (supp_inters_length<c_k))
      return(NA)
  }
  dist=unlist(mapply(.diss_d0_d1_L2,Y_inters_k,MoreArgs=list(v,w,alpha)))
  return(sum(dist*(p_k^m),na.rm=TRUE))
}

probKMA <- function(Y0,Y1=NULL,standardize=FALSE,K,c,c_max=Inf,P0=NULL,S0=NULL,
                    diss='d0_L2',alpha=NULL,w=1,m=2,
                    iter_max=1000,stop_criterion='max',quantile=NULL,tol=1e-8,
                    iter4elong=10,tol4elong=1e-3,max_elong=0.5,trials_elong=10,deltaJk_elong=0.05,max_gap=0.2,
                    iter4clean=50,tol4clean=1e-4,quantile4clean=1/K,
                    return_options=TRUE,return_init=TRUE,worker_number=NULL){
  # Probabilistic k-mean with local alignment to find candidate motifs.
  # Y0: list of N vectors, for univariate curves y_i(x), or
  #     list of N matrices with d columns, for d-dimensional curves y_i(x),
  #     with the evaluation of curves (all curves should be evaluated on a uniform grid).
  #     When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
  # Y1: list of N vectors, for univariate derivative curves y'_i(x), or
  #     list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
  #     with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
  #     When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
  #     Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
  # standardize: if TRUE, each dimension is standardized (Z score on all the regions together).
  # K: number of motifs.
  # c: minimum motif lengths. Can be an integer (or a vector of K integers).
  # c_max: maximum motif lengths. Can be an integer (or a vector of K integers).
  # P0: initial membership matrix, with N row and K column (if NULL, a random P0 is choosen).
  # S0: initial shift warping matrix, with N row and K column (if NULL, a random S0 is choosen).
  # diss: dissimilarity. Possible choices are 'd0_L2', 'd1_L2', 'd0_d1_L2'. 
  # alpha: when diss='d0_d1_L2', weight coefficient between d0_L2 and d1_L2 (alpha=0 means d0_L2, alpha=1 means d1_L2).
  # w: vector of weights for the dissimilarity index in the different dimensions (w>0).
  # m>1: weighting exponent in least-squares functional.
  # iter_max: the maximum number of iterations allowed.
  # stop_criterion: criterion to stop iterate, based on the Bhattacharyya distance between memberships in subsequent iterations. 
  #                 Possible choices are: 'max' for the maximum of distances in the different motifs;
  #                                       'mean' for the average of distances in the different motifs;
  #                                       'quantile' for the quantile of distances in the different motifs (in this case, quantile must be provided).
  # quantile: probability in [0,1] to be used if stop.criterion='quantile'.
  # tol: method tolerance (method stops if the stop criterion <tol).
  # iter4elong: motifs elongation is performed every iter4elong iterations (if iter4elong>iter.max, no elongation is done).
  # tol4elong: tolerance on the Bhattacharyya distance (with the choice made in stop.criterion) for performing motifs elongation.
  # max_elong: maximum elongation allowed in a single iteration, as percentage of the motif length.
  # trials_elong: number of (equispaced) elongation trials at each side of the motif in a single iteration.
  # deltaJk_elong: maximum relative objective function increasing allowed in each motif elongation (for gaps and each side).
  # max_gap: maximum gap allowed in each alignment (percentage of motif length).
  # iter4clean: motif cleaning is performed every iter4clean iterations (if iter4clean>iter_max, no cleaning is done).
  # tol4clean: tolerance on the Bhattacharyya distance (with the choice made in stop_criterion) for performing motif cleaning.
  # quantile4clean: dissimilarity quantile to be used in motif cleaning.
  # return_options: if TRUE, the options K,c,diss,w,m are returned by the function.
  # return_init: if TRUE, P0 and S0 are returned by the function.
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores -1). If worker_number=1, the function is run sequentially. 
  
  ### set parallel jobs #############################################################################
  start=proc.time()
  core_number <- detectCores()
  # check worker number
  if(!is.null(worker_number)){
    if(!is.numeric(worker_number)){
      warning('Worker number not valid. Selecting default number.')
      worker_number=NULL
    }else{
      if((worker_number%%1!=0)|(worker_number<1)|(worker_number>core_number)){
        warning('Worker number not valid. Selecting default number.')
        worker_number=NULL
      }
    }
  }
  if(is.null(worker_number))
    worker_number <- core_number-1
  rm(core_number)
  if(worker_number>1){
    cl_probKMA=makeCluster(worker_number,timeout=60*60*24*30)
    clusterExport(cl_probKMA,c('.diss_d0_d1_L2','.domain','.select_domain'))
    on.exit(stopCluster(cl_probKMA))
  }else{
    cl_probKMA=NULL
  }
  end=proc.time()
  #message('set parallel jobs: ',round((end-start)[3],2))
  
  ### check input ####################################################################################
  start=proc.time()
  # check dissimilarity
  if(!(diss %in% c('d0_L2','d1_L2','d0_d1_L2')))
    stop('Dissimilarity not valid. Possible choices are: \'d0_L2\', \'d1_L2\' and \'d0_d1_L2\'.')
  # check Y0
  if(missing(Y0))
    stop('Y0 must be specified.')
  if(!is.null(Y0))
    if(!is.list(Y0))
      stop('Y0 should be a list of vectors or matrices.')
  if((FALSE %in% lapply(Y0,is.matrix))&&(FALSE %in% lapply(Y0,is.vector)))
    stop('Y0 should be a list of vectors or matrices.')
  N=length(Y0) # number of curves
  if(N<5)
    stop('More curves y_i(x) needed.')
  Y0=lapply(Y0,as.matrix)
  d=unique(unlist(lapply(Y0,ncol),use.names=FALSE)) # dimension of curves
  if(length(d)>1)
    stop('Curves have not all the same dimension.')
  Y0_NA=lapply(Y0,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
  if(FALSE %in% Y0_NA){
    warning('y_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
    Y0=lapply(Y0,
              function(y,d){
                y[rowSums(!is.na(y))!=d,]=NA
                return(y)},
              d)
  }
  rm(Y0_NA)
  if(standardize){
    Y0_tot=Reduce(rbind,Y0)
    Y0_mean=colMeans(Y0_tot,na.rm=TRUE)
    Y0_sd=apply(Y0_tot,2,sd,na.rm=TRUE)
    rm(Y0_tot)
    Y0=lapply(Y0,function(y0,m,s) t((t(y0)-m)/s*100),m=Y0_mean,s=Y0_sd)
  }
  if(diss=='d0_d1_L2'){
    # check required input
    if(!is.numeric(alpha)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    }
    if((alpha<0)|(alpha>1)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    } else if(alpha==0){
      diss='d0_L2'
    } else if(alpha==1){
      diss='d1_L2'
    }
  }
  if(diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(Y0,function(y0) list(y0=y0,y1=NULL))
  }else{
    # check Y1
    if(is.null(Y1))
      stop(paste0('Y1 must be specified with diss=\'',diss,'\'.'))
    if(!is.list(Y1))
      stop('Y1 should be a list of vectors or matrices.')
    if((FALSE %in% lapply(Y1,is.matrix))&&(FALSE %in% lapply(Y1,is.vector)))
      stop('Y1 should be a list of vectors or matrices.')
    if(N!=length(Y1))
      stop('The number of curves in Y0 and derivatives Y1 should be the same.')
    Y1=lapply(Y1,as.matrix)
    d1=unique(unlist(lapply(Y1,ncol),use.names=FALSE)) # dimension of derivatives
    if(length(d1)>1)
      stop('Derivatives have not all the same dimension.')
    if(d!=d1)
      stop('The dimension of curves in Y0 and derivatives Y1 should be the same')
    Y1_NA=lapply(Y1,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
    if(FALSE %in% Y1_NA){
      warning('y\'_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
      Y1=lapply(Y1,
                function(y,d){
                  y[rowSums(!is.na(y))!=d,]=NA
                  return(y)},
                d)
    }
    rm(Y1_NA)
    Y0_NA=lapply(Y0,function(y) rowSums(is.na(y))!=0)
    Y1_NA=lapply(Y1,function(y) rowSums(is.na(y))!=0)
    diff_NA=mapply(function(y0,y1) y0!=y1,Y0_NA,Y1_NA)
    index_diff_NA=which(unlist(lapply(diff_NA,sum))!=0)
    if(length(index_diff_NA)>0){
      warning('y_j(x) and y\'_j(x) are not both defined, for some x. Putting NA in that case.')
      same_NA=mapply(function(y0,y1,diff_NA){
        y0[diff_NA]=NA
        y1[diff_NA]=NA
        return(list(y0,y1))
      },Y0[index_diff_NA],Y1[index_diff_NA],diff_NA[index_diff_NA])
      Y0[index_diff_NA]=same_NA[1,]
      Y1[index_diff_NA]=same_NA[2,]
    }
    rm(Y0_NA,Y1_NA,diff_NA)
    if(standardize){
      Y1=lapply(Y1,function(y1,s) t(t(y1)/s*100),s=Y0_sd)
    }
    if(diss=='d1_L2'){
      alpha=1
      use0=FALSE
      use1=TRUE
      Y=lapply(Y1,function(y1) list(y0=NULL,y1=y1))
    }
    if(diss=='d0_d1_L2'){
      use0=TRUE
      use1=TRUE
      Y=mapply(function(y0,y1) list(y0=y0,y1=y1),Y0,Y1,SIMPLIFY=FALSE)
    }
  }
  # check required input
  if(missing(K))
    stop('K must be specified')
  if(missing(c))
    stop('c must be specified')
  # check K
  if(length(K)!=1)
    stop('Number of motifs K not valid.')
  if(K%%1!=0)
    stop('Number of motifs K should be an integer.')
  if(K<1)
    stop('Number of motifs K should be at least 1.')
  # check c
  if(!(length(c) %in% c(1,K)))
    stop('Minimum motif length c not valid.')
  if(sum(c%%1!=0))
    stop('Minimum motif lengths should be integers.')
  if(sum(c<=1))
    stop('Minimum motif lengths should be at least 2.')
  c=rep(c,length.out=K)
  # find all intervals contained in supp(y_i) and their lengths
  Y_intervals=lapply(Y0,
                     function(y,d){
                       y_not_NA=(rowSums(!is.na(y))==d) # find supp(y)
                       intervals=which((y_not_NA[2:length(y_not_NA)]-y_not_NA[1:(length(y_not_NA)-1)])==1)+1
                       if(y_not_NA[1])
                         intervals=c(1,intervals)
                       intervals_lengths=rle(y_not_NA)$lengths[rle(y_not_NA)$values==TRUE]
                       intervals_all=data.frame(start=intervals,
                                                end=intervals+intervals_lengths-1,
                                                length=intervals_lengths)
                       return(intervals_all)
                     },d)
  # check minimum motif length compatibility
  min_length=unlist(lapply(Y_intervals,
                           function(y_intervals,c){
                             return(sum(y_intervals$length>=c))},
                           max(c)))
  if(0 %in% min_length)
    stop('Minimum motifs length is incompatible with supplied curves. Choose a smaller minimum motifs length.')
  # check c_max
  if(!(length(c_max) %in% c(1,K)))
    stop('Maximum motif length c_max not valid.')
  c_max=rep(c_max,length.out=K)
  for(k in 1:K){
    if(c_max[k]!=Inf){
      if(c_max[k]%%1!=0)
        stop('Maximum motif lengths should be integers.')
      if(c_max[k]<=1)
        stop('Maximum motif length should be at least 2.')
      # check maximum motif length compatibility
      max_length=unlist(lapply(Y_intervals,
                               function(y_intervals,c_max){
                                 return(sum(floor((1+max_gap)*y_intervals$length)>=c_max))},
                               c_max[k]))
      if(0 %in% max_length){
        warning('Maximum motif length is incompatible with supplied curves. Choosing default maximum motif length for motif ',k,'.')
        c_max[k]=Inf
      }
    }
    if(c_max[k]==Inf){
      c_max[k]=floor((1+max_gap)*min(unlist(lapply(Y_intervals,function(y_intervals) max(y_intervals$length)))))
    }
  }
  # check P0
  if(!is.null(P0)){
    if(sum(dim(P0)!=c(N,K))!=0){
      warning('Membership matrix dimensions not valid. Choosing random initial membership matrix.')
      P0=NULL
    }else{
      if(sum((P0<0)|(P0>1))){
        warning('Memberships should be non-negative and <=1. Choosing random initial membership matrix.')
        P0=NULL
      }else{
        if(sum(rowSums(P0)!=1)){
          warning('Memberships of each curve should sum to 1. Choosing random initial membership matrix.')
          P0=NULL
        }else{
          if(sum(colSums(P0)==0)){
            warning('Sum of memberships of each cluster should be positive. Choosing random initial membership matrix.')
            P0=NULL
          }
        }
      }
    }
  }
  # check S0
  if(!is.null(S0)){
    if(sum(dim(S0)!=c(N,K))!=0){
      warning('Shift warping matrix dimensions not valid. Choosing random initial shift warping matrix.')
      S0=NULL
    }else{
      Y_segments=mapply(function(y,S_i,c){
        return(mapply(function(s,c) NA %in% y[s+seq_len(c)-1],S_i,c))},
        Y0,lapply(seq_len(N),function(i) S0[i,]),MoreArgs=list(c),SIMPLIFY=TRUE)
      if(sum(Y_segments)){
        warning('Shift warping matrix not valid. Choosing random initial shift warping matrix.')
        S0=NULL
      }
    }
  }
  # check weigths
  if(!is.vector(w)|!is.numeric(w))
    stop('Weights w not valid.')
  if(!(length(w) %in% c(1,d)))
    stop('Weights w not valid.')
  if(TRUE %in% (w<=0))
    stop('Weights w not valid.')
  w=rep(w,length.out=d)
  # check weighting exponent
  if(!is.numeric(m)|(length(m)>1))
    stop('Weighting exponent m not valid.')
  if(m<=1)
    stop('Weighting exponent m should be >1')
  # check maximum number of iterations
  if(!is.numeric(iter_max)|(length(iter_max)!=1))
    stop('Maximum number of iterations iter_max not valid.')
  if(((iter_max%%1)!=0)|(iter_max<=0))
    stop('Maximum number of iterations iter_max not valid.')
  # check stop criterion
  if(!(stop_criterion %in% c('max','mean','quantile')))
    stop('Stop criterion not valid. Possible choices are: \'max\', \'mean\', \'quantile\'.')
  if(stop_criterion=='quantile')
    if((!is.numeric(quantile))||(length(quantile)!=1)||(quantile<0)||(quantile>1))
      stop('quantile should be a number in [0,1].')
  # check tolerance
  if((!is.numeric(tol))||(length(tol)!=1)||(tol<=0))
    stop('Tolerance should be a positive number.')
  # check elongation input
  if((!is.numeric(iter4elong))||(length(iter4elong)!=1))
    stop('iter4elong not valid.')
  if(((iter4elong%%1)!=0)|(iter4elong<=0))
    stop('iter4elong not valid.')
  if((!is.numeric(tol4elong))||(length(tol4elong)!=1)||(tol4elong<=0))
    stop('tol4elong should be a positive number.')
  if((!is.numeric(max_elong))||(length(max_elong)!=1)||(max_elong<=0))
    stop('max_elong should be a positive number.')
  if((!is.numeric(trials_elong))||(length(trials_elong)!=1))
    stop('trials_elong not valid.')
  if(((trials_elong%%1)!=0)|(trials_elong<=0))
    stop('trials_elong not valid.')
  if((!is.numeric(deltaJk_elong))||(length(deltaJk_elong)!=1)||(deltaJk_elong<=0))
    stop('deltaJk_elong should be a positive number.')
  if((!is.numeric(max_gap))||(length(max_gap)!=1)||(max_gap<0))
    stop('max_gap should be a non-negative number.')
  # check cleaning input
  if((!is.numeric(iter4clean))||(length(iter4clean)!=1))
    stop('iter4clean not valid.')
  if(((iter4clean%%1)!=0)|(iter4clean<=0))
    stop('iter4clean not valid.')
  if((!is.numeric(tol4clean))||(length(tol4clean)!=1)||(tol4clean<=0))
    stop('tol4clean should be a positive number.')
  if((!is.numeric(quantile4clean))||(length(quantile4clean)!=1)||(N*K*quantile4clean<K)||(quantile4clean>1))
    stop('quantile4clean not valid.')
  end=proc.time()
  #message('check input: ',round((end-start)[3],2))
  
  ### initialize #############################################################################################
  start=proc.time()
  if(is.null(P0)){
    # create random membership matrix, with N rows and k columns
    P0=matrix(c(runif(N*(K-1)),rep(1,N)),nrow=N,ncol=K)
    P0=as.matrix(apply(P0,1,sort))
    if(K>1)
      P0=cbind(P0[1,],Reduce('cbind',lapply(2:K,function(k) P0[k,]-P0[k-1,])))
  }
  colnames(P0)=NULL
  P=P0
  if(is.null(S0)){
    # create shift warping matrix, with N rows and k columns
    S0=matrix(unlist(lapply(Y_intervals,
                            function(y_intervals,c,K){
                              s0=rep(NA,K)
                              for(k in 1:K){
                                y_intervals_k=y_intervals[y_intervals$length>=c[k],]
                                y_starts=unlist(apply(y_intervals_k,1,
                                                      function(y_interval)
                                                        return(y_interval[1]:(y_interval[2]-c[k]+1))),
                                                use.names=FALSE)
                                s0[k]=sample(y_starts,1)
                              }
                              return(s0)
                            },c,K)),
              nrow=N,ncol=K,byrow=TRUE)
  }
  S=S0
  
  # create empty motifs
  V=lapply(c,
           function(c_k,d){
             v=list(v0=NULL,v1=NULL)
             if(use0)
               v$v0=matrix(0,nrow=c_k,ncol=d)
             if(use1)
               v$v1=matrix(0,nrow=c_k,ncol=d)
             return(v)
           },d)
  end=proc.time()
  #message('initialize: ',round((end-start)[3],2))
  
  ### iterate #############################################################################################
  iter=0
  J_iter=c()
  BC_dist_iter=c()
  BC_dist=+Inf
  while((iter<iter_max)&(BC_dist>tol)){
    iter=iter+1
    #message('iter ',iter)
    
    ##### clean motifs ###################################################################################
    start=proc.time()
    P_old=P
    if((iter>1)&&(!(iter%%iter4clean))&&(BC_dist<tol4clean)){
      keep=D<quantile(D,quantile4clean)
      empty_k=which(colSums(keep)==0)
      if(length(empty_k)>0){
        for(k in empty_k)
          keep[which.min(D[,k]),k]=TRUE
      }
      P[keep]=1
      P[!keep]=0
    }
    end=proc.time()
    #message('  clean: ',round((end-start)[3],2))
    
    ##### compute motifs ###################################################################################
    start=proc.time()
    S_k=split(S,rep(seq_len(K),each=N))
    P_k=split(P,rep(seq_len(K),each=N))
    V_dom=lapply(V,.domain,use0)
    V_new=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
    changed_s=which(unlist(lapply(V_new,length))>2)
    for(k in changed_s){
      S[,k]=S[,k]+V_new[[k]]$shift
      V_new[[k]]=V_new[[k]][c('v0','v1')]
    }
    S_k=split(S,rep(seq_len(K),each=N))
    V_dom=lapply(V_new,.domain,use0)
    end=proc.time()
    #message('  compute motifs: ',round((end-start)[3],2))
    
    ##### elongate motifs #################################################################################
    start=proc.time()
    if((iter>1)&&(!(iter%%iter4elong))&&(BC_dist<tol4elong)){
      # fill
      with_gaps=which(unlist(lapply(V_dom,function(v_dom) sum(!v_dom)!=0)))
      if(length(with_gaps)>0){
        V_dom_filled=lapply(V_dom[with_gaps],function(v_dom) rep_len(TRUE,length(v_dom)))
        V_filled=mapply(.compute_motif,V_dom_filled,S_k[with_gaps],P_k[with_gaps],MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
        Jk_before=mapply(.compute_Jk,
                         V_new[with_gaps],S_k[with_gaps],P_k[with_gaps],
                         MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))
        Jk_after=mapply(.compute_Jk,
                        V_filled,S_k[with_gaps],P_k[with_gaps],
                        MoreArgs=list(Y=Y,alpha=alpha,w=w,m=m,use0=use0,use1=use1))
        fill=(Jk_after-Jk_before)/Jk_before<deltaJk_elong
        V_dom[with_gaps[fill]]=V_dom_filled[fill]
        V_new[with_gaps[fill]]=V_filled[fill]
      }
      # elongate
      len_dom=unlist(lapply(V_dom,length))
      len_max_elong=mapply(min,floor(len_dom*max_elong),(c_max-len_dom))
      len_elong=lapply(len_max_elong,
                       function(len_max_elong){
                         if(len_max_elong<=trials_elong){
                           len_elong=seq_len(len_max_elong)
                         }else{
                           len_elong=round(seq(1,len_max_elong,length.out=trials_elong))
                         }
                         return(len_elong)
                       })
      # left and right elongation
      keep=D<quantile(D,0.25)
      empty_k=which(colSums(keep)==0)
      if(length(empty_k)>0){
        for(k in empty_k)
          keep[which.min(D[,k]),k]=TRUE
      }
      
      res_left_right=mapply(function(v_new_k,v_dom_k,s_k,p_k,len_elong_k,keep_k,c){
        if(length(len_elong_k)==0){
          return(list(v_new=v_new_k,
                      v_dom=v_dom_k,
                      s_k=s_k))
        }
        s_k_elong_left_right=rep(lapply(c(0,len_elong_k),function(len_elong_k) s_k-len_elong_k),(length(len_elong_k)+1):1)[-1]
        v_dom_elong_left_right=unlist(lapply(c(0,len_elong_k),
                                             function(len_elong_k_left)
                                               lapply(c(0,len_elong_k[len_elong_k<=(max(len_elong_k)-len_elong_k_left)]),
                                                      function(len_elong_k_right) 
                                                        c(rep_len(TRUE,len_elong_k_left),v_dom_k,rep_len(TRUE,len_elong_k_right)))),
                                      recursive=FALSE)[-1]
        v_elong_left_right=mapply(.compute_motif,v_dom_elong_left_right,s_k_elong_left_right,
                                  MoreArgs=list(p_k,Y,m,use0,use1),SIMPLIFY=FALSE)
        start_with_NA=unlist(lapply(v_elong_left_right,length))>2
        v_elong_left_right=v_elong_left_right[!start_with_NA]
        s_k_elong_left_right=s_k_elong_left_right[!start_with_NA]
        Jk_before=.compute_Jk(v_new_k,s_k,p_k,Y,alpha,w,m,use0=use0,use1=use1)
        c_k_after=floor(unlist(lapply(lapply(v_elong_left_right,.domain,use0),length))*(1-max_gap))
        c_k_after[c_k_after<c]=c
        Jk_after=unlist(mapply(.compute_Jk,v_elong_left_right,s_k_elong_left_right,c_k_after,
                               MoreArgs=list(p_k=p_k,Y=Y,alpha=alpha,w=w,m=m,keep_k=keep_k,use0=use0,use1=use1)))
        best_elong=which.min((Jk_after-Jk_before)/Jk_before)
        if(length(best_elong)>0){
          elongate=((Jk_after-Jk_before)/Jk_before)[best_elong]<deltaJk_elong
        }else{
          elongate=FALSE
        }
        if(elongate){
          return(list(v_new=v_elong_left_right[[best_elong]],
                      v_dom=v_dom_elong_left_right[[best_elong]],
                      s_k=s_k_elong_left_right[[best_elong]]))
        }else{
          return(list(v_new=v_new_k,
                      v_dom=v_dom_k,
                      s_k=s_k))
        }
      },V_new,V_dom,S_k,P_k,len_elong,split(keep,rep(1:K,each=N)),c)
      
      V_new=res_left_right[1,]
      V_dom=res_left_right[2,]
      S_k=res_left_right[3,]
      rm(res_left_right)
      S=matrix(unlist(S_k),ncol=K)
    }
    end=proc.time()
    #message('  elongate: ',round((end-start)[3],2))
    
    ##### find shift warping minimizing dissimilarities ###################################################
    start=proc.time()
    c_k=floor(unlist(lapply(V_new,function(v_new) unlist(lapply(v_new,nrow))[1]))*(1-max_gap))
    c_k[c_k<c]=c
    c_k=rep(c_k,each=length(Y))
    YV=expand.grid(Y,V_new)
    SD=.mapply_custom(cl_probKMA,.find_min_diss,YV[,1],YV[,2],c_k,
                      MoreArgs=list(alpha=alpha,w=w,d=d,use0=use0,use1=use1),SIMPLIFY=TRUE)
    S_new=matrix(SD[1,],ncol=K)
    D_new=matrix(SD[2,],ncol=K)
    browser()
    end=proc.time()
    #message('  find shift: ',round((end-start)[3],2))
    
    ##### compute memberships #############################################################################
    start=proc.time()
    # create membership matrix, with N rows and k columns
    P_new=matrix(0,nrow=N,ncol=K)
    # if dist(y_i,v_k)=0 for some k, set p(i,k)=1 and p(i,h)=0 for h!=k
    D0=apply(D_new,2,'%in%',0)
    for(i in which(rowSums(D0)>1)){
      warning(paste0('Curve ',i,' has dissimilarity 0 from two different motifs. Using only one of them...'))
      D0_select=sample(which(D0[i,]),1)
      D0[i,-D0_select]=FALSE
    }
    D0_index=rowSums(D0)==1
    P_new[D0_index,]=1*D0[D0_index,]
    # if dist(y_i,v_k)>0 for all k
    Dm=as.matrix(D_new[!D0_index,]^(1/(m-1)))
    P_new[!D0_index,]=1/(Dm*rowSums(1/Dm))
    # check degenerate clusters (zero membership)
    for(k in which(colSums(P_new)==0)){
      warning(paste0('Motif ',k,' is degenerate (zero membership). Selecting a new center...'))
      P_new[which.min(D_new[,k]),k]=1
    }
    end=proc.time()
    #message('  compute memberships: ',round((end-start)[3],2))
    
    ##### evaluate objective function #####################################################################
    J_iter[iter]=sum(D_new*(P_new^m),na.rm=TRUE)
    
    ##### compute Bhattacharyya distance between P_old and P_new ##########################################
    BC_dist_k=-log(rowSums(sqrt(P_old*P_new)))
    if(stop_criterion=='max')
      BC_dist=max(BC_dist_k)
    if(stop_criterion=='mean')
      BC_dist=mean(BC_dist_k)
    if(stop_criterion=='quantile')
      BC_dist=quantile(BC_dist_k,prob,type=1)
    BC_dist_iter[iter]=BC_dist    
    
    ##### update ##########################################################################################
    V=V_new
    P=P_new
    S=S_new
    D=D_new
  }
  
  ### prepare output ####################################################################################
  start=proc.time()
  if(iter==iter_max){
    warning('maximum number of iterations reached, method stops')
  }
  # compute motifs
  S_k=split(S,rep(seq_len(K),each=N))
  P_k=split(P,rep(seq_len(K),each=N))
  if(!use0){
    use0=TRUE
    Y=mapply(function(y,y0) list(y0=y0,y1=y$y1),Y,Y0,SIMPLIFY=FALSE)
  }
  V=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
  # compute cleaned motifs
  keep=D<quantile(D,quantile4clean)
  empty_k=which(colSums(keep)==0)
  if(length(empty_k)>0){
    for(k in empty_k)
      keep[which.min(D[,k]),k]=TRUE
  }
  P_clean=P
  P_clean[keep]=1
  P_clean[!keep]=0
  P_k=split(P_clean,rep(seq_len(K),each=N))
  S_clean=S
  V_clean=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
  changed_s=which(unlist(lapply(V_clean,length))>2)
  for(k in changed_s){
    S_clean[,k]=S_clean[,k]+V_clean[[k]]$shift
    V_clean[[k]]=V_clean[[k]][c('v0','v1')]
  }
  S_k=split(S_clean,rep(seq_len(K),each=N))
  V_dom=lapply(V_clean,.domain,use0)
  # compute dissimilarities from cleaned motifs
  D_clean=mapply(function(s_k,v_dom,v_clean,Y){
    Y_in_motifs=mapply(function(y,s){
      if(use0){
        d=ncol(y$y0) # dimension of curves
        y_len=nrow(y$y0)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y0=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y0[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
      }
      if(use1){
        d=ncol(y$y1) # dimension of curves
        y_len=nrow(y$y1)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y1=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y1[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
      }
      y=.select_domain(y,v_dom,use0,use1)
      return(y)
    },Y,s_k,SIMPLIFY=FALSE)
    d=unlist(lapply(Y_in_motifs,.diss_d0_d1_L2,.select_domain(v_clean,v_dom,use0,use1),w,alpha))
    return(d)
  },S_k,V_dom,V_clean,MoreArgs=list(Y))
  output=list(Y0=Y0,Y1=Y1,
              V0=lapply(V,function(v) v$v0),V1=lapply(V,function(v) v$v1),
              V0_clean=lapply(V_clean,function(v) v$v0),V1_clean=lapply(V_clean,function(v) v$v1),
              P=P,P_clean=P_clean,S=S,S_clean=S_clean,
              D=D,D_clean=D_clean,iter=iter,J_iter=J_iter,BC_dist_iter=BC_dist_iter)
  if(return_options){
    output=c(output,list(standardize=standardize,K=K,c=c,c_max=c_max,diss=diss,alpha=alpha,w=w,m=m,
                         iter_max=iter_max,stop_criterion=stop_criterion,quantile=quantile,tol=tol,
                         iter4elong=iter4elong,tol4elong=tol4elong,max_elong=max_elong,trials_elong=trials_elong,deltaJk_elong=deltaJk_elong,max_gap=max_gap,
                         iter4clean=iter4clean,tol4clean=tol4clean))
  }
  if(return_init){
    output=c(output,list(P0=P0,S0=S0))
  }
  end=proc.time()
  #message('output: ',round((end-start)[3],2))
  
  
  ### return output ####################################################################################
  return(output)
}

library(parallel)
z <- probKMA(Y0=Y0,Y1=Y1,standardize=params$standardize,K=params$K,c=params$c,c_max=params$c_max,
             P0=data$P0,S0=data$S0,
             diss=diss,alpha=params$alpha,w=params$w,m=params$m,iter_max=params$iter_max,
             stop_criterion=params$stopCriterion,
             quantile=params$quantile,tol=params$tol,iter4elong=params$iter4elong,
             tol4elong=params$tol4elong,max_elong=params$max_elong,
             trials_elong=params$trials_elong,deltaJk_elong=params$deltaJK_elong,
             max_gap=params$max_gap,params$iter4clean,params$tol4clean,
             params$quantile4clean,params$return_options,TRUE,NULL)
