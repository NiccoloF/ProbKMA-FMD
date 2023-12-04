
standardize = TRUE
diss = 'd0_d1_L2' 
alpha = 0.5
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
iter_max=100
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
return_init=TRUE


load("../TempForRcpp/Y_data.Rdata")


params <- list(standardize=standardize, K=K,c = c,c_max = c_max,iter_max = iter_max,
               quantile = quantile,stopCriterion = stop_criterion,tol = tol,
               iter4elong = iter4elong,tol4elong = tol4elong,max_elong = max_elong,
               trials_elong = trials_elong,
               deltaJK_elong = deltaJk_elong,max_gap = max_gap,iter4clean = iter4clean,
               tol4clean = tol4clean,
               quantile4clean = quantile4clean,return_options = return_options,
               return_init = return_init,m = m,w = w,alpha = alpha)

#library(ProbKMAcpp)
Y0_f <- function(Y_i)
{
  return(Y_i$y0)
}
Y1_f <- function(Y_i)
{
  return(Y_i$y1)
}


Y0 <- lapply(Y,Y0_f)
Y1 <- lapply(Y,Y1_f)

a <- ProbKMAcpp::initialChecks(Y0,Y1,P0,S0,params,diss,alpha,w)
params <- a$Parameters
data <- a$FuncData

prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"H1")
dissimilarity <- new(ProbKMAcpp::H1,w,alpha)
motif <-new(ProbKMAcpp::Motif_H1) 
b <- prok$probKMA_run()


l2 = new(L2, c(2,2))
#y <- Y[[2]]
#y <- list(y$y0[1:40],y$y1[1:40])
#v <- Y[[2]]
#v <- list(v$y0[41:80],v$y1[41:80])
#l2$compute(y,v)
# fare test della dichiarazione di ProbKma
V <- list(list(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2),matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)),list(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2),matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)))
mot <- new(Motif_L2,2)
prok = new(ProbKMA,Y,V,params,P0,S0)
prok$set_parameters(params)
prok$run(l2,mot)


######## TEST UP TO ELONGATE MOTIF ####################

prok <- new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0)
dissim <- new(ProbKMAcpp::H1,a$Extra$w,a$Extra$alpha)
mot <-new(ProbKMAcpp::Motif_H1) 

b <- prok$probKMA_run(dissim,mot)


.domain <- function(v,use0){
  if(use0){
    rowSums(!is.na(v[[1]]))!=0
  }else{
    rowSums(!is.na(v[[2]]))!=0
  }
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
    browser()
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

probKMA_prof <- function(Y0=Y0,Y1=Y1,standardize=standardize,K=K,c=c,c_max=c_max,
                         diss=diss,alpha=alpha,w=w,m=m,
                         iter_max=iter_max,stop_criterion=stop_criterion,
                         quantile=quantile,tol=tol,
                         iter4elong=iter4elong,tol4elong=tol4elong,max_elong=max_elong,
                         trials_elong=trials_elong,deltaJk_elong=deltaJk_elong,max_gap=max_gap,
                         iter4clean=iter4clean,tol4clean=tol4clean,
                         quantile4clean=quantile4clean,
                         return_options=TRUE,return_init=TRUE,P0=NULL,S0=NULL){
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
  browser()
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
    browser()
    V_new=mapply(.compute_motif,V_dom,S_k,P_k,MoreArgs=list(Y,m,use0,use1),SIMPLIFY=FALSE)
    changed_s=which(unlist(lapply(V_new,length))>2)
    for(k in changed_s){
      S[,k]=S[,k]+V_new[[k]]$shift
      V_new[[k]]=V_new[[k]][c('v0','v1')]
    }
    S_k=split(S,rep(seq_len(K),each=N))
    V_dom=lapply(V_new,.domain,use0)
    return(list("S_k"=S_k,"P_k"=P_k,"V_dom"=V_dom))
  }
}

use0 <- TRUE
use1 <- TRUE
prof <- probKMA_prof(Y0,Y1,standardize,K,c,c_max,
                     diss,alpha,w,m,
                     iter_max,stop_criterion,
                     quantile,tol,
                     iter4elong,tol4elong,max_elong,
                     trials_elong,deltaJk_elong,max_gap,
                     iter4clean,tol4clean,
                     quantile4clean,
                     return_options=TRUE,return_init=TRUE,P0=NULL,S0=NULL)







