setwd("C:/Users/buldo/OneDrive/Desktop/progetto pacs/probKMA/ProbKMA-FMD/ProbKMAcpp")
devtools::load_all()
diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0 # no gaps allowed
iter4elong = 1 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 71 # maximum motif length 70
### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(61, 51, 41) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

probKMA_silhouette <- function(Y0, Y1, params, diss, probKMA_results,align=FALSE,plot=TRUE){
  ### compute silhouette #####################################################################################
  if(diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(Y0,function(y0) list(y0=y0,y1=NULL))
  }
  if(diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(Y1,function(y1) list(y0=NULL,y1=y1))
  }
  if(diss=='d0_d1_L2'){
    alpha=params$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),Y0,Y1,SIMPLIFY=FALSE)
  }
  w=params$w
  d=ncol(Y0[[1]])
  N=nrow(probKMA_results$P_clean)
  K=ncol(probKMA_results$P_clean)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))
  
  # extract pieces of curves that fall in the different motifs
  curves_in_motifs=apply(probKMA_results$P_clean,2,function(P_clean_k) which(P_clean_k==1))
  if(!is.null(ncol(curves_in_motifs))){
    curves_in_motifs=split(curves_in_motifs,rep(seq_len(ncol(curves_in_motifs)),each=nrow(curves_in_motifs)))
  }
  curves_in_motifs_number=colSums(probKMA_results$P_clean)
  S_clean_k=mapply(function(s_k,curves_in_motif) s_k[curves_in_motif],
                   split(probKMA_results$S_clean,rep(seq_len(K),each=N)),curves_in_motifs,SIMPLIFY=FALSE)
  
  # compute distances between pieces of curves
  Y_in_motifs=unlist(mapply(function(curves_in_motif,s_k,v_dom){
    Y_in_motif=mapply(function(y,s){
      if(use0){
        d=ncol(y$y0) # dimension of curves
        y_len=nrow(y$y0)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y0=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y0[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
        y$y0[!v_dom,]=NA
      }
      if(use1){
        d=ncol(y$y1) # dimension of curves
        y_len=nrow(y$y1)
        index=max(1,s)-1+seq_len(length(v_dom)-max(0,1-s))
        y$y1=rbind(matrix(NA,nrow=max(0,1-s),ncol=d),
                   matrix(y$y1[index[index<=y_len],],ncol=d),
                   matrix(NA,nrow=sum(index>y_len),ncol=d))
        y$y1[!v_dom,]=NA
      }
      return(y)
    },Y[curves_in_motif],s_k,SIMPLIFY=FALSE)
  },curves_in_motifs,S_clean_k,V_dom,SIMPLIFY=FALSE),recursive=FALSE)
  Y_motifs=rep.int(1:K,curves_in_motifs_number)
  YY=combn(Y_in_motifs,2,simplify=FALSE)
  YY=array(unlist(YY,recursive=FALSE),dim=c(2,length(YY)))
  YY_lengths=combn(V_length[Y_motifs],2)
  swap=YY_lengths[1,]<YY_lengths[2,]
  YY[,swap]=rbind(YY[2,swap],YY[1,swap])
  if(align){
    # find distance between the two pieces of curves
    # no alignment for pieces corresponding to motifs with the same length
    # alignment for pieces corresponding to motifs with different lengths (requiring one piece inside the other)
    equal_length=(YY_lengths[1,]==YY_lengths[2,])
    SD=mapply(.find_diss,YY[1,],YY[2,],equal_length,
              MoreArgs=list(alpha=alpha,w=w,d,use0,use1),SIMPLIFY=TRUE)
  }else{
    # find minimum distance between the two pieces of curves, allowing alignment
    # minimum overlap required: minimum motif length
    cc_motifs=apply(combn(probKMA_results$c[Y_motifs],2),2,min)
    SD=mapply(.find_min_diss,YY[1,],YY[2,],cc_motifs,
              MoreArgs=list(alpha=alpha,w=w,d,use0,use1),SIMPLIFY=TRUE)
  }
  YY_D=matrix(0,nrow=length(Y_motifs),ncol=length(Y_motifs))
  YY_D[lower.tri(YY_D)]=SD[2,]
  YY_D=YY_D+t(YY_D)
  
  # compute intra-cluster distances
  intra=Reduce(cbind,lapply(Y_motifs,function(motif) Y_motifs==motif))
  diag(intra)=FALSE
  a=colSums(intra*YY_D)/(curves_in_motifs_number[Y_motifs]-1)
  
  # compute inter-cluster distances
  b_k=Reduce(rbind,lapply(1:(K-1),
                          function(k){
                            inter=Reduce(cbind,lapply(Y_motifs,
                                                      function(motif) Y_motifs==ifelse(motif+k>K,(motif+k)%%K,motif+k)))
                            b_k=colSums(inter*YY_D)/(curves_in_motifs_number[ifelse(Y_motifs+1>K,(Y_motifs+1)%%K,Y_motifs+1)])
                            return(b_k)}))
  if(is.matrix(b_k)){
    b=apply(b_k,2,min)
  }else{
    b=b_k
  }
  
  # compute silhouette
  silhouette=(b-a)/pmax(a,b)
  silhouette[is.nan(silhouette)]=0
  
  # compute average silhouette per cluster
  silhouette_average=rep(NA,K)
  for(k in seq_len(K)){
    silhouette_k=silhouette[Y_motifs==k]
    curves_in_motifs[[k]]=curves_in_motifs[[k]][order(silhouette_k,decreasing=TRUE)]
    silhouette[Y_motifs==k]=sort(silhouette_k,decreasing=TRUE)
    silhouette_average[k]=mean(silhouette_k)
  }
  
  ### plot silhouette ########################################################################################
  if(plot){
    n=length(silhouette)
    sil=rev(silhouette)
    y=barplot(sil,space=c(0,rev(diff(Y_motifs))),xlab='Silhouette index',names='',
              xlim=c(min(0,min(sil)-0.05),1.2),horiz=TRUE,las=1,mgp=c(2.5,1,0),col='gray')
    text(ifelse(sil>=0,sil+0.03,sil-0.03),ifelse(sil>0,y,y+0.2),labels=rev(unlist(curves_in_motifs)),cex=0.5)
    title(main='Silhouette plot')
    title(sub=paste("Average silhouette width:",round(mean(silhouette_average),2)),adj=0)
    mtext(paste("n =",n),adj=0)
    mtext(substitute(K~~"motifs",list(K=K)),adj=1)
    mtext(expression(paste(j," : ",n[j]," | avg ",s[i])),adj=1.04,line=-1.2)
    y=rev(y)
    for(k in seq_len(K)){
      y_k=mean(y[Y_motifs==k])
      text(1,y_k,
           paste0(k,": ",curves_in_motifs_number[k]," | ",format(silhouette_average[k],digits=1,nsmall=2)),xpd=NA,adj=0.1)
    }
  }
  
  return(list(silhouette=silhouette,motifs=Y_motifs,curves=unlist(curves_in_motifs),
              silhouette_average=silhouette_average))
  
}


probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                       iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                       return_options = TRUE, return_init = TRUE,
                       diss = diss, alpha = alpha)

params = list(standardize=probKMA_options$standardize,
              c_max = probKMA_options$c_max,
              iter_max = probKMA_options$iter_max,
              iter4elong = probKMA_options$iter4elong,
              trials_elong = probKMA_options$trials_elong,
              return_options = probKMA_options$return_options ,
              quantile = 0.25,
              stopCriterion = 'max',
              tol = 1e-8,
              tol4elong = 1e-3,
              max_elong = 0.5,
              deltaJK_elong = 0.05,
              max_gap = 0.2,
              iter4clean = 50,
              tol4clean = 1e-4,
              quantile4clean = 1/2,
              m = 2,
              w = 1,
              alpha = probKMA_options$alpha,
              seed = 1, 
              K = 2,
              c = 40) 

load("../TempForRcpp/len200_sd0.1.RData")


Y0_f <- function(y0)
{
  return(as.matrix(y0))
}
Y1_f <- function(y1)
{
  return(as.matrix(y1))
}

Y0 <- lapply(Y0,Y0_f)
Y1 <- lapply(Y1,Y1_f)


library(ProbKMAcpp)

a <- ProbKMAcpp::initialChecks(Y0,Y1,matrix(),matrix(),params,diss,1)
params <- a$Parameters
data <- a$FuncData

prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"H1")

c = 61
K = 2
devtools::load_all()
params$c = c
params$K = K
params$quantile4clean = 1/K
# add an option for the seed, maybe P0 and S0 have to change
params$c_max = probKMA_options$c_max
checked_data <- ProbKMAcpp::initialChecks(Y0,Y1,matrix(),matrix(),params,diss,1)
params <- checked_data$Parameters
data <- checked_data$FuncData
prok$reinit_motifs(params$c,ncol(as.matrix(Y0[[1]])))
initial_motifs <-prok$get_motifs()
prok$set_P0(data$P0)
prok$set_S0(data$S0)
prok$set_parameters(params)
#params <- prok$get_parameters()
silhouette_align = FALSE

devtools::load_all()
i_c_K=expand.grid(seq_len(n_init),c,K)
results=mapply(function(K,c,i,params,prok,data,Y0,Y1,diss,probKMA_options,silhouette_align){ 
  iter=iter_max=1
  while(iter==iter_max){
    start=proc.time()
    # the next line commented is the old version
    # probKMA_results=do.call(probKMA,c(list(Y0=Y0,Y1=Y1,K=K,c=c),probKMA_options))
    params$c = c
    params$K = K
    params$w = 1
    params$quantile4clean = 1/K
    # add an option for the seed, maybe P0 and S0 have to change
    params$c_max = probKMA_options$c_max
    checked_data <- ProbKMAcpp::initialChecks(Y0,Y1,matrix(),matrix(),params,diss,1)
    params <- checked_data$Parameters
    data <- checked_data$FuncData
    prok$reinit_motifs(params$c,ncol(as.matrix(Y0[[1]])))
    prok$set_P0(data$P0)
    prok$set_S0(data$S0)
    prok$set_parameters(params)
    probKMA_results = prok$probKMA_run() # new run for probKMA with updated parameters
    end=proc.time()
    time=end-start
    iter=probKMA_results$iter
    iter_max=params$iter_max
    if(iter==iter_max)
      warning('Maximum number of iteration reached. Re-starting.')
  }
  #pdf(paste0(name,"_K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
  #probKMA_plot(Y0, Y1, probKMA_results,ylab=names_var,cleaned=FALSE) 
  #dev.off()
  #pdf(paste0(name,"_K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
  #probKMA_plot(Y0, Y1, probKMA_results,ylab=names_var,cleaned=TRUE) 
  #dev.off()
  #pdf(paste0(name,"_K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
  silhouette=probKMA_silhouette(Y0,
                                Y1,
                                params,
                                diss,
                                probKMA_results,
                                align=silhouette_align,
                                plot=TRUE) 
  dev.off()
  #save(probKMA_results,time,silhouette,
  #     file=paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
  return(list(probKMA_results=probKMA_results,
              time=time,silhouette=silhouette))
},i_c_K[,3],i_c_K[,2],i_c_K[,1],SIMPLIFY=FALSE,MoreArgs = list(params,prok,data,Y0,Y1,diss,probKMA_options,silhouette_align))
