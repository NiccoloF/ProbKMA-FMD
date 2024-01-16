library(ProbKMAcpp)

mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}

diss_d0_d1_L2 <- function(y,v,w,alpha){
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

domain <- function(v,use0){
  if(use0){
    rowSums(!is.na(v[[1]]))!=0
  }else{
    rowSums(!is.na(v[[2]]))!=0
  }
}

select_domain <- function(v,v_dom,use0,use1){
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}

find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1){
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


find_diss <- function(y,v,alpha,w,aligned,d,use0,use1){
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


compute_motif <- function(v_dom,s_k,p_k,Y,m,use0,use1){
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


compute_Jk <- function(v,s_k,p_k,Y,alpha,w,m,use0,use1,c_k=NULL,keep_k=NULL){
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

probKMA_plot <- function(Y0,Y1,probKMA_results,ylab='',cleaned=FALSE){
  # my idea was passing data instead of Y0 and Y1 but in this way when data suppresses Y0 or Y1 it does not work
  d=ncol(Y0[[1]])
  N=nrow(probKMA_results$P0)
  K=ncol(probKMA_results$P0)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  S_k=split(probKMA_results$S0,rep(seq_len(K),each=N))
  P_k=split(probKMA_results$P0,rep(seq_len(K),each=N)) 
  
  ### plot motifs with matched curves ########################################################################
  if(cleaned){
    S_clean_k=split(probKMA_results$S_clean,rep(seq_len(K),each=N))
    P_clean_k=split(probKMA_results$P_clean,rep(seq_len(K),each=N))
    if(is.null(probKMA_results$V1[[1]])){
      mapply(function(v,v_dom,s_k,p_clean_k,k){
        keep=which(p_clean_k==1)
        Y_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }else{
      mapply(function(v0,v1,v_dom,s_k,p_clean_k,k){
        keep=which(p_clean_k==1)
        Y0_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          Y1[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v0[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y1_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0_clean,probKMA_results$V1_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }
  }else{
    if(is.null(probKMA_results$V1[[1]])){
      mapply(function(v,v_dom,s_k,p_k,k){
        Y_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,V_dom,S_k,P_k,seq_len(K))
    }else{
      mapply(function(v0,v1,v_dom,s_k,p_k,k){
        Y0_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          Y1,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v0[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y1_inters_k))
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,probKMA_results$V1,V_dom,S_k,P_k,seq_len(K))
    }
  }
  layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
  
  ### plot motifs ############################################################################################
  if(cleaned){
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0_clean,nrow))
             plot(probKMA_results$V0_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0_clean),na.rm=TRUE),max(unlist(probKMA_results$V0_clean),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0_clean[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }else{
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0,nrow))
             plot(probKMA_results$V0[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0),na.rm=TRUE),max(unlist(probKMA_results$V0),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }
  if(!is.null(probKMA_results$V1[[1]])){
    layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
    if(cleaned){
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1_clean,nrow))
               plot(probKMA_results$V1_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1_clean),na.rm=TRUE),max(unlist(probKMA_results$V1_clean),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1_clean[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }else{
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1,nrow))
               plot(probKMA_results$V1[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1),na.rm=TRUE),max(unlist(probKMA_results$V1),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }
  }
  
  ### plot memberships #######################################################################################
  par(mfrow=c(K,1),mar=c(3,4,4,2)+0.1)
  if(cleaned){
    mapply(function(p_k,p_clean_k,k){
      col=rep('lightgray',N)
      col[p_clean_k==1]='gray35'
      barplot(p_k,names.arg=seq_len(N),col=col,las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif',k))
    },P_k,P_clean_k,seq_len(K))
  }else{
    mapply(function(p_k,k){
      barplot(p_k,names.arg=seq_len(N),col='gray',las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif',k))
    },P_k,seq_len(K))
  }
  
  ### plot objective function and Bhattacharyya distance #####################################################
  par(mfrow=c(1,1))
  plot(seq_len(probKMA_results$iter),probKMA_results$J_iter,type='l',xlab='iter',ylab='objective function',main='Objective function Jm')
  plot(seq_len(probKMA_results$iter),probKMA_results$BC_dist_iter,type='l',xlab='iter',ylab='distance between memberships',main='Bhattacharyya distance between memberships')
  
  return()
}

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

find_candidate_motifs <- function(Y0,Y1=NULL,K,c,n_init=10,name='results',names_var='',
                                  probKMA_options=NULL,silhouette_align=FALSE,plot=TRUE,worker_number=NULL){
  ### check input #############################################################################################
  # check required input
  if(missing(K))
    stop('K must be specified')
  if(missing(c))
    stop('c must be specified')
  # check K
  if(!is.vector(K))
    stop('K should be a vector.')
  # check c
  if(!is.vector(c))
    stop('c should be a vector.')
  # check n_init
  if(n_init%%1!=0)
    stop('Number of initializations n_init should be an integer.')
  if(n_init<1)
    stop('Number of initializations n_init should be at least 1.')
  # check name
  if(!is.character(name))
    stop('name should be a string.')
  if(grepl(' ',name)){
    warning('name contains spaces. Changing spacing in underscores.')
    name=gsub(" ","_",name,fixed=TRUE)
  }
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
  
  # set return_options=TRUE
  if(!is.null(probKMA_options$return_options)){
    if(probKMA_options$return_options==FALSE){
      warning('Setting return_option=TRUE')
      probKMA_options$return_options=TRUE
    }
  }
  
  params = list(standardize=probKMA_options$standardize,
                c_max = probKMA_options$c_max,
                iter_max = probKMA_options$iter_max,
                iter4elong = probKMA_options$iter4elong,
                trials_elong = probKMA_options$trials_elong,
                return_options = probKMA_options$return_options ,
                quantile = 0.25, #argument
                stopCriterion = 'max', #argument
                tol = 1e-8, #argument
                tol4elong = 1e-3, #argument
                max_elong = 0.5, #argument
                deltaJK_elong = 0.05, #argument
                max_gap = 0.2, #argument
                iter4clean = 50, #argument
                tol4clean = 1e-4, #argument
                quantile4clean = 1/2, #argument
                m = 2, #argument
                w = 1, #argument
                alpha = probKMA_options$alpha,#argument
                seed = 1, #argument
                K = 2, #argument
                c = 40) #argument
  
  checked_data <- ProbKMAcpp::initialChecks(Y0,Y1,
                                            matrix(),matrix(),
                                            params,
                                            probKMA_options$diss,
                                            1)
  
  params <- checked_data$Parameters
  data <- checked_data$FuncData
  
  if (alpha == 0 || alpha == 1){
    prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"L2")
  } else {
    prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"H1")
  }
  
  ### set parallel jobs #############################################################################
  core_number <- parallel::detectCores()
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
  len_mean=mean(unlist(lapply(Y0,nrow)))
  c_mean=mean(c)
  if(((len_mean/c_mean>10)&(length(K)*length(c)*n_init<40))|(length(K)*length(c)*n_init==1)){
    ## parallelism inside probKMA is managed by C++ and not by R, the next 2 lines not necessary
    ## if(is.null(probKMA_options$worker_number))
    ##  probKMA_options$worker_number=worker_number
    cl_find=NULL
  }else{
    probKMA_options$worker_number=1
    if(worker_number>1){
      cl_find=parallel::makeCluster(worker_number,timeout=60*60*24*30)
      parallel::clusterEvalQ(cl_find, library(ProbKMAcpp))
      parallel::clusterExport(cl_find,c('name','names_var','Y0','Y1','probKMA_options',
                                        'probKMA_plot','probKMA_silhouette','compute_motif', 
                                        'mapply_custom','diss_d0_d1_L2','domain',
                                        'select_domain','find_min_diss','compute_Jk',
                                        'params','prok','data'),envir=environment()) 
      parallel::clusterCall(cl_find,function()library(parallel,combinat)) # here is the problem
      on.exit(parallel::stopCluster(cl_find))
    }else{
      cl_find=NULL
    }
  }
  
  browser()
  
  ### run probKMA ##########################################################################################
  i_c_K=expand.grid(seq_len(n_init),c,K)
  results=mapply_custom(NULL,function(K,c,i){ #cl_find
    #dir.create(paste0(name,"_K",K,"_c",c),showWarnings=FALSE)
    #files=list.files(paste0(name,"_K",K,"_c",c))
    #message("K",K,"_c",c,'_random',i)
    #if(paste0('random',i,'.RData') %in% files){
    #  load(paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
    #  return(list(probKMA_results=probKMA_results,
    #              time=time,silhouette=silhouette))
    #}else{
      iter=iter_max=1
      while(iter==iter_max){
        start=proc.time()
        # the next line commented is the old version
        # probKMA_results=do.call(probKMA,c(list(Y0=Y0,Y1=Y1,K=K,c=c),probKMA_options))
        params$c = c
        params$K = K
        params$quantile4clean = 1/K
        # add an option for the seed, maybe P0 and S0 have to change
        params$c_max = probKMA_options$c_max
        checked_data <- ProbKMAcpp::initialChecks(Y0,Y1,matrix(),matrix(),params,probKMA_options$diss,1)
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
      probKMA_plot(Y0, Y1, probKMA_results,ylab=names_var,cleaned=FALSE) 
      dev.off()
      #pdf(paste0(name,"_K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
      probKMA_plot(Y0, Y1, probKMA_results,ylab=names_var,cleaned=TRUE) 
      dev.off()
      #pdf(paste0(name,"_K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
      silhouette=probKMA_silhouette(Y0,
                                    Y1,
                                    params,
                                    probKMA_options$diss,
                                    probKMA_results,
                                    align=silhouette_align,
                                    plot=TRUE) 
      dev.off()
      save(probKMA_results,time,silhouette,
           file=paste0(name,"_K",K,"_c",c,'/random',i,'.RData'))
      return(list(probKMA_results=probKMA_results,
                  time=time,silhouette=silhouette))
    }
  #}
  ,i_c_K[,3],i_c_K[,2],i_c_K[,1],SIMPLIFY=FALSE)
  
  browser() 
  
  results=split(results,list(factor(i_c_K[,2],c),factor(i_c_K[,3],K)))
  results=split(results,rep(K,each=length(c)))
  
  ### plot silhouette average #################################################################################
  silhouette_average_sd=lapply(results,
                               function(results){
                                 silhouette_average=lapply(results,
                                                           function(results){
                                                             silhouette_average=numeric(n_init)
                                                             silhouette_sd=numeric(n_init)
                                                             for(i in seq_len(n_init)){
                                                               silhouette_average[i]=mean(results[[i]]$silhouette$silhouette_average)
                                                               silhouette_sd[i]=sd(results[[i]]$silhouette$silhouette_average)
                                                             }
                                                             return(cbind(silhouette_average,silhouette_sd))
                                                           })
                               })
  if(plot){
    pdf(paste0(name,'_silhouette.pdf'),width=7,height=5)
    for(i in seq_along(K)){
      silhouette_average_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),1])),ncol=length(c))
      silhouette_sd_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),2])),ncol=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
              silhouette_average_plot,type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,1),xaxt='n',
              xlab='',ylab='Silhouette average',main=paste0('K=',K[i]))
      shift=seq(-0.1,0.1,length.out=length(c))
      for(ii in seq_along(c)){
        segments(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                 1:n_init+shift[ii],silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
        segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],
                 1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
        segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                 1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],col=ii+1)
        text(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
             silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
    }
    dev.off()
  }
  
  ### plot processing time ####################################################################################
  times=lapply(results,
               function(results){
                 times=lapply(results,
                              function(results)
                                times=unlist(lapply(results,function(results) results$time[3])))
               })
  if(plot){
    pdf(paste0(name,'_times.pdf'),width=7,height=5)
    y_max=max(unlist(times))*1.2
    times_plot=vector('list',length(K))
    for(i in seq_along(K)){
      times_plot[[i]]=Reduce(cbind,
                             lapply(seq_along(c),
                                    function(j)
                                      as.matrix(times[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]))
      )
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
              times_plot[[i]],type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,y_max),xaxt='n',
              xlab='',ylab='Time',main=paste0('K=',K[i]))
      shift=seq(-0.1,0.1,length.out=length(c))
      for(ii in seq_along(c)){
        text(1:n_init+shift[ii],times_plot[[i]][,ii],silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
      }
      legend('topleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
    }
    for(i in seq_along(K)){
      boxplot(times_plot[[i]],col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Times',main=paste0('K=',K[i]))
    }
    dev.off()
  }
  
  ### plot dissimilarities ####################################################################################
  if(plot){
    D=lapply(results,
             function(results){
               D=lapply(results,
                        function(results){
                          D=lapply(results,
                                   function(results){
                                     D=as.vector(results$probKMA_results$D)
                                   })
                        })
             })
    pdf(paste0(name,'_dissimilarities.pdf'),width=7,height=5)
    y_max=max(unlist(D))
    for(i in seq_along(K)){
      D_plot=matrix(unlist(D[[i]]),ncol=length(c))
      boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
    }
    for(i in seq_along(K)){
      for(j in seq_along(c)){
        silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
        D_plot=D[[i]][[j]][silhouette_order]
        boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
    }
    dev.off()
    
    D_clean=lapply(results,
                   function(results){
                     D_clean=lapply(results,
                                    function(results){
                                      D_clean=lapply(results,
                                                     function(results){
                                                       D_clean=as.vector(results$probKMA_results$D_clean)
                                                     })
                                    })
                   })
    pdf(paste0(name,'_dissimilarities_clean.pdf'),width=7,height=5)
    y_max=max(unlist(D_clean))
    for(i in seq_along(K)){
      D_plot=matrix(unlist(D_clean[[i]]),ncol=length(c))
      boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
              xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
    }
    for(i in seq_along(K)){
      for(j in seq_along(c)){
        silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
        D_plot=D_clean[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]
        boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
    }
    dev.off()
  }
  
  ### plot motif lengths ######################################################################################
  if(plot){
    motif_length=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0,nrow))
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    pdf(paste0(name,'_lengths.pdf'),width=7,height=5)
    motif_length_plot=lapply(motif_length,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
      abline(h=c,col=1+seq_along(c))
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    motif_clean_length=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0_clean,nrow))
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    pdf(paste0(name,'_lengths_clean.pdf'),width=7,height=5)
    motif_length_plot=lapply(motif_clean_length,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
      abline(h=c,col=1+seq_along(c))
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    pdf(paste0(name,'_lengths_perc.pdf'),width=7,height=5)
    motif_length_perc=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0,nrow))/results$probKMA_results$c*100
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    motif_length_plot=lapply(motif_length_perc,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
      abline(h=100)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
    
    pdf(paste0(name,'_lengths_clean_perc.pdf'),width=7,height=5)
    motif_clean_length_perc=mapply(function(results){
      motif_length=mapply(function(results){
        motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                   function(results){
                                                     unlist(lapply(results$probKMA_results$V0_clean,nrow))/results$probKMA_results$c*100
                                                   })))
        return(as.matrix(motif_length))
      },results,SIMPLIFY=FALSE)
    },results,SIMPLIFY=FALSE)
    motif_length_plot=lapply(motif_clean_length_perc,
                             function(motif_length){
                               lapply(motif_length,
                                      function(motif_length){
                                        motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=10^-1),ncol=n_init)
                                      })
                             })
    ymax=max(unlist(motif_length_plot))
    for(i in seq_along(K)){
      plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
      abline(h=100)
      axis(1,1:n_init,labels=rep('',n_init))
      mtext("Init",side=1,line=1)
      shift=seq(-0.1,0.1,length.out=length(c))
      silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
      for(ii in seq_along(c)){
        points(rep(1:n_init+shift[ii],each=K[i]),
               motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
        text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
             rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
      }
      legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
    }
    dev.off()
  }
  
  ### output ##################################################################################################
  return(list(name=name,K=K,c=c,n_init=n_init,silhouette_average_sd=silhouette_average_sd,times=times))
}

#############################
### SIMULATION SCENARIO 1 ###
#############################
setwd("C:/Users/buldo/OneDrive/Desktop/progetto pacs/probKMA/ProbKMA-FMD/ProbKMAcpp")
devtools::load_all()

load("../TempForRcpp/len200_sd0.1.RData")

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


find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,
                                                      name = './results/len200_sd0.1', names_var = 'x(t)',
                                                      probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                             iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                             return_options = TRUE, return_init = TRUE,
                                                                             diss = diss, alpha = alpha),
                                                      plot = TRUE, worker_number = NULL)
