# set the directory of the package, otherwise initialChecks could give problems
setwd("../ProbKMAcpp")
devtools::load_all()

require(fda)

plot_P <- function(probKMA_results,labels,col,names.arg){
  d=ncol(probKMA_results$Y0[[1]])
  N=nrow(probKMA_results$P)
  K=ncol(probKMA_results$P)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  S_k=split(probKMA_results$S,rep(seq_len(K),each=N))
  P_k=split(probKMA_results$P,rep(seq_len(K),each=N))
  
  ### plot memberships #######################################################################################
  par(mfrow=c(K,1),mar=c(9,4,3,2)+0.1)
  mapply(function(p_k,k){
    col_k=rep('lightgray',N)
    col_k[labels==k]=col[k]
    barplot(p_k,names.arg=names.arg,col=col_k,las=2,ylim=c(0,1),ylab='Memberships',main=paste('Cluster',k))
  },P_k,seq_len(K))
  
  par(mfrow=c(K,1),mar=c(9,4,3,2)+0.1)
  mapply(function(p_k,k){
    col_k=rep('lightgray',N)
    col_k[probKMA_results$P_clean[,k]==1]=col[k]
    barplot(p_k,names.arg=names.arg,col=col_k,las=2,ylim=c(0,1),ylab='Memberships',main=paste('Cluster',k))
  },P_k,seq_len(K))
  
  return()
}


##############################
##### ISTAT SURPLUS DATA #####
##############################

load('../Example/istat_mortality_rates_smoothed.Rdata')

### change scale
for(i in 1:3){
  istat_decessi_reg_smooth$fd[[i]]$coefs = istat_decessi_reg_smooth$fd[[i]]$coefs * 100000
  istat_decessi_reg_smooth$mat[[i]] = istat_decessi_reg_smooth$mat[[i]] * 100000
  istat_decessi_reg_smooth$mat_deriv[[i]] = istat_decessi_reg_smooth$mat_deriv[[i]] * 100000
}
ncurves=ncol(istat_decessi_reg_smooth$mat$deceduti_diff)


#### prob k-mean with local alignment, L2 distance between curves (d0), K=2, c=65 ####
c_min <- 65 # max 10 days of shift between each pair of curves

## K=2 (try 10 different initializations, select the best clustering)
seed <- 13
set.seed(seed)
K=2

# additional parameters 
diss = 'd0_L2'
P0 = matrix()
S0 = matrix()
params <- list(standardize=FALSE, K=2,c = c_min,c_max = c_min,iter_max = 1000,
               quantile = 0.25,stopCriterion = 'max',tol = 1e-8,
               iter4elong = 2,tol4elong = 1e-3,max_elong = 0.5,
               trials_elong = 10,
               deltaJK_elong = 0.05,max_gap = 0.2,iter4clean = 1000,
               tol4clean = 1e-4,
               quantile4clean = 1/2,return_options = TRUE,
               m = 2,w = 1, alpha = 0.0,seed = seed,exe_print = TRUE, 
               set_seed= FALSE, n_threads = 7) 

to_matrix <- function(y0){
  return(as.matrix(y0))
}

Y0 = lapply(lapply(1:ncurves, function(i) istat_decessi_reg_smooth$mat$deceduti_diff[,i]),to_matrix)

# checks the parameters
a <- ProbKMAcpp::initialChecks(Y0,NULL,P0,S0,params,diss,seed)
params <- a$Parameters
data <- a$FuncData

prok = new(ProbKMAcpp::ProbKMA,data$Y,data$V,params,data$P0,data$S0,"L2")

J_K2_d0_c65 <- rep(NA, 10)
probKMA_K2_d0_c65_all <- vector('list', 10)

for(i in 1:10){ 
  probKMA_K2_d0_c65_all[[i]] <- prok$probKMA_run()
  J_K2_d0_c65[i] <- probKMA_K2_d0_c65_all[[i]]$J_iter[probKMA_K2_d0_c65_all[[i]]$iter]
  a <- ProbKMAcpp::initialChecks(Y0,NULL,matrix(),matrix(),params,diss,seed)
  data <- a$FuncData
  prok$reinit_motifs(params$c,ncol(Y0[[1]]))
  prok$set_P0(data$P0)
  prok$set_S0(data$S0)
}
save(J_K2_d0_c65, probKMA_K2_d0_c65_all, file="../Example/results/probKMA_K2_d0_c65_results.RData")

load("../Example/results/probKMA_K2_d0_c65_results.RData")
probKMA_K2_d0_c65 <- probKMA_K2_d0_c65_all[[which.min(J_K2_d0_c65)]]
probKMA_K2_d0_c65_labels <- colSums(t(probKMA_K2_d0_c65$P > 0.5) * 1:2) # Dicotomize membership based on probability>0.5
probKMA_K2_d0_c65_P <- round(probKMA_K2_d0_c65$P[, 1], 2)
col_K2_d0_c65_P <- rgb(colorRamp(c("red", "blue"))( probKMA_K2_d0_c65_P ), maxColorValue = 255)

### barplot of probabilitic memberships and distance from mean
pdf('../Example/results/probKMA_K2_d0_c65_memberships.pdf', width = 7, height = 5.5)
plot_P(probKMA_K2_d0_c65, labels=probKMA_K2_d0_c65_labels, col=c("blue", "red"), names.arg = colnames(istat_decessi_reg_smooth$mat$deceduti_diff))
dev.off()

### silhouette plot
pdf('../Example/results/probKMA_K2_d0_c65_silhouette.pdf', width = 7, height = 8.5)
par(mfrow=c(1,1))
ProbKMAcpp::probKMA_silhouette(Y0, NULL, params, diss, probKMA_K2_d0_c65,plot=TRUE)

### silhouette plot after dicotomizing membership based on probability>0.5
probKMA_K2_d0_c65$P_clean <- (probKMA_K2_d0_c65$P > 0.5) * 1
ProbKMAcpp::probKMA_silhouette(Y0, NULL, params, diss, probKMA_K2_d0_c65,plot=TRUE)
dev.off()

### distances
pdf('../Example/results/probKMA_K2_d0_c65_distances.pdf', width = 7, height = 5.5)
hist(probKMA_K2_d0_c65$D, breaks = 30, xlab = "Distances portions to cluster centers", main = "probKMA d0 c=65 days", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
abline(v = median(probKMA_K2_d0_c65$D), col = "red", lty = 2, lwd = 2)

boxplot(list("Cluster 1"=probKMA_K2_d0_c65$D[probKMA_K2_d0_c65_labels == 1, 1], 
             "Cluster 2"=probKMA_K2_d0_c65$D[probKMA_K2_d0_c65_labels == 2, 2]), 
        col = c("blue", "red"), ylab = 'Distances portions to cluster centers', main = "probKMA d0 c=65 days", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
text(2, probKMA_K2_d0_c65$D[9, 2], labels = "Lombardia", pos = 4, cex = 1.5)
dev.off()

### plot clusters
probKMA_K2_d0_c65_shift <- rep(NA, length(probKMA_K2_d0_c65_labels))
probKMA_K2_d0_c65_shift[probKMA_K2_d0_c65_labels == 1] <- probKMA_K2_d0_c65$S[probKMA_K2_d0_c65_labels == 1, 1]
probKMA_K2_d0_c65_shift[probKMA_K2_d0_c65_labels == 2] <- probKMA_K2_d0_c65$S[probKMA_K2_d0_c65_labels == 2, 2]
probKMA_K2_d0_c65_lengths <- unlist(lapply(probKMA_K2_d0_c65$V0, nrow))

pdf('../Example/results/probKMA_K2_d0_c65_results.pdf', width = 7, height = 5.5)
par(mar = c(5, 5, 4, 4) + 0.1)
x=seq_len(nrow(istat_decessi_reg_smooth$mat$deceduti_diff))
matplot(x,istat_decessi_reg_smooth$mat$deceduti_diff, type = 'l', lty = 1, col = 'darkgray', 
        ylab = expression("Excess mortality rate ("%*%"10"^-5*")"), xlab='', xaxt='n', main = 'probKMA d0 c=65 days', 
        lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
abline(h = 0, lty=2, col='darkgray')
abline(v = c(14,14+31)+0.5, lty=2, col='darkgray')
abline(v = 14+9-0.5) # March 9 = national lock down
abline(v = 14+23-0.5) # Mar 23 = closing of all non-essential economic activities
#abline(v = 14+31+30+4-0.5) # May 4 = partial reopening
axis(1, at = x, labels = c(paste(16:29,'Feb'), paste(1:31,'Mar'), paste(1:30,'Apr')), las=2, cex.axis = 1.5)
for(i in seq_along(probKMA_K2_d0_c65_labels)){
  index <- probKMA_K2_d0_c65_shift[i] + seq_len(probKMA_K2_d0_c65_lengths[probKMA_K2_d0_c65_labels[i]]) - 1
  index <- index[(index > 0) & (index <= length(x))]
  lines(x[index], istat_decessi_reg_smooth$mat$deceduti_diff[index, i], col = col_K2_d0_c65_P[i], lty = i, lwd = 2)
}
legend('topright', legend = c("Cluster 1", "Cluster 2"), col = c("blue", "red"), lty = 1, cex = 1.5)

### Plot of cluster 1 vs cluster 2 aligned
i_cluster1 <- which(probKMA_K2_d0_c65_labels == 1)
i_cluster2 <- which(probKMA_K2_d0_c65_labels == 2)
Y0_cluster1 <- c()
for(i in i_cluster1)
  Y0_cluster1 <- cbind(Y0_cluster1, istat_decessi_reg_smooth$mat$deceduti_diff[ , i][probKMA_K2_d0_c65_shift[i] + seq_len(probKMA_K2_d0_c65_lengths[1]) - 1])
Y0_cluster1_mean <- rowMeans(Y0_cluster1, na.rm = TRUE)
Y0_cluster2 <- c()
for(i in i_cluster2)
  Y0_cluster2 <- cbind(Y0_cluster2, istat_decessi_reg_smooth$mat$deceduti_diff[ , i][probKMA_K2_d0_c65_shift[i] + seq_len(probKMA_K2_d0_c65_lengths[2]) - 1])
Y0_cluster2_mean <- rowMeans(Y0_cluster2, na.rm = TRUE)

par(mar = c(6.5, 5, 4, 4) + 0.1)
matplot(x[median(probKMA_K2_d0_c65_shift[i_cluster1]) - 1 + seq_len(probKMA_K2_d0_c65_lengths[1])], 
        Y0_cluster1, type = 'l', col = col_K2_d0_c65_P[i_cluster1], xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        ylab = expression("Excess mortality rate ("%*%"10"^-5*")"), xlab='', xaxt='n', main = 'probKMA d0 c=65 days - Cluster 1',
        lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_K2_d0_c65_shift[i_cluster1]) - 1 + seq_len(probKMA_K2_d0_c65_lengths[1])], 
      Y0_cluster1_mean, col = "black", lwd = 3)
matplot(x[median(probKMA_K2_d0_c65_shift[i_cluster2]) - 1 + seq_len(probKMA_K2_d0_c65_lengths[2])], 
        Y0_cluster2, type = 'l', col = col_K2_d0_c65_P[i_cluster2], xlim = range(x), ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        ylab = expression("Excess mortality rate ("%*%"10"^-5*")"), xlab='', xaxt='n', main = 'probKMA d0 c=65 days - Cluster 2',
        lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(x[median(probKMA_K2_d0_c65_shift[i_cluster2]) - 1 + seq_len(probKMA_K2_d0_c65_lengths[2])], 
      Y0_cluster2_mean, col = "black", lwd = 3)
# together
matplot(seq_len(probKMA_K2_d0_c65_lengths[1]), 
        Y0_cluster1, type = 'l', col = col_K2_d0_c65_P[i_cluster1], ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        ylab = expression("Excess mortality rate ("%*%"10"^-5*")"), xlab='', xaxt='n', main = 'probKMA d0 c=65 days',
        lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
lines(seq_len(probKMA_K2_d0_c65_lengths[1]), 
      Y0_cluster1_mean, col = "black", lwd = 3)
matplot(seq_len(probKMA_K2_d0_c65_lengths[2]), 
        Y0_cluster2, type = 'l', col = col_K2_d0_c65_P[i_cluster2], ylim = range(cbind(Y0_cluster1, Y0_cluster2)),
        ylab = expression("Excess mortality rate ("%*%"10"^-5*")"), xlab='', xaxt='n', main = 'probKMA d0 c=65 days',
        lwd = 2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, add = TRUE)
lines(seq_len(probKMA_K2_d0_c65_lengths[2]), 
      Y0_cluster2_mean, col = "black", lwd = 3)
legend("topleft", legend = paste('Cluster', 1:2), col = c('blue', 'red'), lty = 1, lwd = 2, cex = 1.5, bty = "n")
dev.off()


### barplot of the shifts
# how to interpret: where the part of curve aligned starts
# e.g: trento/bolzano starts 8 days after Lombardia
pdf('../Example/results/probKMA_K2_d0_c65_alignment.pdf', width = 7, height = 5.5)
par(mar = c(6.5, 5, 4, 4) + 0.1)
order=order(probKMA_K2_d0_c65_labels)
barplot(probKMA_K2_d0_c65_shift[order], ylab = 'Start', main = 'Alignment',
        col = c("blue", "red")[probKMA_K2_d0_c65_labels[order]], beside = F, horiz = F, 
        cex.names = 1.2, cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
#rotate 40 degrees (srt = 60)
text(seq(0.6, 0.6 + 1.2*length(probKMA_K2_d0_c65_shift) - 1, by = 1.2), par("usr")[3]-0.25, 
     srt = 40, adj = 1, xpd = TRUE,
     labels = colnames(istat_decessi_reg_smooth$mat$deceduti_diff)[order], cex = 1.2)
dev.off()

#####


