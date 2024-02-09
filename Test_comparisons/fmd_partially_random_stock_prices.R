#This code was used to obtain the candidate motifs
#via functional motif discovery using the partially random initializations

#We do not do the motif search here since we need to generate the candidate
#motifs in order to merge them with the candidate motifs resulting 
#from the random initialization and run the final FMD and search.


source(file ="../Test_comparisons/previous_ProbKMA.R")
#file with all the necessary functions

load("../Test_comparisons/data_wrds.Rdata",smoothed_part2_env <- new.env())



data_smoothed=smoothed_part2_env$data_smoothed



ncurves=length(data_smoothed)

#load("wrds_all_data.Rdata",smoothed_all_stocks<- new.env())
#data_smoothed_all_stocks=smoothed_all_stocks$data_smoothed

#######################
###     RUN FMD     ###
#######################
set.seed(13333)
# use Sobolev-like distance d_0.5
diss = 'd0_d1_L2'
alpha = 0.9

max_gap = 0 # no gaps allowed
iter4elong =1 #5-10 # perform elongation
trials_elong =150 # try all possible elongations 60-70
c_max = 150 # maximum motif length 60-70


### run probKMA multiple times (2x3x10=60 times)

initializations=10
n_init = initializations # number of random initializations to try

V_initt=list()
V_initt[[1]]=list()         #2 clusters
V_initt[[1]][[1]]=list() #2 clusters of length 40
V_initt[[1]][[2]]=list() #2 clusters of length 50
V_initt[[1]][[3]]=list() #2 clusters of length 60


V_initt[[2]]=list()           #3 clusters
V_initt[[2]][[1]]=list()   #3 clusters of length 40
V_initt[[2]][[2]]=list()   #3 clusters of length 50
V_initt[[2]][[3]]=list()   #3 clusters of length 60



for(i in 1:4){
  a40=motifs_init(40,10,40,80,100)
  b40=motifs_init(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a40$motif,ncol=1)),v1=(matrix(a40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b40$motif,ncol=1)),v1=(matrix(b40$motif_derivative,ncol=1))))
  a50=motifs_init(50,10,40,80,100)
  b50=motifs_init(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a50$motif,ncol=1)),v1=(matrix(a50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b50$motif,ncol=1)),v1=(matrix(b50$motif_derivative,ncol=1))))
  a60=motifs_init(60,10,40,80,100)
  b60=motifs_init(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a60$motif,ncol=1)),v1=(matrix(a60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b60$motif,ncol=1)),v1=(matrix(b60$motif_derivative,ncol=1))))
  
  
  aa40=motifs_init(40,10,40,80,100)
  bb40=motifs_init(40,10,40,80,100)
  cc40=motifs_init(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa40$motif,ncol=1),v1=matrix(aa40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb40$motif,ncol=1),v1=matrix(bb40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc40$motif,ncol=1),v1=matrix(cc40$motif_derivative,ncol=1)))
  aa50=motifs_init(50,10,40,80,100)
  bb50=motifs_init(50,10,40,80,100)
  cc50=motifs_init(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa50$motif,ncol=1),v1=matrix(aa50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb50$motif,ncol=1),v1=matrix(bb50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc50$motif,ncol=1),v1=matrix(cc50$motif_derivative,ncol=1)))
  aa60=motifs_init(60,10,40,80,100)
  bb60=motifs_init(60,10,40,80,100)
  cc60=motifs_init(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa60$motif,ncol=1),v1=matrix(aa60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb60$motif,ncol=1),v1=matrix(bb60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc60$motif,ncol=1),v1=matrix(cc60$motif_derivative,ncol=1)))
  
  
}

for(i in 5:8){
  a_2_40=motifs_init_rev(40,10,40,80,100)
  b_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_2_40$motif,ncol=1)),v1=(matrix(a_2_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_40$motif,ncol=1)),v1=(matrix(b_2_40$motif_derivative,ncol=1))))
  a_2_50=motifs_init_rev(50,10,40,80,100)
  b_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_2_50$motif,ncol=1)),v1=(matrix(a_2_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_50$motif,ncol=1)),v1=(matrix(b_2_50$motif_derivative,ncol=1))))
  a_2_60=motifs_init_rev(60,10,40,80,100)
  b_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_2_60$motif,ncol=1)),v1=(matrix(a_2_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_60$motif,ncol=1)),v1=(matrix(b_2_60$motif_derivative,ncol=1))))
  
  
  aa_2_40=motifs_init_rev(40,10,40,80,100)
  bb_2_40=motifs_init_rev(40,10,40,80,100)
  cc_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_2_40$motif,ncol=1),v1=matrix(aa_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_40$motif,ncol=1),v1=matrix(bb_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_40$motif,ncol=1),v1=matrix(cc_2_40$motif_derivative,ncol=1)))
  aa_2_50=motifs_init_rev(50,10,40,80,100)
  bb_2_50=motifs_init_rev(50,10,40,80,100)
  cc_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_2_50$motif,ncol=1),v1=matrix(aa_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_50$motif,ncol=1),v1=matrix(bb_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_50$motif,ncol=1),v1=matrix(cc_2_50$motif_derivative,ncol=1)))
  aa_2_60=motifs_init_rev(60,10,40,80,100)
  bb_2_60=motifs_init_rev(60,10,40,80,100)
  cc_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_2_60$motif,ncol=1),v1=matrix(aa_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_60$motif,ncol=1),v1=matrix(bb_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_60$motif,ncol=1),v1=matrix(cc_2_60$motif_derivative,ncol=1)))
  
  
}

plot(motifs_init(40,10,40,80,100)$motif,type="l",ylab="",xlab="Motif length")
lines(motifs_init(40,10,40,80,100)$motif,col="red")
lines(motifs_init(40,10,40,80,100)$motif,col="green")
lines(motifs_init(40,10,40,80,100)$motif,col="blue")




#Increasing lines
for(i in 9:9){
  a_3_40=motifs_line(40,10,40,80,100)
  b_3_40=motifs_line(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_3_40$motif,ncol=1)),v1=(matrix(a_3_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_40$motif,ncol=1)),v1=(matrix(b_3_40$motif_derivative,ncol=1))))
  a_3_50=motifs_line(50,10,40,80,100)
  b_3_50=motifs_line(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_3_50$motif,ncol=1)),v1=(matrix(a_3_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_50$motif,ncol=1)),v1=(matrix(b_3_50$motif_derivative,ncol=1))))
  a_3_60=motifs_line(60,10,40,80,100)
  b_3_60=motifs_line(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_3_60$motif,ncol=1)),v1=(matrix(a_3_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_60$motif,ncol=1)),v1=(matrix(b_3_60$motif_derivative,ncol=1))))
  
  
  aa_3_40=motifs_line(40,10,40,80,100)
  bb_3_40=motifs_line(40,10,40,80,100)
  cc_3_40=motifs_line(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_3_40$motif,ncol=1),v1=matrix(aa_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_40$motif,ncol=1),v1=matrix(bb_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_40$motif,ncol=1),v1=matrix(cc_3_40$motif_derivative,ncol=1)))
  aa_3_50=motifs_line(50,10,40,80,100)
  bb_3_50=motifs_line(50,10,40,80,100)
  cc_3_50=motifs_line(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_3_50$motif,ncol=1),v1=matrix(aa_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_50$motif,ncol=1),v1=matrix(bb_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_50$motif,ncol=1),v1=matrix(cc_3_50$motif_derivative,ncol=1)))
  aa_3_60=motifs_line(60,10,40,80,100)
  bb_3_60=motifs_line(60,10,40,80,100)
  cc_3_60=motifs_line(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_3_60$motif,ncol=1),v1=matrix(aa_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_60$motif,ncol=1),v1=matrix(bb_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_60$motif,ncol=1),v1=matrix(cc_3_60$motif_derivative,ncol=1)))
  
  
}

#Decreasing lines
for(i in 10:10){
  a_4_40=motifs_line_rev(40,10,40,80,100)
  b_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_4_40$motif,ncol=1)),v1=(matrix(a_4_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_40$motif,ncol=1)),v1=(matrix(b_4_40$motif_derivative,ncol=1))))
  a_4_50=motifs_line_rev(50,10,40,80,100)
  b_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_4_50$motif,ncol=1)),v1=(matrix(a_4_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_50$motif,ncol=1)),v1=(matrix(b_4_50$motif_derivative,ncol=1))))
  a_4_60=motifs_line_rev(60,10,40,80,100)
  b_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_4_60$motif,ncol=1)),v1=(matrix(a_4_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_60$motif,ncol=1)),v1=(matrix(b_4_60$motif_derivative,ncol=1))))
  
  
  aa_4_40=motifs_line_rev(40,10,40,80,100)
  bb_4_40=motifs_line_rev(40,10,40,80,100)
  cc_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_4_40$motif,ncol=1),v1=matrix(aa_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_40$motif,ncol=1),v1=matrix(bb_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_40$motif,ncol=1),v1=matrix(cc_4_40$motif_derivative,ncol=1)))
  aa_4_50=motifs_line_rev(50,10,40,80,100)
  bb_4_50=motifs_line_rev(50,10,40,80,100)
  cc_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_4_50$motif,ncol=1),v1=matrix(aa_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_50$motif,ncol=1),v1=matrix(bb_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_50$motif,ncol=1),v1=matrix(cc_4_50$motif_derivative,ncol=1)))
  aa_4_60=motifs_line_rev(60,10,40,80,100)
  bb_4_60=motifs_line_rev(60,10,40,80,100)
  cc_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_4_60$motif,ncol=1),v1=matrix(aa_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_60$motif,ncol=1),v1=matrix(bb_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_60$motif,ncol=1),v1=matrix(cc_4_60$motif_derivative,ncol=1)))
  
  
}


V_init=V_initt
K = c(2, 3) # number of clusters to try
c=c(40,50)
n_init = 4 # number of partially random initializations to try
# NOTE: rename "results" folder to re-run everything
#(TIME CONSUMING)

# resize V_init
V_init[[1]] <- list(V_init[[1]][[1]],V_init[[1]][[2]])
V_init[[2]] <- list(V_init[[2]][[1]],V_init[[2]][[2]])
V_init[[1]][[1]] <- list(V_init[[1]][[1]][[1]],V_init[[1]][[1]][[2]],V_init[[1]][[1]][[3]],V_init[[1]][[1]][[4]])
V_init[[1]][[2]] <- list(V_init[[1]][[2]][[1]],V_init[[1]][[2]][[2]],V_init[[1]][[2]][[3]],V_init[[1]][[2]][[4]])
V_init[[2]][[1]] <- list(V_init[[2]][[1]][[1]],V_init[[2]][[1]][[2]],V_init[[2]][[1]][[3]],V_init[[2]][[1]][[4]])
V_init[[2]][[2]] <- list(V_init[[2]][[2]][[1]],V_init[[2]][[2]][[2]],V_init[[2]][[2]][[3]],V_init[[2]][[2]][[4]])

#set.seed(13333)
#files = list.files('./motifs_partially_random')
#It is necessary to use the candidate motifs in the folder
#motifs_partially_random in order to obtain the identical results of
#the article.
Y0=Y1=list()
for(i in 1:length(data_smoothed)){
  Y0[[i]]=matrix(data_smoothed[[i]]$v0,ncol=1)
  Y1[[i]]=matrix(data_smoothed[[i]]$v1,ncol=1)#Just the uncorrelated stocks
  
}

#Y0_all=Y1_all=list()
#for(i in 1:length(data_smoothed_all_stocks)){
#  Y0_all[[i]]=matrix(data_smoothed_all_stocks[[i]]$v0,ncol=1)
#  Y1_all[[i]]=matrix(data_smoothed_all_stocks[[i]]$v1,ncol=1)#All the stocks
  
#}

if('motifs_candidate.RData' %in% files){
  # candidate motifs already present, load them
  load('./motifs_partially_random/motifs_candidate.RData')
}else{
  # find candidate motifs
  find_candidate_motifs_results = ProbKMAcpp::find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=V_init,
                                                                    name = '../Test_comparisons/results/our/stock.1', names_var = 'x(t)',
                                                                    probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                                           iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                                           return_options = TRUE, return_init = TRUE,
                                                                                           diss = diss, alpha = alpha,transformed=TRUE),
                                                                    plot = TRUE, exe_print = FALSE, n_threads = 7, set_seed = TRUE, worker_number = NULL)
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=V_init,
                                                        name = '../Test_comparisons/results/prof/stock.1', names_var = 'x(t)',
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha,transformed=TRUE),
                                                        plot = TRUE, worker_number = NULL, set_seed = TRUE)
  save(find_candidate_motifs_results, file = './motifs_partially_random/motifs_candidate.RData')
}

