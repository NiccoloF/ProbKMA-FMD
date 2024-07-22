n_init = 4 # di 10, n_init_motif iterazioni con motif di partenza dato, le restanti non passo il motif
n_init_motif = 2
d = 1 
diss = "d0_d1_L2"
K = c(2, 3) # number of clusters to try
c=c(40,50,60)
#n_init = 4 #
for(i in 1:length(K))
{
  for(j in 1:length(c))
  {
    if(length(V_init[[i]][[j]]) < n_init_motif) # n_init_motif deve dipendere dalla combinazione K e c che fisso ?
    {
      stop('More initial motifs needed')
    }
    else
    {
      if(length(V_init[[i]][[j]]) > n_init_motif)
      {
        # aggiungere un warning per il resize che viene effettuato 
        V_init[[i]][[j]] <- V_init[[i]][[j]][1:n_init_motif]
      }
      for(e in 1:(n_init-n_init_motif))
      {
        V_init[[i]][[j]] <- append(V_init[[i]][[j]], list(NULL))
      }
      for(k in 1:n_init_motif)
      {
        if(length(V_init[[i]][[j]][[k]])!= K[i])
        {
          stop('Number of initial motifs differs from K')
        }
        else
        {
          for(h in 1:K[i])
          {
            if(nrow(V_init[[i]][[j]][[k]][[h]]$v0)!= c[j] || ncol(V_init[[i]][[j]][[k]][[h]]$v0)!= d)
            {
              stop('Uncorrect dimensions of initial motifs provided')
            }
            else
            {
              if(diss == "d0_d1_L2" || diss == "d1_L2")
              {
                if(nrow(V_init[[i]][[j]][[k]][[h]]$v1)!= c[j] || ncol(V_init[[i]][[j]][[k]][[h]]$v0)!= d)
                {
                  stop('Uncorrect dimensions of initial motifs provided')
                }
              }
              V_init_bool = TRUE
            }
          }
        }
      }
    }
  }
}
if(V_init_bool == TRUE){
  print(" V_init it is okay ")
}

if(n_init_motif > 0 && V_init_bool)
{
  i_c_K_motif = expand.grid(seq_len(n_init_motif),c,K)
  i_c_K = expand.grid(seq_len(n_init - n_init_motif),c,K)
  # itero su questi in modo separato, passando in un caso V_init nell'altro no
}
if(n_init_motif == n_init && V_init_bool)
{
  i_c_K_motif = expand.grid(seq_len(n_init_motif),c,K)
  # itero su questi passando solo V_init
}
if(n_init_motif == 0 || !V_init_bool)
{
  i_c_K_motif = expand.grid(seq_len(n_init_motif),c,K)
  # itero su questi senza passare V_init
}
