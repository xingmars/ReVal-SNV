require(dplyr)
require(VGAM)
require(tidyr)

source('VM_model_betabinom.R')
source('VM_funcs.R')

par_df <- read.csv('../results/par.csv',row.names=1, stringsAsFactors = F)

depth.mat <- read.csv('../data/depth.csv',row.names=1,check.names=F) %>%
  data.matrix() %>%
  .[rownames(par_df),]

alt.mat <- read.csv('../data/alt.csv',row.names=1,check.names=F) %>%
  data.matrix() %>%
  .[rownames(par_df),]


pval.mat <- alt.mat
pval.mat[,] <- NA

for (i in 1:nrow(alt.mat)) {
  
  par <- par_df[i,1:2]  %>%
    unlist()
  
  for (j in 1:ncol(alt.mat)) {
    x <- alt.mat[i,j]
    n <- depth.mat[i,j]
    
    pval.mat[i,j] <- 1 - VGAM::pbetabinom(x-1,n,prob=par[1],rho=par[2])
  }
}

write.csv(pval.mat,'../results/pval.csv',quote=F)


pval.mat <- read.csv('../results/pval.csv', stringsAsFactors = F, row.names = 1, check.names=F)

alpha <- 1e-6

call.df <- (pval.mat < alpha) %>%
  as.data.frame(check.names=F) %>%
  mutate(mut_id = rownames(.)) %>%
  gather(-mut_id, key = sample_id, value = call) %>%
  filter(call) %>%
  select(-call)

write.csv(call.df, '../results/positive_calls.csv', row.names=F, quote=F)

