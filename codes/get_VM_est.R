source('VM_model_betabinom.R')
source('VM_funcs.R')

require(dplyr)
require(parallel)
require(doParallel)

if (!dir.exists('../results/par')) dir.create('../results/par')
  
saved_ids <- list.files('../results/par') %>%
  gsub('\\.csv$','',.)

depth.mat <- read.csv('../data/depth.csv',row.names=1,check.names=F) %>%
  data.matrix()

alt.mat <- read.csv('../data/alt.csv',row.names=1,check.names=F) %>%
  data.matrix()

remove.rows <- rownames(depth.mat) %in% saved_ids
depth.mat <- depth.mat[!remove.rows,]
alt.mat <- alt.mat[!remove.rows,]

alpha <- 1e-6


clust <- makeCluster(16)
registerDoParallel(clust)


par_df <- foreach (i=1:nrow(depth.mat), .combine = 'rbind') %dopar% {
  pos_id <- rownames(depth.mat)[i]
  x <- alt.mat[i,]
  n <- depth.mat[i,]
  na.idx <- is.na(x) | is.na(n)
  x <- x[!na.idx]
  n <- n[!na.idx]
  a <- get_a(x,n,bw=0.25)
  fit <- estimate_par(x,n,a=a)
  gof <- GoF(x,n,fit$par,a)
  df <- data.frame(par1 = fit$par[1], par2 = fit$par[2], 
                   vcov1 = fit$vcov[1], vcov2 = fit$vcov[4], vcov12 = fit$vcov[2], 
                   a = a, 
                   ll1 = fit$ll[1], ll2 = fit$ll[2], 
                   convergence = fit$convergence,
                   gof_chisq_stat = gof['gof'], gof_df = gof['df'], gof_pval = gof['pval'])

  write.csv(df, sprintf('../results/par/%s.csv',pos_id),quote=F,row.names=F)
}


stopCluster(clust)


# aggregate individual files into one

par_files <- list.files('../results/par',full.names=T)
saved_ids <- gsub('\\.csv$','',basename(par_files))

par_df <- lapply(par_files,read.csv,stringsAsFactors=F) %>%
  data.table::rbindlist() %>%
  as.data.frame() %>%
  `rownames<-`(saved_ids) %>%
  .[rownames(alt.mat),]

write.csv(par_df, '../results/par.csv',quote=F)

