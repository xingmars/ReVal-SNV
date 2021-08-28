
get_a <- function(x,n,q=0.8,low.lim=1e-3,...) {
  x1 <- x[x>0]
  n1 <- n[x>0]
  r1 <- log10(x1/n1)
  dd <- density(r1,...)
  local.min <- local.max <- NULL
  for (i in 2:(length(dd$y)-1)) {
    if (dd$y[i]<=dd$y[i-1] & dd$y[i]<=dd$y[i+1]) local.min <- c(local.min,dd$x[i])
    if (dd$y[i]>=dd$y[i-1] & dd$y[i]>=dd$y[i+1]) local.max <- c(local.max,dd$x[i])
  }
  noise.peak <- max(local.max[local.max < -2])
  local.min <- c(local.min,max(dd$x))
  c <- min(local.min[local.min>noise.peak])
  a <- 10^quantile(r1[r1<=c],q)
  a <- max(a,low.lim)
  a
}


call_snv <- function(x,n,alpha=1e-6) {
  a <- get_a(x,n,bw=0.25)
  fit <- estimate_par(x,n,a=a)
  par <- fit$par
  pvals <- 1 - sapply(1:length(x),function(i) pfunc(x[i]-1,n[i],par))
  group <- pvals <= alpha
  gof <- GoF(x,n,par,a)
  list(group=group,pval=pvals,par=par,vcov=fit$vcov[c(1,4,2)],a=a,ll=fit$ll,gof=gof)
}


GoF <- function(x,n,par,a,...) {
  x1 <- x[x>0 & x/n < a]
  n1 <- n[x>0 & x/n < a]
  r1 <- x1/n1
  breaks <- quantile(r1,seq(0,1,by=0.1))
  h <- hist(r1,plot=F,breaks=breaks)
  emp <- h$counts/sum(h$counts)
  
  est <- p_interval(x,n,par,breaks)
  est <- est/sum(est)
  # Note the denominator here is emp not est
  gof <- sum((emp-est)^2/emp)*length(x1)
  df <- length(est)-1
  pval <- 1 - pchisq(gof,df)
  c(gof=gof,df=df,pval=pval)
}


plot_GoF <- function(x,n,par,a,...) {
  x1 <- x[x>0 & x/n < a]
  n1 <- n[x>0 & x/n < a]
  r1 <- x1/n1
  h <- hist(r1,plot=F,...)
  breaks <- h$breaks
  breaks.lo <- breaks[-length(breaks)]
  breaks.hi <- breaks[-1]
  breaks.mid <- (breaks.lo+breaks.hi)/2
  bar.width <- (breaks.hi[1]-breaks.lo[1])/5
  
  est <- p_interval(x,n,par,breaks)
  est <- est/sum(est)*length(r1)
  
  # plot(h, freq=T, main=NULL, xlab='VAF')
  # lines((breaks[-1]+breaks[-length(breaks)])/2, est, type = "h", col='orange2', lwd=2)
  
  
  tmp.emp <- data.frame(emp=r1)
  tmp.theo <- data.frame(est=est,break.mid=breaks.mid)
  ggplot() +
    geom_histogram(data=tmp.emp,aes(x=emp,y=..count..),color='white',breaks=breaks) +
    geom_bar(data=tmp.theo,aes(x=break.mid,y=est),stat = "identity",fill='orange2',width=bar.width) +
    xlab('VAF (x/n)') +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank()
    ) 
}


p_interval <- function(x,n,par,breaks) {
  prob.mat <- matrix(NA,length(breaks)-1,length(n))
  
  for (i in 1:nrow(prob.mat)) {
    for (j in 1:ncol(prob.mat)) {
      prob.mat[i,j] <- pfunc(breaks[i+1]*n[j],n[j],par,log=F) - pfunc(breaks[i]*n[j],n[j],par,log=F)
    }
  }
  
  rowMeans(prob.mat)
}
