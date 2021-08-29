dbetabinom <- function(x,n,par,...) {
  prob <- par[1]
  rho <- par[2]
  VGAM::dbetabinom(x,n,prob,rho,...)
}

pbetabinom <- function(x,n,par,...) {
  prob <- par[1]
  rho <- par[2]
  VGAM::pbetabinom(x,n,prob,rho,...)
}

qbetabinom <- function(q,n,par,...) {
  prob <- par[1]
  rho <- par[2]
  for (x in 0:n) {
    cdf <- VGAM::pbetabinom(x,n,prob,rho,...)
    if (cdf>=q) return(x)
  }
}

dtbetabinom <- function(x,n,par,a,log=T) {
  if (length(n)==1) n <- rep(n,length(x))
  prob <- par[1]
  rho <- par[2]
  t <- a*n
  if (log) {
    d <- sapply(1:length(x),function(i) dbetabinom(x[i],n[i],par,log=T) - pbetabinom(t[i],n[i],par,log.p=T))
    d[x>t] <- -Inf
  } else {
    d <- sapply(1:length(x),function(i) dbetabinom(x[i],n[i],par.log=F)/pbetabinom(t[i],n[i],par,log=F))
    d[x>t] <- 0
  }
  d
}

dfunc <- dbetabinom
pfunc <- pbetabinom
qfunc <- qbetabinom
dtfunc <- dtbetabinom

estimate_par <- function(x,n,a) {
  minuslogl <- function(par) {
    d <- dtfunc(x,n,par,a,log=T)
    -sum(d[x/n<a])
  }
  
  fit <- optim(c(1e-4,1e-4), minuslogl,method='Nelder-Mead', control=list(reltol=1e-8),hessian=F)
  fit$vcov <- tryCatch(solve(numDeriv::hessian(minuslogl,fit$par)), error = function(x) NA)
  
  minuslogl1 <- function(par) {
    d <- dtfunc(x,n,par,a,log=T)
    -sum(d[x/n<a & x>0])
  }
  ll1 <- -minuslogl1(fit$par)
  fit$ll <- c(-fit$value,ll1)
  
  fit
}
