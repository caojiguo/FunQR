library(stats)
library(BPST)
# the R package "BPST" can be downloaded from
# https://github.com/funstatpackages/BPST
library(Matrix)
library(fda)
library(lpSolve)
library(R.matlab)
library(MASS)
# quantiles of interest
tau <- seq(0.2,0.8,length.out = 20)
tau[3] <- 0.25
tau[10] <- 0.5 
tau[18] <- 0.75
# load the information about the triangulation
V.est <- readRDS(file='nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')

# load the information about the FPCA result
harmfd <- readRDS(file='harmfd.RDS')
scs <- readRDS(file='scs.RDS')
daty <- readRDS(file='daty.RDS')

# generate bivariate Bernstein polynomials over the given triangulation for the slope function $\beta(t,u)$
newgrids <- seq(0,1,length.out = 299)
d <- 2
r <- 1
phihat <- eval.fd(newgrids,harmfd)
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
Bfull.est <- basis(V.est,Tr.est,d,r,z)
Rough <- energy(V.est,Tr.est,d)
B.est <- Bfull.est$B
ind <- Bfull.est$Ind.inside
z <- z[ind,]

# Approximate the integral of phi(t)*B(t,\tau)
tt <- z[z[,2]==tau[1],1]
l <- length(tt)
if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
bb <- B.est[z[,2]==tau[1],]
p <- t(s*phihat)%*%bb/3*(1/l)

# generate B-spline basis for the intercept
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)
q2 <- Bfull.est$Q2

# The classcal quantile regression can be transfered into a linear programming problem and our proposed penalized quantile regression can be transfered into a quadratic programming problem.

# set up the quadratic programming problem
X <- cbind(scs%*%p%*%q2,matrix(bb0[1,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
lt <- length(tau)
n <- nrow(X)
A <- cbind(X,-X,diag(n),-diag(n),matrix(0,n,(lt-1)*2*n))
for (i in 1:(lt-1))
{
  # Simpson's rule
  tt <- z[z[,2]==tau[i+1],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i+1],]
  p <- t(s*phihat)%*%bb/(3*l)
  X <- cbind(scs%*%p%*%q2,matrix(bb0[i+1,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
  n <- nrow(X)
  tmp <- cbind(X,-X,matrix(0,n,2*i*n),diag(n),-diag(n),matrix(0,n,(lt-1-i)*2*n))
  A <- rbind(A,tmp)
}

######### constraints for monotonicity ##########
B <- c()
for (i in 1:(lt-1))
{
  # Simpson's rule
  tt <- z[z[,2]==tau[i],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i],]
  p <- t(s*phihat)%*%bb/(3*l)
  X1 <- cbind(scs%*%p%*%q2,matrix(bb0[i,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
  # Simpson's rule
  tt <- z[z[,2]==tau[i+1],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i+1],]
  p <- t(s*phihat)%*%bb/(3*l)
  X2 <- cbind(scs%*%p%*%q2,matrix(bb0[i+1,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
  X <- X1-X2
  tmp <- cbind(X,-X,matrix(0,n,lt*2*n))
  B <- rbind(B,tmp)
}

# calculate the other inputs for the equivalent quadratic programming problem
cc <- rep(0,2*ncol(X))
for (i in 1:lt)
{
  cc <- c(cc,tau[i]*rep(1,n),(1-tau[i])*rep(1,n))
}

lc <- length(cc)
y <- rep(daty,lt)
yy <- rep(0,length(daty)*(lt-1))

penalty1 <- t(q2)%*%t(B.est)%*%B.est%*%q2
penalty2 <- t(q2)%*%Rough%*%q2

# save the information about the bivariate splines
saveRDS(object=B.est,file='B.est.RDS')
saveRDS(object=q2,file='q2.RDS')
