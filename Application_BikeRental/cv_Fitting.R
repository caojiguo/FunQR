library(stats)
library(BPST)
library(Matrix)
library(fda)
library(lpSolve)
library(R.matlab)
library(MASS)
# When using clusters to perform the cross validition, we label the folds by the job id on the clusters.
print(id)

# quantiles of interest
tau <- seq(0.1,0.9,length.out = 30)
tau[5] <- 0.2
tau[16] <- 0.5

# time points
newgrids <- seq(0,1,length.out = 299)

# load the information about the triangulation and choose the parmeters for the bivariate splines
V.est <- readRDS(file = 'nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')
d <- 2
r <- 1

# load the information about the FPCA result
scs <- readRDS('scs.RDS')
daty <- readRDS('daty.RDS')
harmfd <- readRDS(file='harmfd.RDS')
phihat <- eval.fd(newgrids,harmfd)


# seperate the data into 10 folds and select one of them
nfolds <- 10
ly <- length(daty)
ncross <- floor(ly/nfolds)
fold <- c(1:ly)[-c(1:ncross+id*ncross)]

# generate B-spline basis for the intercept
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

# generate bivariate Bernstein polynomials over the given triangulation for the slope function $\beta(t,u)$
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
Bfull.est <- basis(V.est,Tr.est,d,r,z)
Rough <- Bfull.est$K
B.est <- Bfull.est$B
ind <- Bfull.est$Ind.inside
z <- z[ind,]

# approximate the integral of $\int \phi(t)*b_j(t,u)dt$ by Simpson's rule and calculate the matrix P(u_r) for all r=1, ..., 30
tt <- z[z[,2]==tau[1],1]
l <- length(tt)
if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
bb <- B.est[z[,2]==tau[1],]
lng <- length(newgrids)
p <- t(s*phihat)%*%bb/3*(1/lng)
q2 <- Bfull.est$Q2

# approximate the integral of $\int \phi(t)*b_j(t,u)dt$ by Simpson's rule and calculate the matrix P(u_r) for all r=1, ..., 30

cvscs <- scs[c(1:ncross+id*ncross),] # "id" is the index for the folds during the cross validtion, which starts from 0

X <- cbind(cvscs%*%p%*%q2,matrix(bb0[1,],nrow=nrow(cvscs),ncol=ncol(bb0),byrow=T))
lt <- length(tau)
n <- nrow(X)
A <- X
for (i in 1:(lt-1))
{
  # Simpson's rule
  tt <- z[z[,2]==tau[i+1],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i+1],]
  p <- t(s*phihat)%*%bb/(3*lng)
  X <- cbind(cvscs%*%p%*%q2,matrix(bb0[i+1,],nrow=nrow(cvscs),ncol=ncol(bb0),byrow=T))
  n <- nrow(X)
  tmp <- cbind(X,-X,matrix(0,n,2*i*n),diag(n),-diag(n),matrix(0,n,(lt-1-i)*2*n))
  A <- rbind(A,X)
}
cvp <- A

tt <- z[z[,2]==tau[1],1]
l <- length(tt)
if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
bb <- B.est[z[,2]==tau[1],]
lng <- length(newgrids)
p <- t(s*phihat)%*%bb/3*(1/lng)
scs <- scs[fold,]

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
  p <- t(s*phihat)%*%bb/(3*lng)
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
  p <- t(s*phihat)%*%bb/(3*lng)
  X1 <- cbind(scs%*%p%*%q2,matrix(bb0[i,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
  # Simpson's rule
  tt <- z[z[,2]==tau[i+1],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i+1],]
  p <- t(s*phihat)%*%bb/(3*lng)
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
y <- rep(daty[fold],lt)
yy <- rep(0,length(daty[fold])*(lt-1))

penalty1 <- t(q2)%*%t(B.est)%*%B.est%*%q2
penalty2 <- t(q2)%*%Rough%*%q2

# save the information about the bivariate splines
saveRDS(object=B.est,file='B.est.RDS')
saveRDS(object=q2,file='q2.RDS')
