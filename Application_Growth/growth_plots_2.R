library(quantreg)
library(lattice)
library(BPST)
library(latex2exp)
library(fda)
setwd('~/Documents/fqr/CODE4FunQuanR/Application_Growth/results_of_model_fitting/')
# quantiles of interest
tau <- seq(0.2,0.8,length.out = 20)
tau[3] <- 0.25
tau[10] <- 0.5
tau[18] <- 0.75

# load the information about the triangulation
V.est <- readRDS(file='nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')
d <- 2 # degree of the polynomial
r <- 1 # smoothness condition

# load the information about the FPCA result
harmfd <- readRDS(file='harmfd.RDS')
scs <- readRDS(file='scs.RDS')
daty <- readRDS(file='daty.RDS')
phihat <- eval.fd(newgrids,harmfd)
# time points
newgrids <- seq(1,12,length.out = 299)

# compare the monotonicity of the quantiles
betahat <- c()

# generate B-spline basis for the intercept
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

# load the information about the bivariate splines
B.est <- readRDS('B.est.RDS')
q2 <- readRDS('q2.RDS')
k <- ncol(q2)

# load the model fitting results
gama <- read.csv('output.txt',header=F)

# estimated coefficients of the bivariate basis over triangulation for the slope function \beta(t,u)
g <- gama$V1[1:k] - gama$V1[(1:k+k+nb)]

# estimate coefficients of the B-splines for the intercept
g0 <- gama$V1[(k+1):(k+nb)] - gama$V1[(k*2+nb+1):(k*2+nb*2)]

# calculate quantile estimations for the proposed method
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
lng <- length(newgrids)
quan_est <- matrix(0,93,20)
for (i in 1:20)
{
  tt <- z[z[,2]==tau[i],1]
  l <- length(tt)
  if (l%%2==0) s <- c(1,rep(c(4,2),l/2-1),1)
  if (l%%2==1) s <- c(1,4,rep(c(2,4),(l-3)/2),1)
  bb <- B.est[z[,2]==tau[i],]
  p <- t(s*phihat)%*%bb/(3*lng)
  X1 <- cbind(scs%*%p%*%q2,matrix(bb0[i,],nrow=nrow(scs),ncol=ncol(bb0),byrow=T))
  quan_est[,i] <- as.numeric(X1%*%c(g,g0))
}

# calculate quantile estimations for the conventional linear functional quantile regression method
quan_est2 <- matrix(0,93,20)
for (k in 1:20)
{
  fit <- rq(daty~scs, tau[k])
  betahat <- rbind(betahat, fit$coefficients[-1]%*%t(phihat))
  quan_est2[,k] <- fit$fitted.values
}

# plot the estimated quantiles of the 37th and the 67th observations based on the conventional method and the proposed method
par(mfrow=c(2,2))
par(mar=c(5,5.5,4,1))
plot(tau, quan_est2[37,],type='l',xlab='u',ylab=TeX('$Q^*_{37}(u)$'),cex.axis=1.5,cex.lab=2,lwd=2)
plot(quan_est[37,],type='l',xlab='u',ylab=TeX('$Q_{37}(u)$'),cex.axis=1.5,cex.lab=2,lwd=2)
plot(quan_est2[67,],type='l',xlab='u',ylab=TeX('$Q^*_{67}(u)$'),cex.axis=1.5,cex.lab=2,lwd=2)
plot(quan_est[67,],type='l',xlab='u',ylab=TeX('$Q_{67}(u)$'),cex.axis=1.5,cex.lab=2,lwd=2)

# Draw the heatmap for the slope function estimation obtained from the conventional method
lt <- length(tau)
lg <- length(newgrids)
bhat.mtx <- matrix(betahat,nrow=lt,ncol=lg)
bhat.mtx <- t(bhat.mtx)
uu <- seq(1,12,length.out = 299)
vv <- tau
x <- rep(uu,each=length(vv))
y <- rep(vv,length(uu))
grid <- data.frame(cbind(x,y))
zz <- betahat

levelplot(zz~x*y,grid,col.regions=heat.colors(20),xlab=list(label='Age',cex=2),ylab=list(label='Quantile',cex=2),scales=list(x=list(at=c(1:12),cex=2),y=list(at=c(1:9)/10,cex=2)),colorkey = list(labels=list(cex=2)))