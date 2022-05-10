library(fda)
library(MASS)
mtx <- readRDS('mtx.RDS')
gcv <-c()
# i is the number of knots; i+4-2 is the number of basis functions
t_grid <- seq(0.01,0.99,length.out = 199) 
grids <- t_grid
for (i in 2:10){
  # order of basis function
  norder = 5
  knots = quantile(grids, seq(0,1,length = i))
  basisobj = create.bspline.basis(rangeval = c(min(grids), max(grids)), norder = norder, breaks = knots)
  fdParobj = fdPar(fdobj=basisobj)
  # Without penalty
  output = smooth.basis(grids, t(mtx), fdParobj)
  # GCV output
  tmp = output$gcv
  gcv <- rbind(gcv,c(i,tmp))
}

# Choose the number of knots based on GCV (elbow point or the first local minimum)
plot(gcv[,1], gcv[,2], type = "l", xlab = "# of Basis Functions", ylab = "GCV criterion")

# apply fpca on the selected observations
knots = quantile(grids, seq(0, 1, length = 5))
nknots = length(knots)
nbasis = nknots + norder - 2
genbasis = create.bspline.basis(rangeval = c(min(grids), max(grids)), nbasis = nbasis,norder = norder,breaks = knots)
harmLfd = int2Lfd(2)
xfdPar = fdPar(genbasis,harmLfd,10^(-1))
xfd = smooth.basis(grids,t(mtx),xfdPar)
xfpc = pca.fd(xfd$fd, nharm = 3)
harmfd = xfpc$harmonics
newgrids <- seq(0.01,0.99,by=0.003)
harmvals = eval.fd(newgrids,harmfd)
mharmvals = eval.fd(grids,harmfd)

# estimated mean curve
mean_est <- eval.fd(newgrids,xfpc$meanfd)
mean_mtx <- eval.fd(grids,xfpc$meanfd)
meanmtx <- matrix(mean_mtx,nrow(mtx),length(grids),byrow=T)

# estimated FPCs
phihat <- matrix(0,length(newgrids),ncol(harmvals))
for (i in 1:3)
{
  phi <- function(t) (-1)^eval(i+1)*2^0.5*cos(i*pi*t)
  if (mean(mharmvals[,i]*phi(grids))>0) {
    phihat[,i] <- harmvals[,i]
  }else{phihat[,i] <- -harmvals[,i]
  mharmvals[,i] = -mharmvals[,i]}
}

# estimated scores
scs <- c()
for (i in 1:3)
{
  scs <- rbind(scs,mharmvals[,i]%*%t(mtx-meanmtx)/ncol(mtx))
}
scs <- t(scs)
saveRDS(object=phihat,file='phihat.RDS')
saveRDS(object=scs,file='scs_est.RDS')