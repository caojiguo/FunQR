library(fda)
library(MASS)
# load the data
hgtm <- read.csv('growth_male.csv',header=T)
hgtf <- read.csv('growth_female.csv',header=T)
n <- nrow(hgtm)
daty <- as.numeric(c(hgtm[n,-1],hgtf[n,-1]))
mtx <- cbind(hgtm[1:19,-1],hgtf[1:19,-1])
mtx <- t(mtx)
gcv <-c()
# i is the number of knots; i+4-2 is the number of basis functions
grids <- c(0:18/18)
for (i in 2:10 ){
  # order of basis function
  norder = 5
  knots = quantile(grids, seq(0,1,length = i))
  basisobj = create.bspline.basis(rangeval = c(min(grids), max(grids)), norder = norder, breaks = knots)
  fdParobj = fdPar(fdobj=basisobj)
  output = smooth.basis(grids,t(mtx), fdParobj)
  tmp = output$gcv
  gcv <- rbind(gcv,c(i,tmp))
}
plot(gcv[,1], gcv[,2], type = "l", xlab = "# of Basis Functions", ylab = "GCV criterion")

knots = quantile(grids, seq(0, 1, length = 9))
nknots = length(knots)
nbasis = nknots + norder - 2
genbasis = create.bspline.basis(rangeval = c(min(grids), max(grids)), nbasis = nbasis,norder = norder,breaks = knots)
harmLfd = int2Lfd(2)
xfdpar = fdPar(genbasis,harmLfd,10^(-1))
xfd = smooth.basis(grids,t(mtx),xfdpar)

# apply FPCA on the data
xfpca = pca.fd(xfd$fd, nharm = 5)
harmfd = xfpca$harmonics
newgrids <- seq(0,1,length.out = 99)
harmvals = eval.fd(newgrids,harmfd)
mharmvals = eval.fd(grids,harmfd)
mean_mtx <- eval.fd(grids,xfpca$meanfd)
mmtx <- matrix(mean_mtx,nrow(mtx),length(grids),byrow=T)

# estimate FPC scores
scs <- c()
for (i in 1:5)
{
  scs <- rbind(scs,mharmvals[,i]%*%t(mtx-mmtx)/ncol(mtx))
}
scs <- t(scs)
saveRDS(object=daty, file='daty.RDS')
saveRDS(object=scs,file='scs.RDS')
saveRDS(object = harmfd,file='harmfd.RDS')