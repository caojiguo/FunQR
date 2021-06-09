library(fda)
# candidate values for tuning parameters
lam1 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
lam2 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
# quantiles of interest
tau <- seq(0.2,0.8,length.out = 20)
tau[3] <- 0.25
tau[10] <- 0.5 
tau[18] <- 0.75
# time points
newgrids <- seq(0,1,length.out = 299)
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
qloss <- function(tau,y)
{
  y*(tau-ifelse(y<0,1,0))
}
# number of folds for cross validation
nfolds <- 10
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
# number of B-spline basis functions
nb = nknots + norder - 2
qls <- array(0,dim=c(length(lam1),nfolds,length(tau)))
daty <- readRDS('daty.RDS')
ly <- length(daty)
ncross <- floor(ly/nfolds)
# load the information of the bivariate splines
q2 <- readRDS('q2.RDS')
B.est <- readRDS('B.est.RDS')
k <- ncol(q2)
dirpath <- getwd()
for (w in 1:nfolds)
{
    foldid <- w-1
    tmpp <- c()
    setwd(paste0(dirpath,'/fold_',w))
    A <- readRDS(file = 'cvp.RDS')
    A <- as.matrix(A)
    for (l in 1:length(lam1))
    {
      # load the results of model fitting
      gama <- read.csv(paste0("output",l,".txt"),header=F)
      # calculate the estimated coefficients of the bivariate splines
      g <- gama$V1[1:(k+nb)] - gama$V1[(1:(k+nb)+k+nb)]
      cvfold <- c(1:ncross+foldid*ncross)
      for (m in 1:length(tau))
      {
        tt <- tau[m]
        qls[l,foldid+1,m] <- sum(qloss(tt,(daty[cvfold]-A[c(1:ncross+(m-1)*ncross),]%*%g)))
      }
    }
}
setwd(dirpath)
saveRDS(object=qls,file ='cv_check.RDS')
cv_qls <- apply(apply(qls,c(1,2),sum),1,mean)

# decide the value for the tuning paramter of the first penalty
lam.index <- which(cv_qls==min(cv_qls))
saveRDS(object=lam.index,file ='lam.index.RDS')
