library(fda)
# candidates for tuning parameters
lam1 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
lam2 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)

# quantiles of interest
tau <- seq(0.1,0.9,length.out = 30)
tau[5] <- 0.2
tau[16] <- 0.5

# time points
newgrids <- seq(0,1,length.out = 299)
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))

qloss <- function(tau,y)
{
  y*(tau-ifelse(y<0,1,0))
}
# number of folds
nfolds <- 10
qls <- array(0,dim=c(length(lam2),nfolds,length(tau)))

# load the information about the bivariate splines
q2 <- readRDS('q2.RDS')
B.est <- readRDS('B.est.RDS')

daty <- readRDS('daty.RDS')
ly <- length(daty)
ncross <- floor(ly/nfolds)

# generate B-splines
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

# summarize the results of cross validation to obtain the "best" value for the tuning parameter of the second penalty
dirpath <- getwd()
for (w in 1:nfolds)
{
    foldid <- w-1
    tmpp <- c()
    setwd(paste0(dirpath,'/second_fold_',w))
    A<- readRDS(file = 'cvp.RDS')
    A <- as.matrix(A)
    k <- ncol(q2)
    for (l in 1:length(lam2))
    {
      gama <- read.csv(paste0("output",l,".txt"),header=F)
      g <- gama$V1[1:(k+nb)] - gama$V1[(1:(k+nb)+k+nb)]
      cvfold <- c(1:ncross+foldid*ncross)
      for (m in 1:length(tau))
      {
        tt <- tau[m]
        qls[l,foldid+1,m] <-  sum(qloss(tt,(daty[cvfold]-A[c(1:ncross+(m-1)*ncross),]%*%g)))
      }
    }
}
setwd(dirpath)
saveRDS(object=qls,file ='cv_check2.RDS')
cv_qls <- apply(apply(qls,c(1,2),sum),1,mean)

# decide the value of the tuning parameter of the second penalty.
lam.index2 <- which(cv_qls==min(cv_qls))

# Based on the chosen values of two tuning parameters, fit the model again on the whole data set.
source('rd-Fitting.R')
writeMat('penalty1.mat',penalty1=penalty1)
writeMat('penalty2.mat',penalty2=penalty2)
writeMat('c.mat',cc = cc)
writeMat('Ales.mat',Ales = B)
writeMat('Aeq.mat',Aeq = A)
writeMat('beq.mat',beq=y)
writeMat('bles.mat',bles=yy)
writeMat('lb.mat',lb=rep(0,length(cc)))
writeMat('nb.mat',nb=nb)
lam.index1 <- readRDS('lam.index.RDS')
lam1 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
lam2 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
writeMat('p1.mat',p1=lam1[lam.index1])
writeMat('p2.mat',p2=lam2[lam.index2])

# generate the MATLAB script for the quadratic programming problem that needs to be submitted to MATLAB/2018b to solve 
zztmp <- paste0('sim.m')
zz <- file(zztmp, open = "wt")
texttmp <- paste0("load('Ales.mat');
load('c.mat');
load('Aeq.mat');
load('penalty1.mat');
load('penalty2.mat');
load('beq.mat');
beq=double(beq);
load('bles.mat');
load('lb.mat');
load('p1.mat');
load('p2.mat');
load('nb.mat')
[nr,nc]=size(penalty1);
tmp1 = [penalty1,zeros(nr,nb),-penalty1,zeros(nr,nb)];
[nr,nc]=size(tmp1);
tmp2=zeros(nb,nc);
RP=[tmp1;tmp2];
RP = [RP;-RP];
[nr,nc]=size(RP);
lc = length(lb);
RP = [RP,zeros(nr,lc-nc)];
[nr,nc]=size(RP);
RP1=[RP;zeros(lc-nr,nc)];
RP1=p1*RP1;
[nr,nc]=size(penalty2);
tmp1 = [penalty2,zeros(nr,nb),-penalty2,zeros(nr,nb)];
[nr,nc]=size(tmp1);
tmp2=zeros(nb,nc);
RP=[tmp1;tmp2];
RP = [RP;-RP];
[nr,nc]=size(RP);
lc = length(lb);
RP = [RP,zeros(nr,lc-nc)];
[nr,nc]=size(RP);
RP2=[RP;zeros(lc-nr,nc)];
RP2=p2*RP2;
RP = RP1 +RP2;
                  options=optimoptions('quadprog','linearsolver','sparse');
                  [x,fval] = quadprog(RP,cc/2,Ales,bles,Aeq,beq,lb,[],[],options);
                  dlmwrite('output.txt', x, 'delimiter', ',', 'precision', 7);")
write(texttmp, file=zz)
close(zz)
