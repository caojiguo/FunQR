mdir <- getwd()
nfolds <- 10
for (id in 1:nfolds)
{
  subdir <- paste0("fold_", id)
  id <- id -1
  source('cv_Fitting.R')
  dir.create(file.path(mdir, subdir), showWarnings = FALSE)
  setwd(file.path(mdir, subdir))
  saveRDS(object=cvp, file='cvp.RDS')
  
  # load the coefficients of the quadratic programming problem
  writeMat('penalty1.mat',penalty1=penalty1)
  writeMat('penalty2.mat',penalty2=penalty2)
  writeMat('c.mat',cc = cc)
  writeMat('Ales.mat',Ales = B)
  writeMat('Aeq.mat',Aeq = A)
  writeMat('beq.mat',beq=y)
  writeMat('bles.mat',bles=yy)
  writeMat('lb.mat',lb=rep(0,length(cc)))
  # quantiles of interest
  tau <- seq(0.2,0.8,length.out = 20)
  tau[3] <- 0.25
  tau[10] <- 0.5 
  tau[18] <- 0.75
  knots <- quantile(tau,seq(0,1,length=31))
  norder= 2
  nknots = length(knots)
  # number of B-spline basis functions
  nb = nknots + norder - 2
  writeMat('nb.mat',nb=nb)
  # candidate values for the tuning paramter of the first penalty
  lam1 <- c(10^-6,10^-5.5,10^-5,10^-4.5,10^-4)
  
  # generate the MATLAB scripts for the quadratic programming problems that need to be submitted to MATLAB/2018b to solve
  for (i in 1:length(lam1)){
    zztmp <- paste0('sim',i,'.m')
    zz <- file(zztmp, open = "wt")
    texttmp <- paste0("load('Ales.mat');
                  load('c.mat');
                  load('Aeq.mat');
                  load('penalty1.mat');
                  load('penalty2.mat');
                  load('beq.mat');
                  beq = double(beq);
                  load('bles.mat');
                  load('lb.mat');
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
                  options=optimoptions('quadprog','linearsolver','sparse');
                  x = quadprog(",lam1[i],"*RP1,cc/2,Ales,bles,Aeq,beq,lb,[],[],options);
                  dlmwrite('output",i,".txt', x, 'delimiter', ',', 'precision', 7);")
    write(texttmp, file=zz)
    close(zz)
  }
  setwd(mdir)
}
