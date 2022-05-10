# compared with Kato's method
library(psych)
library(quantreg)
wdir <- getwd()
# quantiles of interest
tau <- seq(0.2,0.8,length.out = 30)
# time points
t_grid <- seq(0.01,0.99,length.out = 199) 
newgrids <- seq(0.01,0.99,by=0.003)

tmp2 <- '2^0.5*cos(pi*t)'
for (i in 2:10)
{
  tmp1 <- paste0('(-1)^',eval(i+1),'*2^0.5*cos(',eval(i),'*pi*t)')
  ad <- paste0(eval(4*i^(-2)),'*',tmp1)
  tmp2 <- paste0(tmp2,'+',ad)
}
text <- paste0('function(t)',tmp2)
# rho1(t)
rho1 <- eval(parse(text=text))

# rho2(t) for simple case
rho2 <- function(t) 0

# rho2(t) for complex case
rho2 <- function(t) 2^0.5*cos(pi*t)

# true beta
tbeta <- function(t,alpha) rho1(t) +rho2(t)*qnorm(alpha,0,0.5)

# MSE and MAE
tmp1 <- c()
tmp2 <- c()
tmp3 <- c()
tmp4 <- c()

for (w in 1:100)
{
  setwd(paste0(wdir,'/simulation',w))
  mtx <- readRDS('mtx.RDS')
  daty <- readRDS('daty.RDS')
  # load the information about the bivariate splines
  q2 <- readRDS('q2.RDS')
  B.est <- readRDS('B.est.RDS')
  B.est <- as.matrix(B.est)
  k <- ncol(q2)
  # load the results of model fitting
  gama <- read.csv('output.txt',header=F)
  knots = quantile(tau, seq(0, 1, length = 31))
  norder= 2
  nknots = length(knots)
  nbasis = nknots + norder - 2
  # calculate the estimated coefficients of the bivariate splines
  g <- gama$V1[1:k] - gama$V1[(1:k+k+nbasis)]
  source('Estimate_FPCs.R')
  err1 <- c()
  err2 <- c()
  err3 <- c()
  err4 <- c()
  for (j in 1:length(tau)){
    z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
    fit <- rq(daty~scs,tau[j])
    betahat <- fit$coefficients[-1]%*%t(phihat)
    bhat <- B.est[z[,2]==tau[j],]%*%q2%*%g
    err1 <- c(err1,sum((tbeta(newgrids,tau[j])-betahat)^2))
    err2 <- c(err2,sum((tbeta(newgrids,tau[j])-bhat)^2))
    err3 <- c(err3,max(abs(tbeta(newgrids,tau[j])-betahat)))
    err4 <- c(err4,max(abs(tbeta(newgrids,tau[j])-bhat)))
  }
  tmp1 <- c(tmp1, mean(err1))
  tmp2 <- c(tmp2, mean(err2))
  tmp3 <- c(tmp3, mean(err3))
  tmp4 <- c(tmp4, mean(err4))
  print(w)
}
# MSE
mean(tmp2)/length(newgrids)
mean(tmp1)/length(newgrids)
mean(tmp2)/mean(tmp1)

# MAE
mean(tmp4)
mean(tmp3)
mean(tmp4)/mean(tmp3)

par(mfrow=c(1,2))

# draw the boxplot of MSE
boxplot(tmp2/length(newgrids),tmp1/length(newgrids),names=c('MSE of SFQR','MSE of FQR'),cex.axis=1.5,cex.lab=1.5)
# draw the boxplot of MAE
boxplot(tmp4,tmp3,names=c('MAE of SFQR','MAE of FQR'),cex.axis=1.5,cex.lab=1.5)
