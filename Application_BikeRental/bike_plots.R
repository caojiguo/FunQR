library(lattice)
# quantiles of interest
tau <- seq(0.1,0.9,length.out = 30)
tau[5] <- 0.2
tau[16] <- 0.5

# time points
newgrids <- seq(0,1,length.out = 299)

# load the information about the bivariate splines
q2 <- readRDS('q2.RDS')
B.est <- readRDS('B.est.RDS')
B.est <- as.matrix(B.est)
k <- ncol(q2)

# generate B splines for the intercept
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2

# load the model fitting results
gama <- read.csv('output.txt',header=F)
# calculate the estimated coefficients of the bivariate splines
g <- gama$V1[1:k] - gama$V1[(1:k+k+nb)]

z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))

par(mfrow=c(2,2))
# Plot the 10%-quantile
bhat <- B.est[z[,2]==tau[1],]%*%q2%*%g
plot(seq(7,17,length.out = 299),bhat,ylim=c(-1000,600),type='l',xlab='Time',ylab='Estimation for 10%-quantile',cex.axis=2,cex.lab=1.5,lwd=2,xaxt='n')
axis(side=1,at=7:17,cex.axis=1.5)
lines(c(1:length(newgrids)),rep(0,length(newgrids)),col='red',lty=2)

# Plot the 20%-quantile
bhat <- B.est[z[,2]==tau[5],]%*%q2%*%g
plot(seq(7,17,length.out = 299),bhat,ylim=c(-1000,600),type='l',xlab='Time',ylab='Estimation for 20%-quantile',cex.axis=2,cex.lab=1.5,lwd=2,xaxt='n')
axis(side=1,at=7:17,cex.axis=1.5)
lines(c(1:length(newgrids)),rep(0,length(newgrids)),col='red',lty=2)

# Plot the 50%-quantile
bhat <- B.est[z[,2]==tau[16],]%*%q2%*%g
plot(seq(7,17,length.out = 299),bhat,ylim=c(-1000,600),type='l',xlab='Time',ylab='Estimation for 50%-quantile',cex.axis=2,cex.lab=1.5,lwd=2,xaxt='n')
axis(side=1,at=7:17,cex.axis=1.5)
lines(c(1:length(newgrids)),rep(0,length(newgrids)),col='red',lty=2)

# Plot the 90%-quantile
bhat <- B.est[z[,2]==tau[30],]%*%q2%*%g
plot(seq(7,17,length.out = 299),bhat,ylim=c(-1000,600),type='l',xlab='Time',ylab='Estimation for 90%-quantile',cex.axis=2,cex.lab=1.5,lwd=2,xaxt='n')
axis(side=1,at=7:17,cex.axis=1.5)
lines(c(1:length(newgrids)),rep(0,length(newgrids)),col='red',lty=2)

# Draw the heatmap
lt <- length(tau)
lg <- length(newgrids)
bhat <- B.est%*%q2%*%g
bhat.mtx <- matrix(bhat,nrow=lt,ncol=lg)
bhat.mtx <- t(bhat.mtx)
uu <- seq(7,17,length.out = 299)
vv <- tau
x <- rep(uu,each=length(vv))
y <- rep(vv,length(uu))
grid <- data.frame(cbind(x,y))
zz <- bhat

levelplot(zz~x*y,grid,col.regions=heat.colors(20),xlab=list(label='Time',cex=2),ylab=list(label='Quantile',cex=2),scales=list(x=list(at=c(7:17),cex=2),y=list(at=c(1:9)/10,cex=2)),colorkey = list(labels=list(cex=2)))

