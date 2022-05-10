library(lattice)
# quantiles of interest
tau <- seq(0.2,0.8,length.out = 20)
tau[3] <- 0.25
tau[10] <- 0.5
tau[18] <- 0.75

# time points
newgrids <- seq(1,12,length.out = 299)

# load the information of the bivariate splines
q2 <- readRDS('q2.RDS')
B.est <- readRDS('B.est.RDS')
B.est <- as.matrix(B.est)
k <- ncol(q2)

# load the results of model fitting
gama <- read.csv('output.txt',header=F)
knots <- quantile(tau,seq(0,1,length=31))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
# calculate the estimated coefficients of the bivariate splines
g <- gama$V1[1:k] - gama$V1[(1:k+k+nb)]
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))

par(mar=c(5,4.5,4,1))
# Plot the 10%-quantile
bhat <- B.est[z[,2]==tau[1],]%*%q2%*%g
plot(seq(1,12,length.out = 299),bhat,type='l',ylim=c(-60,20),xlab='Age (Years)',ylab='Estimation of the Slope Function for Multiple Quantiles',lwd=3,cex.axis = 2,cex.lab=2,xaxt='n')
axis(side=1,at=1:12,cex.axis=2)
# Plot the 25%-quantile
bhat <- B.est[z[,2]==tau[3],]%*%q2%*%g
lines(seq(1,12,length.out = 299),bhat,type='l',lty=2,col=2, lwd=3,cex.axis = 1.5,cex.lab=1.5)
# Plot the 50%-quantile
bhat <- B.est[z[,2]==tau[10],]%*%q2%*%g
lines(seq(1,12,length.out = 299),bhat,type='l',lty=3,col=3, lwd=3,cex.axis = 1.5,cex.lab=1.5)
# Plot the 75%-quantile
bhat <- B.est[z[,2]==tau[18],]%*%q2%*%g
lines(seq(1,12,length.out = 299),bhat,type='l',lty=4,col=4,lwd=3,cex.axis = 1.5,cex.lab=1.5)
# Plot the 90%-quantile
bhat <- B.est[z[,2]==tau[20],]%*%q2%*%g
lines(seq(1,12,length.out = 299),bhat,type='l',lty=5,col=5,lwd=3,cex.axis = 1.5,cex.lab=1.5)
legend("bottomleft",legend=c("20%","25%","50%","75%","80%"),
       col=c(1:5), lty=1:5, cex=1.5, lwd=3,
       title="Quantile Levels", text.font=4, ncol=3)
abline(0,0,lty =6, col=6)

# plot the estimated slope function at t=5
bhat <- B.est[z[,1]==newgrids[110],]%*%q2%*%g
plot(tau,as.numeric(bhat),type='l',xlim=c(0.2,0.8),xlab='Quantile',ylab='Estimation of the slope function at t=5',cex.axis=2,cex.lab=2,lwd=3)

# Draw the heatmap
newgrids <- seq(1,12, length.out = 299)
lt <- length(tau)
lg <- length(newgrids)
bhat <- B.est%*%q2%*%g
bhat.mtx <- matrix(bhat,nrow=lt,ncol=lg)
bhat.mtx <- t(bhat.mtx)
uu <- newgrids
vv <- tau
x <- rep(uu,each=length(vv))
y <- rep(vv,length(uu))
grid <- data.frame(cbind(x,y))
zz <- bhat

levelplot(zz~x*y,grid,col.regions=heat.colors(100),xlab=list(label='Age (Years)',cex=2),ylab=list(label='Quantile',cex=2),scales=list(x=list(at=c(1:12),cex=2),y=list(at=c(1:9)/10,cex=2)),colorkey = list(labels=list(cex=2)))
title(main = "Heat plot of estimated slope function for Berkeley growth data", font.main = 4)
