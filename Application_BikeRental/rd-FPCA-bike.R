library(dplyr)
library(fda)
library(MASS)
# load the dataset
wdir <- getwd()
setwd(paste0(wdir,'/Bike-Sharing-Dataset/'))
tmpy <- read.csv('day.csv',header=T)
tmpx <- read.csv('hour.csv', header=T)

# data cleaning
tmp=tmpx%>%group_by(dteday)%>%summarize(ncount=length(unique(hr)),rain=max(weathersit)>=3)
tmpx=tmpx%>%filter(dteday%in%tmp$dteday[tmp$ncount==24&tmp$rain==F])
tmpy=tmpy%>%filter(dteday%in%tmp$dteday[tmp$ncount==24&tmp$rain==F])
datx <- matrix(tmpx$temp[tmpx$weekday%in%c(0, 6)],ncol=24, byrow=T)
hr_cnt <- matrix(tmpx$cnt[tmpx$weekday%in%c(0, 6)],ncol=24, byrow=T)
t_int <- c(7:18)
datx <- datx*47-8
daty <- apply(hr_cnt[,t_int],1,sum)

# select the observations between 7:00 to 17:00
mtx <- datx[,t_int]

# apply fpca on the selected observations
norder = 5
grids <- seq(0,1,length.out = length(t_int))
knots = quantile(grids, seq(0, 1, length.out = 5))
nknots = length(knots)
nbasis = nknots + norder - 2
genbasis = create.bspline.basis(rangeval = c(min(grids), max(grids)), nbasis = nbasis,norder = norder,breaks = knots)
harmLfd = int2Lfd(2)
xfdPar = fdPar(genbasis,harmLfd,10^-3)
xfd = smooth.basis(grids,t(mtx),xfdPar)
xfpca = pca.fd(xfd$fd, nharm = 3)
harmfd = xfpca$harmonics
newgrids <- seq(0,1,length.out = 299)
harmvals = eval.fd(newgrids,harmfd)
mharmvals = eval.fd(grids,harmfd)
mean_mtx <- eval.fd(grids,xfpca$meanfd)
mmtx <- matrix(mean_mtx,nrow(mtx),length(grids),byrow=T)

# calculated the estimated scores
scs <- c()
for (i in 1:3)
{
  scs <- rbind(scs,mharmvals[,i]%*%t(mtx-mmtx)/ncol(mtx))
}
scs <- t(scs)

setwd(wdir)
saveRDS(object=daty, file='daty.RDS')
saveRDS(object=scs,file='scs.RDS')
saveRDS(object = harmfd,file='harmfd.RDS')
