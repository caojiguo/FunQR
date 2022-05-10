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
# time points
t_grid <- seq(0.01,0.99,length.out = 199) 
# number of observations
nobs <- 400

# x(t)
mtx <- matrix(0,nobs,length(t_grid))
# y
daty <- c()

for (j in 1:nobs)
{
  tmp1 <- c()
  tmp2 <- c()
  sc <- c()
  for (i in 1:10)
  {
    z <- runif(1,-3^0.5,3^0.5)
    tmp1 <- paste0(z,'*',eval(i),'^(-1)','*','2^0.5*cos(',eval(i),'*pi*t)')
    tmp2 <- paste0(tmp2,'+',tmp1)
    sc <- c(sc,z)
  }
  tmp <- paste0(tmp2,'+3*6^0.5*cos(pi*t)')
  text <- paste0('function(t)',tmp)
  # X(t)
  X <- eval(parse(text=text))
  mtx[j,] <- X(t_grid)
  tmp1 <- c()
  tmp2 <- sc[1]
  for (i in 2:10)
  {
    tmp1 <- paste0('(-1)^',eval(i+1),'*',eval(4*i^(-2)),'*',eval(i),'^(-1)','*',eval(sc[i]))
    tmp2 <- paste0(tmp2,'+',tmp1)
  }
  # y
  daty <- c(daty,eval(parse(text=tmp2)) +  3*3^0.5+ rnorm(1))
  print(j)
}
saveRDS(object=mtx,file='mtx.RDS')
saveRDS(object=daty,file='daty.RDS')