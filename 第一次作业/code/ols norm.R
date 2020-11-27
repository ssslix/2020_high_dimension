object<-function(m,n){
  beta = rep(1,m)
  X=matrix(nrow=m,ncol=n)
  y = rep(0,n)
  for(j in 1:n)
  {
    X[,j]=runif(m,min=-10,max=10)
    lambda=beta%*%X[,j]
    y[j]=rnorm(1,lambda,1)+rnorm(1,1,0.1)
  }
  mu=rowMeans(X)
  var=cov(t(X-mu))
  lamda =solve(eigen(var)$vectors)%*%var%*%(eigen(var)$vectors)
  lamda_sqrt <- diag(sqrt(diag(lamda)))
  var_sqrt <- (eigen(var)$vectors)%*%lamda_sqrt%*%solve(eigen(var)$vectors)
  sd=solve(var_sqrt)
  Z=sd%*%(X-mu)
  y=y-mean(y)
  EYZ= Z%*%y/n
  EYZ
  return(sd%*%EYZ)
}

object(100,1000)
