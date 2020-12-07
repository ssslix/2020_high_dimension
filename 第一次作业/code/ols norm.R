##主体函数x维度p,样本n
object<-function(p,n){
  beta = rep(1,p)#初始化beta 设为m维全为1
  X=matrix(nrow=p,ncol=n)#初始化样本
  y = rep(0,n)
  for(j in 1:n)#按正态分布生成数据y
  {
    X[,j]=runif(p,min=-1,max=1)
    lambda=beta%*%X[,j]
    y[j]=rnorm(1,lambda,1)+rnorm(1,1,0.1)
  }
  #标准化X得到Z，中心化y
  mu=rowMeans(X)
  var=cov(t(X-mu))
  lamda =solve(eigen(var)$vectors)%*%var%*%(eigen(var)$vectors)
  lamda_sqrt <- diag(sqrt(diag(lamda)))
  var_sqrt <- (eigen(var)$vectors)%*%lamda_sqrt%*%solve(eigen(var)$vectors)
  sd=solve(var_sqrt)
  Z=sd%*%(X-mu)
  y=y-mean(y)
  #计算EYZ,最后返回sd%*%EYZ估计betahat
  EYZ= Z%*%y/n
  EYZ
  return(sd%*%EYZ)
}

object(10,1000)

experiment <- function(dim, n,times){
  true <- rep(1,dim)
  mse <- 0
  for (i in 1:times){
    beta <- object(dim,n)
    mse <- mse + sqrt(sum((beta - true)^2))
  }
  return(mse)
}

experiment(20,1000,100)

# dim :10 <- mse*100 = 1.75999
# dim: 15 <- mse*100 = 2.092107
# dim: 20 <- mse*100 = 2.497289
