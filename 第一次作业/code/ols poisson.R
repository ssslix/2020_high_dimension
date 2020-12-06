#计算目标函数对c一阶偏导与二阶偏导的比值
l<-function(coef,X,y,n,EXY){
  f1=0
  f2=0
  for(i in 1:n){
    f1=f1-2*t(EXY)%*%X[,i]*(y[i]-exp(t(EXY)%*%X[,i]*coef))*exp(t(EXY)%*%X[,i]*coef)
    f2=f2-2*(t(EXY)%*%X[,i])^(2)*y[i]*exp(t(EXY)%*%X[,i]*coef)+4*(t(EXY)%*%X[,i])^(2)*exp(2*t(EXY)%*%X[,i]*coef)
  }
  return(f1/f2)
}
#牛顿迭代主体部分
NEWTON<-function(coef,X,y,n,EXY){
  eta=1e-10
  for(i in 1:100){
    a=coef
    b=coef-l(coef,X,y,n,EXY)
    if(abs(a-b)<=eta){
      break;
    }
    else{
      a=b
      coef=b
    }
  }
  return(c(coef))
}

#初始化beta，设定样本维度，大小，并给一个c的初值进行迭代
object<-function(p,n,coef){
  beta = rep(1,p)
  X=matrix(nrow=p,ncol=n)
  y = rep(0,n)
  for(j in 1:n)
  {
    X[,j]=runif(p,min=-1,max=1)
    lambda=exp(beta%*%X[,j])
    y[j]=rpois(1,lambda)+rnorm(1,0,1)
  }
  EXY=X%*%y/n
  coef=coef
  C=NEWTON(coef,X,y,n,EXY)
  EXY*C
  return(EXY*C)
}

object(10,2000,2)



experiment <- function(dim, n,init,times){
  true <- rep(1,dim)
  mse <- 0
  for (i in 1:times){
    beta <- object(dim,n,init)
    mse <- mse + sqrt(sum((beta - true)^2))
  }
  return(mse)
}

experiment(10,2000,2,100)
