l<-function(coef,X,y,n,EXY){
  f1=0
  f2=0
  for(i in 1:n){
    f1=f1-2*t(EXY)%*%X[,i]*(y[i]-exp(t(EXY)%*%X[,i]*coef))*exp(t(EXY)%*%X[,i]*coef)
    f2=f2-2*(t(EXY)%*%X[,i])^(2)*y[i]*exp(t(EXY)%*%X[,i]*coef)+4*(t(EXY)%*%X[,i])^(2)*exp(2*t(EXY)%*%X[,i]*coef)
  }
  return(f1/f2)
}

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
set.seed(1)
object<-function(p,n,coef){
  beta = rep(1,p)
  X=matrix(nrow=p,ncol=n)
  y = rep(0,n)
  for(j in 1:n)
  {
    X[,j]=runif(p,min=-1,max=1)
    lambda=exp(beta%*%X[,j])/(1+exp(beta%*%X[,j]))
    y[j]=rbinom(1,lambda)+rnorm(1,0,1)
  }
  EXY=X%*%y/n
  coef=coef
  C=NEWTON(coef,X,y,n,EXY)
  EXY*C
  return(EXY*C)
}

object(5,100,2)
