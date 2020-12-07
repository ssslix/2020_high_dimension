l<-function(coef,X,y,n,EXY){
  f1=0
  f2=0
  for(i in 1:n){
   
    
    
    f1 = f1 +2*t(EXY)%*%X[,i]*exp(-t(EXY)%*%X[,i]*coef)/(1+exp(-t(EXY)%*%X[,i]*coef))^2*
      (y[i]-1/(1+exp(-t(EXY)%*%X[,i]*coef)))
    f2 = f2 + 2*(-(t(EXY)%*%X[,i]*exp(-t(EXY)%*%X[,i]*coef)/(1+exp(-t(EXY)%*%X[,i]*coef))^2)^2+
                  (y[i]-1/(1+exp(-t(EXY)%*%X[,i]*coef)))*(((t(EXY)%*%X[,i])^2*exp(-t(EXY)%*%X[,i]*coef)*
                                                             (1-exp(-t(EXY)%*%X[,i]*coef))/(1+exp(-t(EXY)%*%X[,i]*coef))^3)))
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
    y[j]=rbinom(1,1,lambda)+rnorm(1,0,1)
  }
  EXY=X%*%y/n
  coef=coef
  C=NEWTON(coef,X,y,n,EXY)
  EXY*C
  return(EXY*C)
}

object(10,1000,2)
  

experiment <- function(dim, n,init,times){
  true <- rep(1,dim)
  mse <- 0
  for (i in 1:times){
    beta <- object(dim,n,init)
    mse <- mse + sqrt(sum((beta - true)^2))
  }
  return(mse)
}

# experiment(20,2000,2,100)


# dim:10 <- mse*100 = 33.25142
# dim:15 <- mse*100 =  128.2584
# dim:20 <- mse*100 =  169.1181