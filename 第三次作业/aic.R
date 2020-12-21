aicc<-function(data,y,p){
  q=ncol(data)
  n=nrow(data)
  fit=lm(y~.,data=data)
  err=fit$residuals
  aic=mean(err^2) +q*mean(ERR^2)* 2/(n-p) 
  return(aic)
}
t1=proc.time()
library(MASS)
p=6
n=1000
beta = c(0,1,0,0,1,1)
X=matrix(nrow=p,ncol=n)
y = rep(0,n)
for(j in 1:n)
{
  X[,j]=runif(p,min=-1,max=1)
  lambda=exp(beta%*%X[,j])
  y[j]=rnorm(1,lambda,1)+rnorm(1,0,1)
}
shuju<-data.frame(t(rbind(y,X)))
fit=lm(y~.,data=shuju)
ERR=fit$residuals
#step(lm(y~.,data=shuju))
BitMatrix <- function(n){
  # 作用：返还Bit矩阵
  # Args:n：数据长度
  
  set <- 0:(2^n-1)
  rst <- matrix(0,ncol = n,nrow = 2^n)
  for (i in 1:n){
    rst[, i] = ifelse((set-rowSums(rst*rep(c(2^((n-1):0)), each=2^n)))/(2^(n-i))>=1, 1, 0)
  }
  rst
}
bit<-BitMatrix(p)
bit<-cbind(rep(1,2^p),bit)
aic<-rep(0,2^p-1)
for(i in 1:2^p)
{
  a<-(shuju[,bit[i,]!=0])
  if(i!=1){
    aic[i-1]<-aicc(a,y,p)
  }
}
aorder=order(aic)[1]+1
bit[aorder,2:(p+1)]
t2=proc.time()
t=t2-t1
print(paste0('执行时间：',t[3][[1]],'秒'))
