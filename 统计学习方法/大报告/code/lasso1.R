library(msgps)
library(lars)
p=10
n=100
X=matrix(nrow=p,ncol=n)#初始化样本

y = rep(0,n)
y1 = rep(0,n)
beta =  rbinom(p,1,0.4)*rbinom(p,10,0.7)#初始化真实参数beta
for(j in 1:n)#按正态分布生成数据y
{
  X[,j]=rnorm(p,0,1)
  lambda=beta%*%X[,j]
  y[j]=beta%*%X[,j]+rnorm(1,0,0.1)
}
##读取实际数据 糖尿病代码 X10维，p 1维，清空工作区，注释以上代码并注释最后第三行print（beta）
# data(diabetes)
# 
# attach(diabetes)
# X=t(x)
# n=ncol(X)
# p=10

X_new=matrix(nrow=p,ncol=n)
y=y-mean(y)
X=X-rowMeans(X)
coef=solve(X%*%t(X))%*%X%*%y
#给出初始wi，转换X数据
##给出超参数gamma=-1,-0.5,-2
gamma=1

##对转换后的数据做LASSO
##设超参数lambda
lambda=0.01
##采用循环坐标下降法得到估计
cd<-function(){
  #随机给出初值betahat,以及用于迭代的暂时储存估计beta的变量
  beta_hat =rnorm(p,0,1)
  temp1= beta_hat
  temp2=temp1-1
  ##设置循环截止条件
  while (sum(abs(temp1-temp2))>1e-12) {
    temp1= beta_hat
    for(i in 1:p){
      y1=as.vector(y-beta_hat[-i]%*%X[-i,])
      X[i,]=matrix(c(X[i,]),nrow=1,ncol=n)
      beta_i_ols=solve(sum(X[i,]^2))%*%X[i,]%*%y1
      if(abs(beta_i_ols)>lambda){
        beta_hat[i]=sign(beta_i_ols)*(abs(beta_i_ols)-lambda)
      }
      else{
        beta_hat[i]=0
      }
    }
    temp2=beta_hat
  }
  
  return( beta_hat)
}
beta_hat_1=cd()
##R自带函数结果
alasso<-msgps(t(X),y,penalty="alasso",gamma=1,lambda=lambda,intercept=FALSE)
fit<-summary(alasso);
beta_hat_2<-fit$dfcp_result$coef
##结果比较
names(beta) <-  row.names(X)
names(beta_hat_1) <-  row.names(X)
names(beta_hat_2) <-  row.names(X)
print(beta)
print(beta_hat_1)
print(as.vector(beta_hat_2))