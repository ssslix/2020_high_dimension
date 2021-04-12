
p=3
n=100
alpha=0.01
beta = rep(1,p)#初始化beta 设为m维全为1
X=matrix(nrow=p,ncol=n)#初始化样本
y = rep(0,n)
for(j in 1:n)#按正态分布生成数据y
{
  X[,j]=runif(p,min=-1,max=1)
  lambda=beta%*%X[,j]
  y[j]=rnorm(1,lambda,1)+rnorm(1,1,0.1)
}

coef=solve(X%*%t(X)+alpha*diag(p))%*%X%*%y
coef