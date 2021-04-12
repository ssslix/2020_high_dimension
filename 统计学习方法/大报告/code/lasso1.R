set.seed(123)
library(glmnet)
library(msgps)
library(lars)
library(MASS)
X <- model.matrix(medv~.,Boston)[,-1]
y <- Boston$medv
n <- dim(X)[1];p <- (dim(X)[2]+1)


y1 = rep(0,n)


##读取实际数据 糖尿病代码 X10维，p 1维，清空工作区，注释以上代码并注释最后第三行print（beta）
# data(diabetes)
# 
# attach(diabetes)
# X=t(x)
# n=ncol(X)
# p=10
con <- rep(1,n)
X1 = rbind(con,t(X))

#给出初始wi，转换X数据

##对转换后的数据做LASSO
##设超参数lambda
lambda=0.02118502
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
      y1=as.vector(y-beta_hat[-i]%*%X1[-i,])
      X1[i,]=matrix(c(X1[i,]),nrow=1,ncol=n)
      beta_i_ols=solve(sum(X1[i,]^2))%*%X1[i,]%*%y1
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
beta_hat_1
res<-cv.glmnet(x=X,y=y,family="gaussian",alpha = 1,nfolds = 10)
# setEPS()
# postscript("lassocv.eps")
plot(res)
# dev.off()
# lar1 <- lars(t(X),y,type = "lasso")
# beta_hat_1
# beta_hat_2 <- lar1$beta
# summary(lar1)
# beta_hat_2
rmin <- glmnet(x=X,y=y,family = "gaussian",alpha = 1,lambda = res$lambda.min)
coef(rmin)
rmin.fit <-  t(t(coef(rmin))%*%X1)
rmin.mse <-  mean((rmin.fit-y)^2)
ofit <- t(t(beta_hat_1)%*%X1)
omse <- mean((ofit-y)^2)
