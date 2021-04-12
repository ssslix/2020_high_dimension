library(MASS)
library(xtable)
X <- model.matrix(medv~.,Boston)[,-1]
y <- Boston$medv
n <- dim(X)[1];p <- dim(X)[2]
alpha=0.01
X <- cbind(rep(1,n),X)
X <- t(X)



coef=solve(X%*%t(X)+alpha*diag(p+1))%*%X%*%y
coef

lr <- lm.ridge(medv~.,data= Boston,lambda = 0.01,model=TRUE)
lr

