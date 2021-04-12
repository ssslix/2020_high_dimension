library(MASS)
library(glmnet)
x <- model.matrix(medv~.,Boston)[,-1]
y <- Boston$medv

grid=10^seq(-3,10 ,length=100)
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)
coef(ridge.mod)[,100]
coef(ridge.mod)[,1]
plot(ridge.mod)


set.seed(512)
cv.out=cv.glmnet(x,y,alpha=0)#交叉验证
plot(cv.out)
bestlam1=cv.out$lambda.min
bestlam2=cv.out$lambda.1se

ridge.pred=predict(ridge.mod,s=c(bestlam1,bestlam2),newx=x)
mean((ridge.pred[,1]-y)^2)#当lambda=lambda.min
mean((ridge.pred[,2]-y)^2)#当lambda=lambda.lse
predict(ridge.mod,s=c(bestlam1,bestlam2),newx=x,type="coefficients") 
a<-predict(ridge.mod,s=c(bestlam1,bestlam2),newx=x,type="coefficients") 
