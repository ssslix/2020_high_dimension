p=2
m=10
n=1000
#object<-function(p,m,n){
beta  <- diag(rep(1,m))[,1:p]##是个单位阵，取前p个行向量
gamma<-rep(1,p)
X=matrix(nrow=m,ncol=n)
y = rep(0,n)
for(j in 1:n)
{
  X[,j]=rnorm(m,0,1)
  lambda1=sum(t(beta)%*%X[,j])
  lambda2=exp(sum(t(beta)%*%X[,j]))
  lambda3=exp(sum(t(beta)%*%X[,j]))/(1+exp(sum(t(beta)%*%X[,j])))
  #y[j]=rnorm(1,lambda1,1)+rnorm(1,1,0.1) #正态
  #y[j]=rpois(1,lambda2)+rnorm(1,0,0.1)   #泊松
  #y[j]=rbinom(1,1,lambda3)++rnorm(1,0,0.1)  #logistic
  #y[j]=cos(2*c(t(beta[,1])%*%X[,j]))+cos(c(t(beta[,2])%*%X[,j]))+rnorm(1,1,0.1)
  #y[j]=rnorm(1,2*X[2,j],1)+rnorm(1,X[1,j],1)+rnorm(1,0,0.1) 
  #y[j]=exp(X[1,j])+exp(X[2,j])+rnorm(1,0,0.1) 
  y[j]=rpois(1,exp(X[1,j]))+rpois(1,exp(X[2,j]))+rnorm(1,0,0.1)
  y[j]=2*X[2,j]+X[1,j]+rnorm(1,0,0.1) 
}
mu=rowMeans(X)
var=cov(t(X-mu))
lamda =solve(eigen(var)$vectors)%*%var%*%(eigen(var)$vectors)
lamda_sqrt <- diag(sqrt(diag(lamda)))
var_sqrt <- (eigen(var)$vectors)%*%lamda_sqrt%*%solve(eigen(var)$vectors)
sd=solve(var_sqrt)
Z=sd%*%(X-mu)
y=y-mean(y)
y=c(y)
H1=matrix(0,nrow=m,ncol=m)
H2=matrix(0,nrow=m,ncol=m)
a=Z%*%y/n
H1=t(t(Z)*(y-mean(y)))%*%t(Z)/n
H2=t(t(Z)*(y-mean(y)-c(t(a)%*%Z)))%*%t(Z)/n
alpha1=eigen(H1)$vectors[,order(abs(eigen(H1)$values),decreasing = T)]
 alpha2=eigen(H2)$vectors[,order(abs(eigen(H2)$values),decreasing = T)]
c1=sd%*%alpha1[,1:p]
c2=sd%*%alpha2[,1:p]
ZZ1=t(c1)%*%X
ZZ2=t(c2)%*%X

#object(1,5,1000)