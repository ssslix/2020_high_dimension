#object<-function(p,m,n){
#m初始化x维数,p初始化t(beta)作用在X上用的有效维度 n样本量
p=2
m=4
n=1000
h=ceiling(sqrt(n))
#初始化beta矩阵m*p维
beta  <- diag(rep(1,m))[,1:p]##是个单位阵，取前p个行向量
X=matrix(nrow=m,ncol=n)
y = rep(0,n)
for(j in 1:n)
{
  #按标准正态分布生成每个样本点Xi的每个维度的值 并按不同模型生成y
  X[,j]=rnorm(m,0,1)
  lambda1=sum(t(beta)%*%X[,j])
  lambda2=exp(sum(t(beta)%*%X[,j]))
  lambda3=exp(sum(t(beta)%*%X[,j]))/(1+exp(sum(t(beta)%*%X[,j])))
  y[j]=rnorm(1,lambda1,1)+rnorm(1,1,0.1) #正态
  #y[j]=rpois(1,lambda2)+rnorm(1,0,0.1)   #泊松
  #y[j]=rbinom(1,1,lambda3)++rnorm(1,0,0.1)  #logistic
  y[j]=sin(2*c(t(beta[,1])%*%X[,j]))+sin(c(t(beta[,2])%*%X[,j]))+rnorm(1,1,0.1)
  #y[j]=rnorm(1,2*X[2,j],1)+rnorm(1,X[1,j],1)+rnorm(1,0,0.1) 
  #y[j]=exp(X[1,j])+exp(X[2,j])+rnorm(1,0,0.1) 
  #y[j]=rpois(1,exp(X[1,j]))+rpois(1,exp(X[2,j]))+rnorm(1,0,0.1)
  #y[j]=2*X[2,j]+X[1,j]
  #y[j]=exp(X[1,j])/(1+exp(X[1,j]))+exp(X[2,j])/(1+exp(X[2,j]))
}
#标准化X
mu=rowMeans(X)
var=cov(t(X-mu))
lamda =solve(eigen(var)$vectors)%*%var%*%(eigen(var)$vectors)
lamda_sqrt <- diag(sqrt(diag(lamda)))
var_sqrt <- (eigen(var)$vectors)%*%lamda_sqrt%*%solve(eigen(var)$vectors)
sd=solve(var_sqrt)
Z=sd%*%(X-mu)
#y从小到大排序，按该次序对Z进行重排
Zord=Z[,order(y)]
y1=rep(0,n)
pn=floor(n/h)
y1[1:pn]=1
for(i in 1:(h-2)){
  y1[(pn*i+1):(pn*(i+1))]=i+1
}
y1[(pn*(h-1)+1):n]=h
#计算M
g=matrix(nrow=m,ncol=h)
prob=rep(0,h)
for(i in 1:h){ 
  prob[i]=length(y1[y1==i])/n
 g[,i]=apply(Zord[,y1==i],1,mean)
}
M=g%*%diag(prob)%*%t(g)
alpha=eigen(M)$vector
beta_hat=sd%*%alpha[,1:p]
beta_hat
#得到beta矩阵估计值
# beta_hat1=sd%*%alpha1[,1:p]
# beta_hat2=sd%*%alpha2[,1:p]
# 
# ZZ1=t(beta_hat1)%*%X
# ZZ2=t(beta_hat2)%*%X

###计算 估计的beta和真实beta张成空间的距离
##把beta和beta_hat施密特正交化
# span_beta=qr.Q(qr(beta))
# span_beta_hat1=qr.Q(qr(beta_hat1))
# ##计算距离
# dis=p
# for(i in 1:p){
#   dis=dis-sum((t(span_beta[,i])%*%span_beta_hat1)^2)
# }
# dis=sqrt(dis)
# #return(dis)
# #}