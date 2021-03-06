library(stargazer)
library(ggplot2)
object <- function(p, q, n) {
  #p初始化x维数,q初始化t(beta)作用在X上用的有效维度 n样本量
  
  
  
  #初始化beta矩阵p*q维
  beta  <- diag(rep(1, p))[, 1:q]##是个单位阵，取前p个行向量
  X = matrix(nrow = p, ncol = n)
  y = rep(0, n)
  for (j in 1:n)
  {
    #按标准正态分布生成每个样本点Xi的每个维度的值 并按不同模型生成y
    X[, j] = rnorm(p, 0, 1)
    lambda1 = sum(t(beta) %*% X[, j])
    lambda2 = exp(sum(t(beta) %*% X[, j]))
    lambda3 = exp(sum(t(beta) %*% X[, j])) / (1 + exp(sum(t(beta) %*% X[, j])))
    #y[j]=rnorm(1,lambda1,1)+rnorm(1,0,0.1) #正态
    #y[j]=rpois(1,lambda2)+rnorm(1,0,0.1)   #泊松
    #y[j]=rbinom(1,1,lambda3)++rnorm(1,0,0.1)  #logistic
     y[j] = cos(2 * c(t(beta[, 1]) %*% X[, j])) + sin(c(t(beta[, 2]) %*% X[, j])) +
      +rnorm(1, 1, 0.1)
    #y[j]=rnorm(1,2*X[2,j],1)+rnorm(1,X[1,j],1)+rnorm(1,0,0.1)
    #y[j]=exp(X[1,j])+exp(X[2,j])+rnorm(1,0,0.1)
    #y[j]=rpois(1,exp(X[1,j]))+rpois(1,exp(X[2,j]))+rnorm(1,0,0.1)
    #y[j]=2*X[2,j]+X[1,j]
    #y[j]=exp(X[1,j])/(1+exp(X[1,j]))+exp(X[2,j])/(1+exp(X[2,j]))
  }
  #中心化标准化
  mu = rowMeans(X)
  var = cov(t(X - mu))
  lamda = solve(eigen(var)$vectors) %*% var %*% (eigen(var)$vectors)
  lamda_sqrt <- diag(sqrt(diag(lamda)))
  var_sqrt <-
    (eigen(var)$vectors) %*% lamda_sqrt %*% solve(eigen(var)$vectors)
  sd = solve(var_sqrt)
  Z = sd %*% (X - mu)
  y = y - mean(y)
  y = c(y)
  #计算H1，H2
  H1 = matrix(0, nrow = p, ncol = p)
  H2 = matrix(0, nrow = p, ncol = p)
  a = Z %*% y / n
  H1 = t(t(Z) * (y - mean(y))) %*% t(Z) / n
  H2 = t(t(Z) * (y - mean(y) - c(t(a) %*% Z))) %*% t(Z) / n
  alpha1 = eigen(H1 %*% H1)$vectors
  alpha2 = eigen(H2 %*% H2)$vectors
  #得到beta矩阵估计值
  beta_hat1 = sd %*% alpha1[, 1:q]
  beta_hat2 = sd %*% alpha2[, 1:q]
  
  ZZ1 = t(beta_hat1) %*% X
  ZZ2 = t(beta_hat2) %*% X
  
  ###计算 估计的beta和真实beta张成空间的距离
  ##把beta和beta_hat施密特正交化
  span_beta = qr.Q(qr(beta))
  span_beta_hat1 = qr.Q(qr(beta_hat1))
  ##计算距离
  dis = q
  for (i in 1:q) {
    dis = dis - sum((t(span_beta[, i]) %*% span_beta_hat1) ^ 2)
  }
  dis = sqrt(dis)
  return( beta_hat1)
}
object(10,2,1000)

experiment <-function(dim_x,dim_r,n){
  dis <- rep(0,dim_x-dim_r)
  for (i in 1:(dim_x-dim_r)){
    for(j in 1:100){
      dis[i] <-dis[i]+ object(i+dim_r,dim_r,n)/100
    }
  }
  dim <- rep((dim_r+1):dim_x)
  dis <- data.frame(dim,dis)
  return(dis)
}

dis <- experiment(20,2,1000)


ggplot(data=dis,aes(x=dim,y=dis))+geom_line(linetype="dotted")+geom_point(size=4, shape=20)
