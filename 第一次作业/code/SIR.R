library(ggplot2)
object <- function(p, q, n) {
  #p初始化x维数,q初始化t(beta)作用在X上用的有效维度 n样本量
  
  h = ceiling(sqrt(n))
  #初始化beta矩阵m*p维
  beta  <- diag(rep(1, p))[, 1:q]##是个单位阵，取前q个列向量
  X = matrix(nrow = p, ncol = n)
  y = rep(0, n)
  for (j in 1:n)
  {
    #按标准正态分布生成每个样本点Xi的每个维度的值 并按不同模型生成y
    X[, j] = rnorm(p, 0, 1)
    lambda1 = sum(t(beta) %*% X[, j])
    lambda2 = exp(sum(t(beta) %*% X[, j]))
    lambda3 = exp(sum(t(beta) %*% X[, j])) / (1 + exp(sum(t(beta) %*% X[, j])))
    #y[j] = rnorm(1, lambda1, 1) + rnorm(1, 0, 0.1) #正态
    #y[j]=rpois(1,lambda2)+rnorm(1,0,0.1)   #泊松
    #y[j]=rbinom(1,1,lambda3)++rnorm(1,0,0.1)  #logistic
     y[j] = sin(2 * c(t(beta[,1]) %*% X[, j]))  + sin(c(t(beta[, 2]) %*% X[, j])) 
    + rnorm(1, 0, 0.1)
    #y[j]=rnorm(1,2*X[2,j],1)+rnorm(1,X[1,j],1)+rnorm(1,0,0.1)
    #y[j]=exp(X[1,j])+exp(X[2,j])+rnorm(1,0,0.1)
    #y[j]=rpois(1,exp(X[1,j]))+rpois(1,exp(X[2,j]))+rnorm(1,0,0.1)
    #y[j]=2*X[2,j]+X[1,j]
    #y[j]=exp(X[1,j])/(1+exp(X[1,j]))+exp(X[2,j])/(1+exp(X[2,j]))
  }
  #标准化X
  mu = rowMeans(X)
  var = cov(t(X - mu))
  lamda = solve(eigen(var)$vectors) %*% var %*% (eigen(var)$vectors)
  lamda_sqrt <- diag(sqrt(diag(lamda)))
  var_sqrt <-
    (eigen(var)$vectors) %*% lamda_sqrt %*% solve(eigen(var)$vectors)
  sd = solve(var_sqrt)
  Z = sd %*% (X - mu)
  #y从小到大排序，按该次序对Z进行重排
  Zord = Z[, order(y)]
  y1 = rep(0, n)
  pn = floor(n / h)
  y1[1:pn] = 1
  for (i in 1:(h - 2)) {
    y1[(pn * i + 1):(pn * (i + 1))] = i + 1
  }
  y1[(pn * (h - 1) + 1):n] = h
  #计算M
  g = matrix(nrow = p, ncol = h)
  prob = rep(0, h)
  for (i in 1:h) {
    prob[i] = length(y1[y1 == i]) / n
    g[, i] = apply(Zord[, y1 == i], 1, mean)
  }
  M = g %*% diag(prob) %*% t(g)
  alpha = eigen(M)$vector
  beta_hat = sd %*% alpha[, 1:q]
  
  
  ###计算 估计的beta和真实beta张成空间的距离
  ##把beta和beta_hat施密特正交化
  
  span_beta = qr.Q(qr(beta))
  span_beta_hat = qr.Q(qr(beta_hat))
  # ##计算距离
  dis = q
  for (i in 1:q) {
    dis = dis - sum((t(span_beta[, i]) %*% span_beta_hat) ^ 2)
  }
  dis = sqrt(dis)
  return(dis)
}
##估计原本beta对应的q的值
"
fixq<-function(eigenvalue,n,p){
  cn=1/n^(1/3)
  t=rep(0,p-1)
  q_hat=0
  for(i in 1:(p-1)){
    if((eigenvalue[i+1]+cn)/(eigenvalue[i]+cn)<=0.5)
    {
      q_hat=i
    }
  }
  return(q_hat)
}
q_hat=fixq(eigen(M)$value,n,p)

"
# object(10,2,1000)
experiment <-function(dim_x,dim_r,n){
  dis <- rep(0,11)
  for (i in 1:11){
    for(j in 1:100){
      dis[i] <-dis[i]+  object(i+9,dim_r,n)/100
    }
    
  }
  
  dis <- data.frame(dim,dis)
  return(dis)
}
dim <- rep(10:20)
dis <- experiment(20,2,1000)
#dis_cos <- experiment(20,1,1000)
#dis_sin <- experiment(20,1,1000)
#dim <- rep(10:20)
#dis <- data.frame(dim,dis_sin,dis_cos)
# ggplot(data=dis)+geom_point(aes(dim,dis_cos,colour='cos'))+geom_point(aes(dim,dis_sin,colour="sin"))

ggplot(data=dis,aes(x=dim,y=dis))+geom_line(linetype="dotted")+geom_point(size=4, shape=20)
