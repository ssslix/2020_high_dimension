library(MASS)

x <- model.matrix(medv~.,Boston)[,-1]
y <- Boston$medv