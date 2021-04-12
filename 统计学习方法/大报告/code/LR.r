library(MASS)
fix(Boston)



lm.fit <- lm(medv~., data = Boston)
summary(lm.fit)
