library(MASS)
library(xtable)
fix(Boston)

head(Boston)

lm.fit <- lm(medv~., data = Boston)
summary(lm.fit)
