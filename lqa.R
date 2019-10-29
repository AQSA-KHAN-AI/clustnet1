library(lqa)


set.seed (1111)

n <- 200
p <- 5
X <- matrix (rnorm (n * p), ncol = p)
X[,2] <- X[,1] + rnorm (n, sd = 0.1)
X[,3] <- X[,1] + rnorm (n, sd = 0.1)
true.beta <- c (1, 2, 0, 0, -1)
y <- drop (X %*% true.beta) + rnorm (n)

obj1 <- lqa (y ~ X, family = gaussian (), penalty = lasso (1.5), 
             control = lqa.control ())
obj1$coef


set.seed (4321)

n <- 25
p <- 5
X <- matrix (rnorm (n * p), ncol = p)
X[,2] <- X[,1] + rnorm (n, sd = 0.1)
X[,3] <- X[,1] + rnorm (n, sd = 0.1)
true.beta <- c (1, 2, 0, 0, -1)

family1 <- binomial ()
eta.true <- drop (X %*% true.beta)
mu.true <- family1$linkinv (eta.true)
prob1 <- sum (as.integer (y > 0)) / n
nvec <- 1 : n
y2 <- sapply (mu.true, function (n.vec) {rbinom (1, 1, mu.true)})

obj2 <- lqa (y2 ~ X, family = binomial (), 
             penalty = fused.lasso (c (0.0001, 0.2)))
obj2$coef
