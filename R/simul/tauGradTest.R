z <- function(x, y)
{
  sin(x) + y^2
}

theta <- function(x, y) {
  x^2 + y
}

zthea <- function(theta, y)
{
  sin(sqrt(theta-y)) + y^2
}

n <- 100
x <- rnorm(n, 20)
y <- runif(n, 10, 50)
theta0 <- theta(x, y)




out0 <- rep(NA, n)
require("numDeriv")
for(i in 1:length(theta0))
{
  out0[i] <- grad(func = zthea,
                  x = theta0[i],
                  y = y[i])

}
