# PREDICTIVE MODELING
# LAB 3

## L3.4
# Transformation
g <- function(x) -x * (x + 1) * (x < 0) + x * (x - 1) * (x > 0)
curve(g(x), from = -2, to = 2)

# Inverse of g1 and derivative of inverse of g1
g1_inv <- function(y) (y <= 1 / 4) * (-1 - sqrt(1 - 4 * y)) / 2
d_g1_inv <- function(y) (y <= 1 / 4) * 1 / sqrt(1 - 4 * y)



# Density of Y
f_Y <- function(y, f_X = function(x) dunif(x, -2, 2)) {
  
  d <- numeric(length(y))
  ind_1 <- y <= -1 / 4
  ind_2 <- y > -1 / 4 & y < 1/4
  ind_3 <- y >= 1 / 4
  y_1 <- y[ind_1]
  y_2 <- y[ind_2]
  y_3 <- y[ind_3]
  d[ind_1] <- # TODO in labs/homework
    d[ind_2] <- # TODO in labs/homework
    d[ind_3] <- # TODO in labs/homework
    return(d)
  
}

# L3.5

# Sample M samples of size n
mu <- 2
sigma <- 2
n <- 100
M <- 1e4
X <- matrix(rnorm(n * M, mean = mu, sd = sigma), nrow = n, ncol = M)

# a:

# Compute the statistics and obtain its histogram
Tn <- (colMeans(X) - mu) / (sigma / sqrt(n))
hist(Tn, probability = TRUE, breaks = 100)
rug(Tn)

# Overlay the distribution's pdf
curve(dnorm(x, mean = 0, sd = 1), col = 2, lwd = 2, add = TRUE)
# x <- seq(-4, 4, l = 200) # Another approach
# lines(x, dnorm(x, mean = 0, sd = 1), col = 2, lwd = 2)

# b:

