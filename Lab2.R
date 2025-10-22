# PREDICTIVE MODELING
# LAB 2



## L2.2
# a:

# Joint pdf
f <- function(x1, x2) {
  exp(-(x1 + x2)) * (x1 > 0 & x2 > 0)
}

# For viridis() palette
library(viridis)

# Plot joint pdf
g <- seq(-1, 3, l = 200)
lev <- seq(0, 1, l = 21)
col <- viridis(length(lev) - 1)
f_g <- outer(g, g, FUN = f)
image(g, g, f_g, col = col, breaks = lev, xlab = "x1", ylab = "x2")
contour(g, g, f_g, add = TRUE, levels = lev)

# b:

# Joint cdf
F <- function(x1, x2) {
  (1 - exp(-x1)) * (1 - exp(-x2)) * (x1 > 0 & x2 > 0)
}

# Plot joint cdf
lev <- seq(0, 1, l = 21)
F_g <- outer(g, g, FUN = F)
image(g, g, F_g, col = col, breaks = lev, xlab = "x1", ylab = "x2")
contour(g, g, F_g, add = TRUE, levels = lev)

# c:

# Marginal pdf and cdf
f_X1 <- function(x1) {
  exp(-x1) * (x1 > 0)
}
F_X1 <- function(x1) {
  (1 - exp(-x1)) * (x1 > 0)
}

# Plot marginal pdf and cdf
plot(g, f_X1(g), xlab = "x1", ylab = "Pdf/cdf", type = "l", col = 2,
     ylim = c(0, 1.25), main = "Marginals of X1")
lines(g, F_X1(g), col = 3)
legend("topleft", legend = c("Pdf", "cdf"), col = 2:3, lwd = 2)

# d:

# Conditional pdf
f_X1_X2 <- function(x1, x2) {
  exp(-x1) * (x1 > 0)
}

# Plot conditional pdf dynamically (if it depended on x2, the plot would update)
library(manipulate)
manipulate({
  
  plot(g, f_X1_X2(x1 = g, x2 = x2), type = "l", ylim = c(0, 1.25),
       main = paste("Pdf of X1 | X2 =", x2), xlab = "x1", ylab = "Density")
  
}, x2 = slider(0, 3, step = 0.01))

## L2.4---------------------------------------------------------------------------

# a:

# Joint pdf
f <- function(x1, x2) {
  2 * (0 < x1 & x1 < x2 & x2 < 1)
}

# Plot joint pdf
g <- seq(-1, 2, l = 200)
lev <- seq(0, 2, l = 21)
col <- viridis(length(lev) - 1)
f_g <- outer(g, g, FUN = f)
image(g, g, f_g, col = col, breaks = lev, xlab = "x1", ylab = "x2")
contour(g, g, f_g, add = TRUE, levels = lev)

# b:

# Joint cdf
F <- function(x1, x2) {
  0 * (x1 < 0 | x2 < 0) +
    (2 * x1 * x2 - x1^2) * (0 <= x1 & x1 < x2 & x2 < 1) +
    x2^2 * (0 <= x2 & x2 <= 1 & x1 >= x2) +
    (2 * x1 - x1^2) * (0 <= x1 & x1 <= 1 & x2 > 1) +
    1 * (x1 > 1 & x2 > 1)
}

# Plot joint cdf
F_g <- outer(g, g, FUN = F)
lev <- seq(0, 1, l = 21)
image(g, g, F_g, col = col, breaks = lev, xlab = "x1", ylab = "x2")
contour(g, g, F_g, add = TRUE, levels = lev)

# c:

# Marginal pdf and cdf
f_X1 <- function(x1) {
  2 * (1 - x1) * (0 < x1 & x1 < 1)
}
F_X1 <- function(x1) {
  x1 * (2 - x1) * (0 < x1 & x1 < 1) + (x1 >= 1)
}

# Plot marginal pdf and cdf
plot(g, f_X1(g), xlab = "x1", ylab = "Pdf/cdf", type = "l", col = 2,
     ylim = c(0, 2), main = "Marginals of X1")
lines(g, F_X1(g), col = 3)
legend("topleft", legend = c("Pdf", "cdf"), col = 2:3, lwd = 2)

# d:

# Marginal pdf and cdf
f_X2 <- function(x2) {
  2 * x2 * (0 < x2 & x2 < 1)
}
F_X2 <- function(x2) {
  x2^2 * (0 < x2 & x2 < 1) + (x2 >= 1)
}

# Plot marginal pdf and cdf
plot(g, f_X2(g), xlab = "x2", ylab = "Pdf/cdf", type = "l", col = 2,
     ylim = c(0, 2), main = "Marginals of X2")
lines(g, F_X2(g), col = 3)
legend("topleft", legend = c("Pdf", "cdf"), col = 2:3, lwd = 2)

# e:

# Conditional pdf of X1 given X2 = x2
f_X1_X2 <- function(x1, x2) {
  1 / x2 * (0 < x1 & x1 < x2 & x2 < 1)
}

# Plot conditional pdf dynamically
manipulate({
  
  plot(g, f_X1_X2(x1 = g, x2 = x2), type = "l", ylim = c(0, 2),
       main = paste("Pdf of X1 | X2 =", x2), xlab = "x1", ylab = "Density")
  
}, x2 = slider(-1, 2, step = 0.01))

# f:

# Conditional pdf of X2 given X1 = x1
f_X2_X1 <- function(x1, x2) {
  1 / (1 - x1) * (0 < x1 & x1 < x2 & x2 < 1)
}

# Plot conditional pdf dynamically
manipulate({
  
  plot(g, f_X2_X1(x1 = x1, x2 = g), type = "l", ylim = c(0, 2),
       main = paste("Pdf of X2 | X1 =", x1), xlab = "x2", ylab = "Density")
  
}, x1 = slider(-1, 2, step = 0.01))

# g:

# Conditional expectations
E_X1_X2 <- function(x2) {
  x2 / 2
}
E_X2_X1 <- function(x1) {
  (1 + x1) / 2
}

# Plot joint pdf and conditional expectations
g <- seq(-1, 2, l = 200)
lev <- seq(0, 2, l = 21)
col <- viridis(length(lev) - 1)
f_g <- outer(g, g, FUN = f)
image(g, g, f_g, col = col, breaks = lev, xlab = "x1", ylab = "x2")
contour(g, g, f_g, add = TRUE, levels = lev)
lines(g, E_X2_X1(x1 = g), col = 2, lwd = 2)
lines(E_X1_X2(x2 = g), g, col = 4, lwd = 2)
legend("topleft",legend = c("E[X2 | X1 = x1]", "E[X1 | X2 = x2]"),
       lwd = 2, col = c(2, 4))


