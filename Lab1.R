# PREDICTIVE MODELING
# LAB 1 


## L1.1
# a:

(exp(2) + sin(2)) / (2 + acos(0.5))

# b:

sqrt(3^2.5 + log(10))

# c:

(2^0.93 - log2(3 + sqrt(2 + sin(1)))) *
  (10^(tan(1 / 3))) * sqrt(3^2.5 + log(10))

## L1.2
# a:

y <- (-123)

# b:

z <- log(y^2)

# c:

y <- (y - z) / (y + z^2)
rm(z)

# d:

y

## L1.3
# a:

x <- c(1, 7, 3, 4)
x

# b:

y <- 100:1
y

# c:

e <- 2:5
z <- c(2^e, 2^5 * 3)
z

# d:

x[2] + y[4]
cos(x[3]) + sin(x[2]) * exp(-y[2])

# e:

x[2] <- 0
y[2] <- -1
y
x[2] + y[4]
cos(x[3]) + sin(x[2]) * exp(-y[2])

# f:

z <- y[x + 1]
z

## L1.4
# a:

mean(y)
median(y)
var(y)

# b:

mean(y + 1)
median(y + 1)
var(y + 1)

# c:

max(y)
which.max(y)

# d:

y.s <- sort(y)
y.s[c(5, 76)]

# e:

cov(y,y)
var(y)

## L1.5
# a:

M <- matrix(c(y[3:5], y[3:5]^2, log(y[3:5])), nrow = 3, byrow = TRUE)
M

# b:

my_df <- data.frame(y = y[3:5], y2 = y[3:5]^2, logy = log(y[3:5]))
my_df

# c:

l <- list(x = x, M = M)
l$x
l$M

# d:

my_df_2 <- my_df^2
my_df_2

# e:

log(my_df + my_df2)

## L1.6
# a:

data(faithful)
help(faithful)

# b:

dim(faithful)
head(faithful)

# c:

faithful$eruptions[5]
faithful[5, 1]

# d:

summary(faithful$waiting)

## L1.7
# a:

x <- 0.3 * (1:4)
x

# b:

seq(0, 1, length = 100)

# c:

x <- c(0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0)
table(x)
table(rev(x))

# d:

seq(0.1, 100.1, by = 1)
seq(0, 100, by = 1) + 0.1
0:100 + 0.1
rev(100:1) + 0.1

## L1.8
# a:

head(faithful)
plot(faithful)

# b:

abline(a = 110, b = -15, col = "red")

# c:

x <- seq(-pi, pi, l = 500)
y <- sin(x)
plot(x, y)
lines(x, sin(2 * x), col = "red")
lines(x, sin(3 * x), col = "blue")
lines(x, sin(4 * x), col = "orange")

## L1.9
# a:

qf(0.9, df1 = 1, df2 = 5)
qf(0.95, df1 = 1, df2 = 5)
qf(0.99, df1 = 1, df2 = 5)

# b:

x <- seq(0, 1, length = 100)
plot(x, punif(x))

# c:

rpois(100, lambda = 5)

# d:

u <- runif(100, min = -1, max = 1)
u
mean(u)

# e:

x <- seq(-4, 4, l = 200)
plot(x, dt(x, df = 1), type = "l", ylim = c(0, 0.4))
lines(x, dt(x, df = 5), col = 2)
lines(x, dt(x, df = 10), col = 3)
lines(x, dt(x, df = 50), col = 4)
lines(x, dt(x, df = 100), col = 5)

## L1.10

## L1.11

# b. Create a function that samples a N(0, 1) and returns the first sampled point that is larger than 4.

b_fun <- function() {
  
  repeat{
    
    x <- rnorm(n = 1, mean = 0, sd = 1)
    if (x > 4) break
    
  }
  return(x)
  
}

b_fun()

# c. Create a function that simulates N samples from the distribution of max(X1, . . . , Xn) where X1, . . . , Xn are iid U(0, 1).

c_fun <-function(N, n) {
  
  maxes <- numeric(N)
  for (i in 1:N) {
    
    # Simulate X1, . . . , Xn are iid U(0, 1)
    x <- runif(n = n, min = 0, max = 1)
    
    # max(X1, . . . , Xn)
    maxes[i] <- max(x)
    
  }
  return(maxes)
  
}

c_fun(N = 5, n = 3)

## L1.12
# a:

# Clear workspace
ls()
rm(list = ls())
ls()

# Read data
credit <- read.table("datasets/creditab.txt", header = TRUE)
head(credit)

# b:

credit[1:10, ]
head(credit, n = 10)

# c:

creditnew <- credit
names(creditnew) <- c("cred", "current", "duration", "amount", "former", "purpose")
write.table(x = creditnew, file = "datasets/creditnew.txt",
            row.names = FALSE)
creditnew <- read.table("datasets/creditnew.txt", header = TRUE)
head(creditnew)

# d:

dim(creditnew)
n <- nrow(creditnew)
n

# e:

# cred
attach(creditnew)
cred

# f:

is.factor(cred)
is.vector(duration)

# g:

table(cred)

# h:

table(cred) / sum(table(cred))
prop.table(table(cred))

# i:

summary(duration)

# j:

summary(cred)
summary(as.factor(cred))

# k:

p <- seq(0, 1, by = 0.1)
q <- quantile(duration, p)
q

# l:

par(mfrow = c(1, 2))
hist(duration)
hist(amount, main = "Histogram of amount", xlab = "Amount", prob = TRUE)
par(mfrow = c(1, 1)) # Restore graphical device, a good practice

# m:

boxplot(amount)

# n:

pie(table(cred))

# o:

barplot(table(cred))

# p:

barplot(100 * prop.table(table(cred)))

## L1.13
# a:

# install.packages("sae")
library(sae)
data(cornsoybean)
?cornsoybean
head(cornsoybean)

# b:

n <- dim(cornsoybean)[1]
n

# c:

table(cornsoybean$County)

# d:

summary(cornsoybean[, -1])

# e:

par(mfrow = c(1, 2))
plot(cornsoybean$CornPix, cornsoybean$CornHec,
     xlab = "Corn Pixels", ylab = "Corn Hectares")
plot(cornsoybean$SoyBeansPix, cornsoybean$CornHec,
     xlab = "Soybeans Pixels", ylab = "Corn Hectares")
par(mfrow = c(1, 1)) # Restore graphical device, a good practice

# f:

plot(cornsoybean$SoyBeansPix, cornsoybean$SoyBeansHec,
     xlab = "Soybeans Pixels", ylab = "Soybeans Hectares")
plot(cornsoybean$CornPix, cornsoybean$SoyBeansHec,
     xlab = "Corn Pixels", ylab = "Soybeans Hectares")





# L1.1.14
# for the simple integrals
fa = function(x){
  return(cos(x))
}

integral = function(a, b){
  steps = seq(a, b, (b-a)/1000)
  sum = 0
  for (i in steps){
    sum = sum + fa(i)
  }
  
  return (sum*(b-a)/1000)
}

a = integral(1,5)
sin(5) - sin(1)

# for the double integrals
fb = function(a, b){
  return(sin(a*b))
}


integral_double = function(a, b, c, d){
  steps1 = seq(a, b, (b-a)/100)
  steps2 = seq(c, d, (d-c)/100)
  grid = expand.grid(steps1, steps2)
  sum = 0
  for(i in grid){
    sum = sum + fb(i[1], i[2])    
  }
  return(sum*(b-a)*(d-c)/10000)
}

integral_double(-1, 1, -1, 1)

# special integral


# L1.1.15



# L1.1.16


 