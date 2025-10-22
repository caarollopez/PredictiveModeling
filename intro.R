
## INTRO R
##-------------------------------------------------------------------

## Simple computations

# The console can act as a simple calculator
1.0 + 1.1
2 * 2
3 / 2
2^3
1 / 0
0 / 0

# Use ";" for performing several operations in the same line
(1 + 3) * 2 - 1; 3 + 2

# Elemental mathematical functions
sqrt(2); 2^0.5
exp(1)
log(10) # Natural or neperian logarithm
log10(10); log2(10) # Logs in base 10 and 2
sin(pi); cos(0); asin(0)
tan(pi/3)
sqrt(-1)

# Remember to close the parenthesis -- errors below
# 1 +
# (1 + 3

## Variables and assignment

# Any operation that you perform in R can be stored in a variable
# (or object) with the assignment operator "<-"
x <- 1

# To see the value of a variable, simply type it
x

# A variable can be overwritten
x <- 1 + 1

# Now the value of x is 2 and not 1, as before
x

# Capitalization matters
X <- 3
x; X

# See what are the variables in the workspace
ls()

# Remove variables
rm(X)
# X # Error

## Vectors

# We combine numbers with the function "c"
c(1, 3)
c(1.5, 0, 5, -3.4)

# A handy way of creating integer sequences is the operator ":"
1:5
5:1

# Storing some vectors
myData <- c(1, 2)
myData2 <- c(-4.12, 0, 1.1, 1, 3, 4)
myData
myData2

# Entrywise operations
myData + 1
myData^2
exp(myData)

# If you want to access a position of a vector, use [position]
myData[1]
myData2[6]

# You can also change elements
myData[1] <- 0
myData

# Think on what you want to access...
myData2[7]
myData2[0]

# If you want to access all the elements except a position,
# use [-position]
myData2[-1]
myData2[-2]

# Also with vectors as indexes
myData2[1:2]
myData2[myData]

# And also
myData2[-c(1, 2)]

# But do not mix positive and negative indexes!
# myData2[c(-1, 2)] # Error

# Remove the first element
myData2 <- myData2[-1]

## Some functions

# Functions take arguments between parenthesis and transform them
# into an output
sum(myData)
prod(myData)

# Summary of an object
summary(myData)

# Length of the vector
length(myData)

# Mean, standard deviation, variance, covariance, correlation
mean(myData)
var(myData)
cov(myData, myData^2)
cor(myData, myData * 2)
quantile(myData)

# Maximum and minimum of vectors
min(myData)
which.min(myData)

# Usually the functions have several arguments, which are set
# by "argument = value" In the next case, the second argument is
# a logical flag to indicate the kind of sorting
sort(myData) # If nothing is specified, decreasing = FALSE is
# assumed
sort(myData, decreasing = TRUE)

# Do not know what are the arguments of a function? Use args
# and help!
args(mean)
?mean

## Matrices, data frames, and lists

# A matrix is an array of vectors
A <- matrix(1:4, nrow = 2, ncol = 2)
A

# Another matrix
B <- matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
B

# Matrix is a vector with dimension attributes
dim(A)
nrow(A)
ncol(A)

# Binding by rows or columns
rbind(1:3, 4:6)
cbind(1:3, 4:6)

# Entrywise operations
A + 1
A * B

# Accessing elements
A[2, 1] # Element (2, 1)
A[1, ] # First row -- this is a vector
A[, 2] # First column -- this is a vector

# Obtain rows and columns as matrices (and not as vectors)
A[1, , drop = FALSE]
A[, 2, drop = FALSE]

# Matrix transpose
t(A)

# Matrix multiplication
A %*% B
A %*% B[, 1]
A %*% B[1, ]

# Care is needed
# A %*% B[1, , drop = FALSE] # Incompatible product

# Compute inverses with "solve"
solve(A) %*% A

# A data frame is a matrix with column names
# Useful when you have multiple variables
myDf <- data.frame(var1 = 1:2, var2 = 3:4)
myDf

# You can change names
names(myDf) <- c("newname1", "newname2")
myDf

# You can access variables by its name with the "$" operator
myDf$newname1

# And create new variables also (they have to be of the same
# length as the rest of variables)
myDf$myNewVariable <- c(0, 1)
myDf

# Matrix to data frame
myDfA <- as.data.frame(A)
class(myDfA)

# Matrix to data frame
mymatDf <- as.matrix(myDf)
class(mymatDf)

# Structure
str(myDfA)
str(mymatDf)

# A list is a collection of arbitrary variables
myList <- list(myData = myData, A = A, myDf = myDf)

# Access elements by names
myList$myData
myList$A
myList$myDf

# Reveal the structure of an object
str(myList)
str(myDf)

# A less lengthy output
names(myList)

## More on data frames

# The iris dataset is already imported in R
# (beware: locfit has also an iris dataset, with different names
# and shorter)

# The beginning of the data
head(iris)

# "names" gives you the variables in the data frame
names(iris)

# So we can access variables by "$" or as in a matrix
iris$Sepal.Length[1:10]
iris[1:10, 1]
iris[3, 1]

# Information on the dimension of the data frame
dim(iris)

# "str" gives the structure of any object in R
str(iris)

# Recall the species variable: it is a categorical variable
# (or factor), not a numeric variable
iris$Species[1:10]

# Factors can only take certain values
levels(iris$Species)

# If a file contains a variable with character strings as
# observations (either encapsulated by quotation marks or not),
# the variable will become a factor when imported into R

## Vector-related functions

# The function "seq" creates sequences of numbers equally separated
seq(0, 1, by = 0.1)
seq(0, 1, length.out = 5)

# You can short the latter argument
seq(0, 1, l = 5)

# Repeat number
rep(0, 5)

# Reverse a vector
myVec <- c(1:5, -1:3)
rev(myVec)

# Another way
myVec[length(myVec):1]

# Count repetitions in your data
table(iris$Species)

## Logical conditions and subsetting

# Relational operators: x < y, x > y, x <= y, x >= y, x == y, x!= y
# They return TRUE or FALSE

# Smaller than
0 < 1

# Greater than
1 > 1

# Greater or equal to
1 >= 1 # Remember: ">="" and not "=>"" !

# Smaller or equal to
2 <= 1 # Remember: "<="" and not "=<"" !

# Equal
1 == 1 # Tests equality. Remember: "=="" and not "="" !

# Unequal
1 != 0 # Tests iequality

# TRUE is encoded as 1 and FALSE as 0
TRUE + 1
FALSE + 1

# In a vector-like fashion
x <- 1:5
y <- c(0, 3, 1, 5, 2)
x < y
x == y
x != y

# Subsetting of vectors
x
x[x >= 2]
x[x < 3]

# Easy way of work with parts of the data
data <- data.frame(x = c(0, 1, 3, 3, 0), y = 1:5)
data

# Data such that x is zero
data0 <- data[data$x == 0, ]
data0

# Data such that x is larger than 2
data2 <- data[data$x > 2, ]
data2

# Other possibility
subset(x = data, subset = x == 0, select = y)

# Problem -- what happened?
data[x > 2, ]

# AND operator "&"
TRUE & TRUE
TRUE & FALSE
FALSE & FALSE

# OR operator "|"
TRUE | TRUE
TRUE | FALSE
FALSE | FALSE

# Both operators are useful for checking for ranges of data
y
index1 <- (y <= 3) & (y > 0)
y[index1]
index2 <- (y < 2) | (y > 4)
y[index2]

## Plotting functions

# "plot" is the main function for plotting in R
# It has a different behaviour depending on the kind of object
# that it receives

# How to plot some data
plot(iris$Sepal.Length, iris$Sepal.Width,
     main = "Sepal.Length vs. Sepal.Width")

# Change the axis limits
plot(iris$Sepal.Length, iris$Sepal.Width, xlim = c(0, 10),
     ylim = c(0, 10))

# How to plot a curve (a parabola)
x <- seq(-1, 1, l = 50)
y <- x^2
plot(x, y, main = "A red thick parabola", type = "l",
     col = "red", lwd = 3)

# Plotting a more complicated curve between -pi and pi
x <- seq(-pi, pi, l = 50)
y <- (2 + sin(10 * x)) * x^2
plot(x, y, type = "l") # Kind of rough...

# Remember that we are joining points for creating a curve!
# More detailed plot
x <- seq(-pi, pi, l = 500)
y <- (2 + sin(10 * x)) * x^2
plot(x, y, type = "l")

# For more options in the plot customization see
?plot
?par

# "plot" is a first level plotting function. That means that
# whenever is called, it creates a new plot. If we want to
# add information to an existing plot, we have to use a second
# level plotting function such as "points", "lines" or "abline"

plot(x, y) # Create a plot
lines(x, x^2, col = "red") # Add lines
points(x, y + 10, col = "blue") # Add points
abline(a = 5, b = 1, col = "orange", lwd = 2) # Add a straight
# line y = a + b * x
abline(h = 5, col = "green", lwd = 2)
abline(v = 0, col = "yellow", lwd = 2)

# curve
rm(x)
curve(dnorm(x), from = -2, to = 2, n = 200)

## Distributions

# R allows to sample [r], compute density/probability mass
# functions [d], compute distribution function [p], and compute
# quantiles [q] for several continuous and discrete distributions.
# The format employed is [rdpq]name, where name stands for:
# - norm -> Normal
# - unif -> Uniform
# - exp -> Exponential
# - t -> Student's t
# - f -> Snedecor's F
# - chisq -> Chi squared
# - pois -> Poisson
# - binom -> Binomial
# More distributions:
?Distributions

# Sampling from a Normal -- 5 random points from a N(0, 1)
rnorm(n = 5, mean = 0, sd = 1)

# If you want to have always the same result, set the seed of
# the random number generator
set.seed(45678)
rnorm(n = 5, mean = 0, sd = 1)

# Plotting the density of a N(0, 1) -- the Gaussian bell
x <- seq(-4, 4, l = 100)
y <- dnorm(x = x, mean = 0, sd = 1)
plot(x, y, type = "l")

# Plotting the distribution function of a N(0, 1)
x <- seq(-4, 4, l = 100)
y <- pnorm(q = x, mean = 0, sd = 1)
plot(x, y, type = "l")

# Computing the 95% quantile for a N(0, 1)
qnorm(p = 0.95, mean = 0, sd = 1)

# All distributions have the same syntax: rname(n,...),
# dname(x,...), dname(p,...) and qname(p,...), but the
# parameters in ... change. Look them in ?Distributions
# For example, here is que same for the uniform distribution

# Sampling from a U(0, 1)
set.seed(45678)
runif(n = 10, min = 0, max = 1)

# Plotting the density of a U(0, 1)
x <- seq(-2, 2, l = 100)
y <- dunif(x = x, min = 0, max = 1)
plot(x, y, type = "l")

# Computing the 95% quantile for a U(0, 1)
qunif(p = 0.95, min = 0, max = 1)

# Sampling from a Bi(10, 0.5)
set.seed(45678)
samp <- rbinom(n = 200, size = 10, prob = 0.5)
table(samp) / 200

# Plotting the probability mass of a Bi(10, 0.5)
x <- 0:10
y <- dbinom(x = x, size = 10, prob = 0.5)
plot(x, y, type = "h") # Vertical bars

# Plotting the distribution function of a Bi(10, 0.5)
x <- 0:10
y <- pbinom(q = x, size = 10, prob = 0.5)
plot(x, y, type = "h")

## Functions

# A function is a way of encapsulating a block of code so it can
# be reused easily. They are useful for simplifying repetitive
# tasks and organize analyses

# This is a silly function that takes x and y and returns its sum
# Note the use of "return" to indicate what should be returned
add <- function(x, y) {
  z <- x + y
  return(z)
}

# Calling add -- you need to run the definition of the function
# first!
add(x = 1, y = 2)
add(1, 1) # Arguments names can be omitted

# A more complex function: computes a linear model and its
# posterior summary. Saves us a few keystrokes when computing a
# lm and a summary
lmSummary <- function(formula, data) {
  model <- lm(formula = formula, data = data)
  summary(model)
}
# If no return(), the function returns the value of the last
# expression

# Usage
lmSummary(Sepal.Length ~ Petal.Width, iris)

# Recall: there is no variable called model in the workspace.
# The function works on its own workspace!
# model # Error

# Add a line to a plot
addLine <- function(x, beta0, beta1) {
  lines(x, beta0 + beta1 * x, lwd = 2, col = 2)
}

# Usage
x <- 0:10
y <- pbinom(q = x, size = 10, prob = 0.5)
plot(x, y)
addLine(x, beta0 = 0.1, beta1 = 0)

# The function "sapply" allows to sequentially apply a function
sapply(1:10, sqrt)
sqrt(1:10) # The same

# The advantage of "sapply" is that you can use with any function
myFun <- function(x) c(x, x^2)
sapply(1:10, myFun) # Returns a 2 x 10 matrix

# "sapply" is useful for plotting non-vectorized functions
sumSeries <- function(n) sum(1:n)
plot(1:10, sapply(1:10, sumSeries), type = "l")

# "apply" applies iteratively a function to rows (1) or columns
# (2) of a matrix or data frame
A <- matrix(1:10, nrow = 5, ncol = 2)
A
apply(A, 1, sum) # Applies the function by rows
apply(A, 2, sum) # By columns
rowSums(A)
colSums(A)

# With other functions
apply(A, 1, sqrt)
apply(A, 2, function(x) x^2)

## Control structures

# The "for" statement allows to create loops that run along a
# given vector
# Print 3 times a message (i varies in 1:3)

for (i in 1:3) {
  print(i)
}

# Another example

a <- 0
for (i in 1:3) {
  a <- i + a
}
a

# Nested loops are possible
A <- matrix(0, nrow = 2, ncol = 3)
for (i in 1:2) {
  for (j in 1:3) {
    A[i, j] <- i + j
  }
}

# The "if" statement allows to create conditional structures of
# the forms:
# if (condition) {
#  # Something
# } else {
#  # Something else
# }
# These structures are thought to be inside functions

# Is the number positive?
isPositive <- function(x) {
  if (x > 0) {
    print("Positive")
  } else {
    print("Not positive")
  }
}
isPositive(1)
isPositive(-1)

# A loop can be interrupted with the "break" statement
# Stop when x is above 100
x <- 1
for (i in 1:1000) {
  x <- (x + 0.01) * x
  print(x)
  if (x > 100) {
    break
  }
}

## More advanced stuff

rm(a)
# Assignment functions
f1 <- function(x) {
  a <- x # Assign in current environment
  return(0)
}
f2 <- function(x) {
  a <<- x # Assign in upper environment
  return(0)
}

# Check
f1(x = 2)
a
f2(x = 4)
a

# If you want to see the evaluation path
debug(f1)
f1(x = 2)
f1(x = 3)
undebug(f1)
f1(x = 2)