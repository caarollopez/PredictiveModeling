rm(list=ls())

# LAB 9

# EXERCISE 1
data(iris)

# a:

# Standardization to avoid scale effects
pca_iris <- princomp(iris[, -5], cor = TRUE)

# Accumulated variance with 2 PCs: 0.9581321
summary(pca_iris)
cumsum(pca_iris$sdev^2) / sum(pca_iris$sdev^2)

# b:

biplot(pca_iris)
pca_iris$loadings

# Interpretation of components:
# - PC1 captures Petal.* information and, secondarily, Sepal.Length
# - PC2 captures Sepal.Width, primarily, and Sepal.Length, secondarily
# - PC1 separates the two clusters of points, PC2 does not

# Other insights:
# - Sepal.Width and Sepal.Length are almost uncorrelated
# - Sepal.Width and Petal.* are lowly correlated

# c:

# PC1 separates much better the clusters
car::scatterplotMatrix(x = pca_iris$scores, groups = iris$Species)

# Compare it with the original data
car::scatterplotMatrix(x = iris[, -5], groups = iris$Species)

# d:

# Better to load the library as we will call several of its functions
library(pls)

# PCR with l = 2
mod_pcr <- pcr(Petal.Width ~ . - Species, data = iris, scale = TRUE,
               ncomp = 2)
summary(mod_pcr)

# Coefficients for the two components
mod_pcr$coefficients[, , 2]

# Select l by LOOCV
mod_pcr_cv <- pcr(Petal.Width ~ . - Species, data = iris, scale = TRUE,
                  validation = "LOO")
summary(mod_pcr_cv)
validationplot(mod_pcr_cv, val.type = "MSEP")
# One could take l = 1: simpler without "too much" predictive loss, even though
# the minimum is attained for l = 3 components.

# e:

# PLS with l = 2
mod_plsr <- plsr(Petal.Width ~ . - Species, data = iris, scale = TRUE,
                 ncomp = 2)
summary(mod_plsr)

# Coefficients for the two components
mod_plsr$coefficients[, , 2]

# Select l by LOOCV
mod_plsr_cv <- plsr(Petal.Width ~ . - Species, data = iris, scale = TRUE,
                    validation = "LOO")
summary(mod_plsr_cv)
validationplot(mod_plsr_cv, val.type = "MSEP")
# Again, one could take l = 1: simpler without "too much" predictive loss,
# even though the minimum is attained for l = 3 components.

# f:

# PC and PLS directions are very similar -- both PC1 and PLS1 are separating
# the groups
car::scatterplotMatrix(x = rbind(mod_pcr$scores), groups = iris$Species)
car::scatterplotMatrix(x = rbind(mod_plsr$scores), groups = iris$Species)

# g:

# PCA on the predictors used in the model in part d
pca_iris_pl <- princomp(subset(iris, select = -c(Petal.Length, Species)),
                        cor = TRUE)

# The plots differ on:
# (1) the signs of the loadings and scores (which are arbitrarily chosen in
# both cases and depend on the data);
# (2) the scaling of both quantities (the functions use different arrangements).
# However, neither of these differences change qualitatively the resulting
# plot nor the relations between its elements: they are just a reflection
# and rescaling of each other.
biplot(pca_iris_pl)
biplot(mod_pcr, var.axes = TRUE, main = "PCA scores")

# h:

# The biplot can also be used to extract the PLS scores!
# (This is the pls::biplot() function, not stats::biplot())
biplot(mod_plsr, var.axes = TRUE, main = "PLS scores")


# -----------------------------------------------------------------------------
# EXERCISE 2
data("Hitters", package = "ISLR")


# Exclude NAs (one possibility)
Hitters <- na.omit(Hitters)
anyNA(Hitters)

# b:

# There are three groups of highly-correlated variables. From the documentation
# ?ISLR::Hitters we can see that the variables associated to the first block
# (AtBat--Walks) are related with metrics pf the player exclusively from 1986.
# The second block (Years--CWalks) represents metrics from the whole player's
# career (before 1987), so it is expected to have milder correlations between
# the two blocks of variables. The third block (PutOuts--Errors) are again
# metrics from 1986, but they are less correlated with the first block because
# they are related to "defensive" (and not "offensive") actions in baseball.
corrplot::corrplot(cor(subset(Hitters,
                              select = -c(League, Division, NewLeague))))

# c:

set.seed(42)
n <- nrow(Hitters)
ind_train <- sample(x = n, size = round(0.70 * n))
train_Hitters <- Hitters[ind_train, ]
test_Hitters <- Hitters[-ind_train, ]

# d:

# PCR
library(pls)
mod_pcr <- pcr(Salary ~ ., data = train_Hitters, scale = TRUE,
               validation = "LOO")
summary(mod_pcr)

# Cross-validation error
validationplot(mod_pcr)


rmsep_pcr <- RMSEP(mod_pcr)
(l_pcr <- which.min(rmsep_pcr$val[1, 1, ]) - 1) # l = 16 minimizes the error objectively

# e:

# Prediction error
mean((predict(mod_pcr, newdata = test_Hitters,
              ncomp = l_pcr) - test_Hitters$Salary)^2)

# f:

# PLSR
mod_plsr <- plsr(Salary ~ ., data = train_Hitters, scale = TRUE,
                 validation = "LOO")
summary(mod_plsr)

# Cross-validation error
validationplot(mod_plsr)
rmsep_plsr <- RMSEP(mod_plsr)
(l_plsr <- which.min(rmsep_plsr$val[1, 1, ]) - 1) # now we obtain 10 components

# Prediction error
mean((predict(mod_plsr, newdata = test_Hitters,
              ncomp = l_plsr) - test_Hitters$Salary)^2)

# the prediction error is larger now 

# g:

M <- 500
l_pcr <- l_plsr <- rep(NA, M)
error_pcr <- error_plsr <- rep(NA, M)
pb <- txtProgressBar(style = 3)
for (i in 1:M) {
  
  # Random splitting
  ind_train <- sample(x = n, size = round(0.70 * n))
  train_Hitters <- Hitters[ind_train, ]
  test_Hitters <- Hitters[-ind_train, ]
  
  # Model fit
  mod_pcr <- pcr(Salary ~ ., data = train_Hitters, scale = TRUE,
                 validation = "LOO")
  mod_plsr <- plsr(Salary ~ ., data = train_Hitters, scale = TRUE,
                   validation = "LOO")
  
  # Best l's
  l_pcr[i] <- which.min(RMSEP(mod_pcr)$val[1, 1, ]) - 1
  l_plsr[i] <- which.min(RMSEP(mod_plsr)$val[1, 1, ]) - 1
  
  # Prediction errors
  error_pcr[i] <- mean((predict(mod_pcr, newdata = test_Hitters,
                                ncomp = l_pcr[i]) - test_Hitters$Salary)^2)
  error_plsr[i] <- mean((predict(mod_plsr, newdata = test_Hitters,
                                 ncomp = l_plsr[i]) - test_Hitters$Salary)^2)
  
  # Update progress
  setTxtProgressBar(pb, i / M)
  
}

# Mean errors and number of components
mean(error_pcr)
mean(error_plsr) # Smaller, but not so much
mean(l_pcr)
mean(l_plsr) # Significantly smaller
# In conclusion: smaller predictive error and a simpler model with PLSR
# than with PCR


# -----------------------------------------------------------------------------
# EXERCISE 3
#install.packages("pls")
library(pls)
data("gasoline")

# a:

# Huge linear dependency between wavelengths
car::scatterplotMatrix(gasoline$NIR[, 1:10])

# The model "can" be fitted... partially until the fit is perfect (which is
# useless because we are overfitting)
mod <- lm(octane ~ ., data = gasoline)
summary(mod)

# The multicollinearity is perfect and therefore the VIFs are not even defined
car::vif(mod)
# Error in vif.default(mod) : there are aliased coefficients in the model
# "aliased coefficients" means coefficients that are perfectly dependent on
# the others

# b:

# Wavelengths
wv <- seq(900, 1700, by = 2)

# Plot of all the NIR curves
matplot(wv, t(gasoline$NIR), type = "l", lty = 1, xlab = "wavelength",
        ylab = "NIR")
# Each wavelength is a new predictor, the NIR being its value. The different
# curves are different observations. We now see why there is such a high
# dependency between predictors: because they are related by a roughly
# continuous wavelength.

# We could magnify the differences with a logarithmic transformation to have
# a better visualization
matplot(wv, log10(t(gasoline$NIR) + 0.1), type = "l", lty = 1,
        xlab = "wavelength", ylab = "log10(NIR + 0.1)")

# c:

# Color the curves by their octane
col <- viridis::viridis(nrow(gasoline))[rank(gasoline$octane)]
matplot(wv, t(gasoline$NIR), type = "l", lty = 1, xlab = "wavelength",
        ylab = "NIR", col = col)

# Better visualization
matplot(wv, log10(t(gasoline$NIR) + 0.1), type = "l", lty = 1,
        xlab = "wavelength", ylab = "log10(NIR + 0.1)", col = col)
# There seems to be regions where the color gradient is more aligned with
# the height of the curves and the presence of "bands" of NIR associated
# to certain octanes. The third peak from the left seems to be nicely
# capturing the color gradient.

# d:

set.seed(42)
ind_test <- sample(60, size = 18)
train_gas <- gasoline[-ind_test, ]
test_gas <- gasoline[ind_test, ]

# e:

mod_pcr <- pcr(octane ~ ., data = train_gas, scale = TRUE)
summary(mod_pcr)
# l = 4 seems to be the optimal. Beyond l = 4 there seem to be only very
# marginal improvements.

# f:

# Leave-one-out cross-validation
mod_pcr_cv <- pcr(octane ~ ., data = train_gas, scale = TRUE,
                  validation = "LOO")

# The minimum does not seem to be located at l = 4
validationplot(mod_pcr_cv, type = "o")
summary(mod_pcr_cv)

# Extract the RMSEP
(cv_pcr <- RMSEP(mod_pcr_cv))
cv_pcr$val[1, 1, ] # First row
which.min(cv_pcr$val[1, 1, ]) - 1
# l = 8 is the minimum

# l = 8 is the optimal number of components according to LOOCV. However, a
# smaller number, like l = 4, seems to be more reasonable, as it gives a
# similar error for half the complexity.

# g:

biplot(mod_pcr, var.axes = TRUE, cex = 0.5, comps = c(1, 2))
biplot(mod_pcr, var.axes = TRUE, cex = 0.5, comps = c(1, 3))
biplot(mod_pcr, var.axes = TRUE, cex = 0.5, comps = c(2, 3))

# h:

# March along a given component
march_comp <- function(score, mod, comp = 1) {
  if (is.null(mod$scale)) mod$scale <- 1
  mod$scale * (mod$Xmeans + score * mod$loadings[, comp])
}

# Dynamic exploration
library(manipulate)
manipulate({
  matplot(wv, t(gasoline$NIR), type = "l", lty = 1, xlab = "wavelength",
          ylab = "NIR", col = gray(0.5, alpha = 0.5), main = paste0("PC", comp))
  lines(wv, march_comp(score = score, mod = mod_pcr, comp = comp),
        col = 2, lwd = 2)
}, score = slider(min = -200, max = 200, step = 10, initial = 0),
comp = slider(min = 1, max = 5, step = 1, initial = 1))

# Static exploration
col <- viridis::viridis(nrow(gasoline))[rank(gasoline$octane)]
matplot(wv, t(gasoline$NIR), type = "l", lty = 1, xlab = "wavelength",
        ylab = "NIR", col = col, main = "Original data")
scores <- seq(-200, 200, by = 25)
old_par <- par(mfrow = c(2, 2))
cols <- viridis::inferno(length(scores))
for (k in 1:4) {
  matplot(wv, t(gasoline$NIR), type = "l", lty = 1, xlab = "wavelength",
          ylab = "NIR", col = gray(0.5, alpha = 0.5), main = paste0("PC", k))
  for (i in seq_along(scores)) {
    lines(wv, march_comp(score = scores[i], mod = mod_pcr, comp = k),
          col = cols[i])
  }
}
par(old_par)

# Interpretations:
# * PC1 is mainly capturing a variation in height of the NIR curves on all the
#   wavelengths simultaneously. The exception is the secondary peak at ~1150
#   and the rightmost slope past 1600, where the small variation in
#   height is *opposed* to the rest of height variations (observe the inversion
#   of the color gradient). Therefore, PC1 is indicative of the overall
#   vertical position of a curve with respect to the mean of the curves.
# * PC2 is mainly capturing the differential height between the two peaks at
#   1200, the steepness of the right slope of that second peak, the height of
#   the pear at 1400 with respect to the two basins at the left and right
#   (roughly constant), and the variation in shape of the landing at ~1650.
#   Therefore, PC2 is indicative of changes in shape of a curve that change
#   the form of their peaks.
# * PC3 mainly focus on capturing the changes between the wavelengths 1200-1400,
#   excluding the valley. It also captures the changes in height of the basin
#   at ~1650 and changes in the shape of the landing on 1650 (differences on
#   the convex/concave shape).
# * Among other effects, PC4 captures an interesting effect on the peak in 1400:
#   an horizontal shift on the curves.

# Notice how a specific feature (e.g., "vertical variation about 1000") can
# be captured by different principal components (PC1, PC2, PC3, PC4) with
# slightly different variations. This does not contradict the fact that the
# PCs are orthogonal: PCs are a linear combination of many variables, so
# several PCs can be related with a single variable and have a similar
# loading/coefficient for it while still retaining uncorrelation between PCs.

# PC1 is not very useful for explaining the octane. The octane therefore does
# not seem to depend too much on the NIR level. It mostly depends on PC2 and
# PC3, which are related with shape variability of the curves.

# i:

octane_hat_pcr <- predict(mod_pcr, newdata = test_gas, ncomp = 4)
mean((octane_hat_pcr - test_gas$octane)^2)


# ------------------------------------------------------------------------------
# EXERCISE 4 (voluntary)
plsr_mod = plsr(octane~., data = train_gas, scale = TRUE)
summary(plsr_mod)

# Leave-one-out cross-validation
mod_plsr_cv <- pcr(octane ~ ., data = train_gas, scale = TRUE,
                  validation = "LOO")

validationplot(mod_plsr_cv, type = "o")
summary(mod_pcr_cv)

# Extract the RMSEP
(cv_plsr <- RMSEP(mod_plsr_cv))
cv_plsr$val[1, 1, ] # First row
which.min(cv_plsr$val[1, 1, ]) - 1




# ------------------------------------------------------------------------------
# EXERCISE 5
data("yarn")

# a: 
# Wavelengths
#wv <- seq(900, 1700, by = 2) 

# Color the curves by their density
col <- viridis::viridis(nrow(yarn))[rank(yarn$density)]
#matplot(wv, t(yarn$NIR), type = "l", lty = 1, xlab = "wavelength",
        ylab = "NIR", col = col)

# b: 
# split in training and testing

train_idx <- which(yarn$train == "train")
test_idx <- which(yarn$train == "test")

train_data <- yarn[train_idx,]
test_data <- yarn[test_idx,]

# c:
# PCR and PLSR density 
pcr_mod <- pcr(density ~ NIR, data = train_data, scale = TRUE, validation = "LOO")
summary(pcr_mod)

plsr_mod <- plsr(density ~ NIR, data = train_data, scale = TRUE, validation = "CV")
summary(plsr_mod)

# d: 
# lhats that give the lowest 10-cross fold validation error

# e: 
# three first PCSs and PLs using march_comp()

# f: 
# predictice mean squared error for PCR and PLSR on the test dataset
# for the selected l, l hats




# -----------------------------------------------------------------------------
# EXERCISE 6

# Implementation part

bic_ci_boot <- function(formula, data, B = 500, alpha = 0.05) {
  
  # Step 1
  bic_0 <- BIC(lm(formula = formula, data = data))
  
  # Step 2
  n <- nrow(data)
  bic_star <- numeric(B)
  for (b in 1:B) {
    
    # Step 2.1
    data_star <- data[sample(n, replace = TRUE), ]
    
    # Step 2.2
    bic_star[b] <- BIC(lm(formula = formula, data = data_star))
    
  }
  
  # Step 3
  ci <- quantile(bic_star, prob = c(alpha / 2, 1 - alpha / 2), names = FALSE)
  return(data.frame("BIC" = bic_0, "lower" = ci[1], "upper" = ci[2]))
  
}

# A simple check
bic_ci_boot(Sepal.Length ~ ., data = iris)

# Validation part


# -----------------------------------------------------------------------------
# EXERCISE 7


