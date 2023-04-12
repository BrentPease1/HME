# Assignment 03 - GLM solution
library(abd)
library(ggplot2)
library(dplyr)
library(nimble)
library(here)
library(randomcoloR)

#Q1 ####
#Modify the swamp rabbit Bayesian model to directly calculate 
# A) the mean abundance/density in bottomland hardwood forests,
# B) mean abundance/density in upland mixed forests, 
# and C) the mean difference in abundance between the two forest types. 
# Report the point estimate and the 95% credible interval with an interpretation of the parameter estimates.


# COUNTS OF RABBITS #

sites <- 20
forest <- gl(n = 2, k = sites, labels = c('BLH', 'UMF'))
nObs <- 2*sites

# mean density of swampers in BLH and UMF is 5 and 2, respectively.
mean_blh <- 5
mean_umf <- 2

# generate intensity
# b0 = log(mean_blh) = 1.609
# b1 = log(mean_umf) = b0 + b1, b1 = log(mean_blh) - log(umf)
lam <- exp(1.609 - 0.9162907*(as.numeric(forest) - 1))


# add poisson noise and generate counts
C <- rpois(n = nObs, lambda = lam)

# observed means
aggregate(C, by = list(forest), FUN = mean) # The observed means
boxplot(C ~ forest, col = "grey", xlab = "Forest Type", ylab = "Swampper count", las = 1)

# fit the model
p1 <- glm(C ~ forest, family = poisson)
summary(p1)
anova(p1, test = "Chisq") # Likelihood ratio test (LRT)

# return true mean swamper density 
exp(abs(coef(p1)))
exp(abs(confint(p1)))



p1b <- nimbleCode({
  
  # priors
  B0 ~ dnorm(mean = 0, sd = 10)
  B1[1] <- 0 # 'corner constraints' so we can est. params (identifiability)
  B1[2] ~ dnorm(mean = 0, sd = 10)
  
  
  # likelihood
  for(i in 1:nObs){
    C[i] ~ dpois(lambda[i])
    log(lambda[i]) <- B0 + B1[forest[i]]
    
    
    # Fit assessments
    # Pearson goodness of fit
    # Residual check for deviation between observed and predicted
    Pres[i] <- (C[i] - lambda[i] / sqrt(lambda[i]))            # pearson residual
    sq.res[i] <- pow(Pres[i],2)                                # squared residuals
    C.new[i] ~ dpois(lambda[i])                                # Replicate Dataset - generated under perfect conditions
    Pres.new[i] <- (C.new[i] - lambda[i] / sqrt(lambda[i]))    # pearson residual for new data
    sq.res.new[i] <- pow(Pres.new[i], 2)
  }
  
  fit <- sum(sq.res[1:nObs])
  fit.new <- sum(sq.res.new[1:nObs])
  
  
  # Derived values
  mean_blh <- exp(B0)
  mean_umf <- exp(B0 + B1[2])
  mean_diff <- mean_blh - mean_umf
  
})

nimData <- list(C = C)
nimConsts <- list(nObs = length(C),
                  forest = ifelse(forest == 'BLH', 1, 2))

nimInits <- list(B0 = rnorm(1,0,10),
                 B1 = c(NA, rnorm(1, 0, 10)))

keepers <- c('mean_blh', 'mean_umf', 'mean_diff')

p1b <- nimbleMCMC(code = p1b,
                  data = nimData,
                  constants = nimConsts,
                  monitors = keepers,
                  inits = nimInits,
                  niter = 6000,
                  nburnin = 1000,
                  thin = 1,
                  nchains = 3,
                  summary = T)

# look at estimates
head(p1b$summary$all.chains)


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Q2 ####

# Modify the moose model to "derive" the regression coefficients as probabilities 
# instead of log odds. 
# Again, for each regression coefficient, report the mean and 95% credible interval

# Moose detection

moose <- read.table(here('Modules/03_linear_models/Data/moose.txt'))


# fit freq model
m1 <- glm(observed ~ voc, data = moose, family = binomial())
summary(m1)


(ci.m1 <- confint(m1))

# report odds instead of log odds
# "odds of detection increases/decreases by..."
exp(coef(m1))
exp(ci.m1)

# these are not probabilities! They are log-odds and odds
# can transform to  probabilities using inverse-logit 
# if logit is log(x) / (1 - log(x)), then inverse of logit is:
inv.logit <- function(x){
  exp(x) / (1+exp(x))
}

# intercept converted to probability scale
inv.logit(coef(m1)[1])

# same as above
plogis(coef(m1)[1])


m2code <- nimbleCode({
  
  # priors
  B0 ~ dnorm(0, sd = 10)
  B1 ~ dnorm(0, sd = 10)
  
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- B0 + B1*voc[i]
  
  }

  B0_prob <- exp(B0) / (1+exp(B0))
  B1_prob <- exp(B1) / (1+exp(B1))
  
  
})


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc)

nimInits <- list(B0 = rnorm(1),
                 B1 = rnorm(1))

keepers <- c('B0', 'B1', 'B0_prob', 'B1_prob')

m2 <- nimbleMCMC(code = m2code,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 inits = nimInits,
                 niter = 6000,
                 nburnin = 1000,
                 thin = 1,
                 nchains = 3,
                 summary = T)

# output
m2$summary$all.chains

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Q3 ####
#Within the moose model, 
# write code to estimate the probability of detecting a moose at 100% vegetation cover. 
# Report the probability estimate with a 95% credible interval

m2code <- nimbleCode({
  
  # priors
  B0 ~ dnorm(0, sd = 10)
  B1 ~ dnorm(0, sd = 10)
  
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- B0 + B1*voc[i]
    
  }
  
  B0_prob <- exp(B0) / (1+exp(B0))
  B1_prob <- exp(B1) / (1+exp(B1))
  
  full_cover <- exp(B0 + B1*100) / (1+exp(B0 + B1*100))
  
  #OR

  logit(p.full) <- B0 + B1 * 100
  
})


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc)

nimInits <- list(B0 = rnorm(1),
                 B1 = rnorm(1))

keepers <- c('B0', 'B1', 'B0_prob', 'B1_prob', 'full_cover', 'p.full')

m2 <- nimbleMCMC(code = m2code,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 inits = nimInits,
                 niter = 6000,
                 nburnin = 1000,
                 thin = 1,
                 nchains = 3,
                 summary = T)

# output
m2$summary$all.chains


#Q4 ####
#For the Poisson example, 
# we generated a dataset and then tested whether we could return truth from a model. 
# In a similar fashion, generate a simple dataset of species occurrence 
# that is a function of at least one ecological covariate
# and write a model in nimble to reflect the data generating process. Does your model return truth?

M <- 200
forest <- seq(-1, 1, l = M)
forest <- rnorm(M, mean = 0, sd = 1)

B0 <- 0.25
B1 <- 1.2

#true occurrence prob
psi <- plogis(B0 + B1 * forest)

# generate true occurrence
y <- rbinom(n = M, size = 1, prob = psi)

m1 <- glm(y ~ forest, family = binomial)
summary(m1)
