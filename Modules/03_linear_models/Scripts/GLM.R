# GLMs
library(abd)
library(ggplot2)
library(dplyr)
library(nimble)
library(here)
library(randomcoloR)


## ------------------ ##
## --POISSON DIST--   ##
## ------------------ ##

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

# POISSON BAYES ####

# Notes on bayesian posterior predictive checks
# Posterior predictive checks are used to evaluate
# how well a model fits the observed data.
# The basic idea is to simulate new data (observations)
# from the model's posterior distribution and compare
# the new data to the observed data

# If the simulated data from the posterior
# are similar to the observed data, then this
# suggests model is reasonable for data.
# however, simulated data that is quite different
# from observed may suggest poor fitting model given
# observed data.

# By generating simulated datasets 
# from the model's posterior predictive 
# distribution, we can see how well the 
# model captures the patterns and variability
# in the observed data, 
# and we can use this information to 
# improve our understanding of the 
# underlying process we're trying to model.

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
  
})

nimData <- list(C = C)
nimConsts <- list(nObs = length(C),
                  forest = ifelse(forest == 'BLH', 1, 2))

nimInits <- list(B0 = rnorm(1,0,10),
                 B1 = c(NA, rnorm(1, 0, 10)))

keepers <- c('B0', 'B1', 'lambda', 'Pres', 'fit', 'fit.new')

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
p1b$summary$all.chains
coef(p1)

# get samples usable
samples_mcmc <- coda::as.mcmc.list(lapply(p1b$samples, coda::mcmc))
samples <- do.call(rbind, samples_mcmc)


# BEFORE LOOKING AT ESTIMATES, two things should be done:
# 1. Assess model convergence
# 2. Assess whether fitted model is adequate for dataset

par(mfrow=c(1,3))
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1\\[2'))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'lambda\\[1\\]'))])

# Check Rhat
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1\\[2'))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'lambda\\[1\\]'))])


# 2a. Plot residuals
dev.off()
plot(apply(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'Pres\\['))], 2, mean), las = 1)
abline(h = 0)

# 2a. Posterior Predictive check
plot(samples[, 'fit'],
     samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.new'))], xlab = 'Discrepancy measure for actual data',
     ylab = 'Discrepancy for perfect data')
abline(0,1, lwd = 2, col = "black")

# 2a. Bayes p-value
# 2a. Proportion of simu;ated datasets that are equally or more extreme than observed.
# 2a. Target is 0.5 = model is performing as expected
# 2a. Deviation from 0.5 is okay, no strict threshold for rejection. 0.3 - 0.7 is fine, use judgement as you approach bounds
mean(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.new'))] > samples[, 'fit'])

# Histogram Colored (blue and red)
hist(exp(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))]), col=rgb(1,0,0,0.5),xlim=c(0,10), ylim=c(0,4000),
     main = 'Effect of veg on swamp rabbit density', xlab = 'Expected Swamper Count')

# get value of b0 + b1
this <- samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))] + samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1\\[2\\]'))]
hist(exp(this), col=rgb(0,0,1,0.5), add=T)
box()
legend("topleft", legend = levels(forest), col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pch = 16)


## ------------------ ###
## --BINOMIAL DIST--  ####
## ------------------ ###


# Moose detection

moose <- read.table(here('Modules/03_linear_models/Data/moose.txt'))
str(moose)
head(moose)

# visualize relationship
ggplot(moose, aes(voc,observed))+theme_bw()+
  geom_point(position = position_jitter(w = 2, h = 0.05), size=3) +
  geom_smooth(colour="red") + xlab("Visual Obstruction") +
  ylab("Detection = 1")


# !! Define data/problem/model



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

# can use predict if still feeling uncomfortable
predict(m1, data.frame(voc = 0), type = 'link')
predict(m1, data.frame(voc = 0), type = 'response')
# model predicts 85% chance of detecting moose when it is completely out in the open (no veg, voc = 0)

ggplot(moose, aes(voc,observed)) + theme_bw() + 
  geom_point(position = position_jitter(w = 2, h = 0.05), size=3) +
  xlab("Visual Obstruction") + 
  stat_smooth(method="glm", method.args = list(family = "binomial") )+
  ylab("Probability of Detection") 


# BINOM BAYES ####

m2code <- nimbleCode({
  
  # priors
  B0 ~ dnorm(0, sd = 10)
  B1 ~ dnorm(0, sd = 10)

  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- B0 + B1*voc[i]

  # Fit assessments
  # Computation of fit statistic (for Bayesian p-value)
  Presi[i] <- abs(y[i]-p[i])	 # Absolute residual
  y.new[i]~ dbern(p[i])
  Presi.new[i] <- abs(y.new[i]-p[i])
}


fit <- sum(Presi[1:nObs])# Discrepancy for actual data set
fit.new <- sum(Presi.new[1:nObs]) 		# Discrepancy for replicate data set


})


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc)

nimInits <- list(B0 = rnorm(1),
                 B1 = rnorm(1))

keepers <- c('B0', 'B1', 'fit', 'fit.new')

m2 <- nimbleMCMC(code = m2code,
                  data = nimData,
                  constants = nimConsts,
                  monitors = keepers,
                  niter = 6000,
                  nburnin = 1000,
                  thin = 1,
                  nchains = 3,
                  summary = T)

# get samples usable
samples_mcmc <- coda::as.mcmc.list(lapply(m2$samples, coda::mcmc))

samples <- do.call(rbind, samples_mcmc)

# Assess model convergence

par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])

coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])

# Assess model fit
# 2a. Posterior Predictive check
plot(samples[, 'fit'],
     samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.new'))], xlab = 'Discrepancy measure for actual data',
     ylab = 'Discrepancy for perfect data')
abline(0,1, lwd = 2, col = "black")

# Bayes P-value
# proportion of simulated datasets that are as or more
# extreme than the observed data
mean(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.new'))] > samples[, 'fit'])

# Inspect model output
m2$summary$all.chains


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Great, let's use the output to learn about the problem
# For example, we can predict response given a new dataset


# make a vector of all possible VOC values for predictions
pred.voc <- seq(0, 100, 1)

# predict from model
m2code <- nimbleCode({
  
  # priors
  B0 ~ dnorm(0, sd = 10)
  B1 ~ dnorm(0, sd = 10)
  
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- B0 + B1*voc[i]
  }
  
  #predictions
  for(k in 1:nPred){
    logit(p.pred[k]) <- B0 + B1 * pred.voc[k]
  }
})


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc,
                  pred.voc = pred.voc,     # bring in predicted values here
                  nPred = length(pred.voc)) 

nimInits <- list(B0 = rnorm(1),
                 B1 = rnorm(1))

keepers <- c('B0', 'B1', 'p.pred') # monitor p.pred

m2 <- nimbleMCMC(code = m2code,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 niter = 6000,
                 nburnin = 1000,
                 thin = 1,
                 nchains = 3,
                 summary = T)


samples_mcmc <- coda::as.mcmc.list(lapply(m2$samples, coda::mcmc))
# pull out predicted ps
predicted_p <- samples_mcmc[, which(stringr::str_detect(string = colnames(samples_mcmc[[1]]), pattern = 'p.pred\\['))]
# see what it looks like
# three chains with 5,000 iterations
# in R, list of 3 matrices with structure of iterations on rows and predicted ps on columns 5,000 x 101
str(predicted_p)

# Visualize predictions with uncertainty
dev.off()
pp <- predicted_p[[1]] # one chain
pp <- predicted_p[[1]][,1] # one chain, the first value of VOC (VOC = 0)
hist(pp)
pp <- predicted_p[[1]][1,] # one chain, first iteration for all sites

for(i in 1:500){
  hold <- predicted_p[[1]][i,] # one chain, one iteration for all sites
  
  if(i == 1){
    plot(x = pred.voc, y = hold, type = 'l', lty = 1, lwd = 2, col = alpha(randomcoloR::randomColor(count = 1), 0.2),
         xlab = 'Visual Obstruction', ylab = 'Detection Probability', ylim = c(0,1))
  } else{
    lines(pred.voc, hold, lty = 1, lwd = 2, col = alpha(randomcoloR::randomColor(count = 1), 0.2))
  }
}

# add in mean and 95% CI
mean_p <- apply(predicted_p[[1]],2,mean)
low <- apply(predicted_p[[1]],2,FUN = function(x) quantile(x, probs = c(0.025)))
high <- apply(predicted_p[[1]],2,FUN = function(x) quantile(x, probs = c(0.975)))

lines(pred.voc, mean_p, lty = 1, lwd = 6, col = 'black')
lines(pred.voc, low, lty = 2, lwd = 4, col = 'black')
lines(pred.voc, high, lty = 2, lwd = 4, col = 'black')
