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

# POISSON BAYES ####
p1b <- nimbleCode({
  
  # priors
  B0 ~ dnorm(mean = 0, sd = 100)
  B1[1] <- 0 # 'corner constraints' so we can est. params (identifiability)
  B1[2] ~ dnorm(mean = 0, sd = 100)
  

  # likelihood
  for(i in 1:nObs){
    C[i] ~ dpois(lambda[i])
    log(lambda[i]) <- B0 + B1[forest[i]]
  
  
  # Fit assessments
  # Pearson goodness of fit
  # Residual check for deviation between observed and predicted
    Pres[i] <- (C[i] - lambda[i] / sqrt(lambda[i]))            # pearson residual
    C.new[i] ~ dpois(lambda[i])                                # Replicate Dataset - generated under perfect conditions
    Pres.new[i] <- (C.new[i] - lambda[i] / sqrt(lambda[i]))    # pearson residual for new data
    sq.res[i] <- pow(Pres[i],2)                                # squared residuals
    sq.res.new[i] <- pow(Pres.new[i], 2)
  }
  
  fit <- sum(sq.res[1:nObs])
  fit.new <- sum(sq.res.new[1:nObs])
  
})

nimData <- list(C = C)
nimConsts <- list(nObs = length(C),
                  forest = ifelse(forest == 'BLH', 1, 2))

nimInits <- list(B0 = rnorm(1,0,100),
                 B1 = c(NA, rnorm(1, 0, 100)))

keepers <- c('B0', 'B1', 'lambda', 'Pres', 'fit', 'fit.new')

p1b <- nimbleMCMC(code = p1b,
                       data = nimData,
                       constants = nimConsts,
                       monitors = keepers,
                       niter = 6000,
                       nburnin = 1000,
                       thin = 1,
                       nchains = 3,
                       summary = T)

p1b$summary$all.chains
p1

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

# visualize relationship
ggplot(exp.m, aes(voc,observed))+theme_bw()+
  geom_point(position = position_jitter(w = 2, h = 0.05), size=3) +
  geom_smooth(colour="red") + xlab("Visual Obstruction") +
  ylab("Detection = 1")


# fit freq model
m1 <- glm(observed ~ voc, data = moose, family = binomial())
summary(m1)

(ci.m1 <- confint(m1))

# report odds instead of log odds
# "odds of detection increases/decreases by..."
exp(coef(m1))
exp(ci.m1)

# these are not probabilities! They are log-odds and odds
# can transform using inverse-logit to get probabilities
inv.logit <- function(x){
  exp(x) / (1+exp(x))
}
inv.logit(coef(m1)[1])
plogis(coef(m1)[1])

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
  }
  
  
})


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc)

nimInits <- list(B0 = rnorm(1),
                 B1 = rnorm(1))

keepers <- c('B0', 'B1')

m2 <- nimbleMCMC(code = m2code,
                  data = nimData,
                  constants = nimConsts,
                  monitors = keepers,
                  niter = 6000,
                  nburnin = 1000,
                  thin = 1,
                  nchains = 3,
                  summary = T)

m2$summary$all.chains


# get samples usable
samples_mcmc <- coda::as.mcmc.list(lapply(m2$samples, coda::mcmc))

samples <- do.call(rbind, samples_mcmc)

# BEFORE LOOKING AT ESTIMATES, two things should be done:
# 1. Assess model convergence
# 2. Assess whether fitted model is adequate for dataset

par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])

coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0'))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])


# PREDICTIONS FROM MODEL

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

# make a vector of all possible VOC values
all_voc <- seq(0, 100, 1)


nimData <- list(y = moose$observed)
nimConsts <- list(nObs = nrow(moose),
                  voc = moose$voc,
                  pred.voc = all_voc,     # bring in predicted values here
                  nPred = length(all_voc)) 

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
predicted_p <- samples_mcmc[, which(stringr::str_detect(string = colnames(samples_mcmc[[1]]), pattern = 'p.pred\\['))]
str(predicted_p)

# learn
dev.off()
test <- predicted_p[[1]] # one chain
test <- predicted_p[[1]][,1] # one chain, the first value of VOC (VOC = 0)
hist(test)
test <- predicted_p[[1]][1,] # one chain, first iteration for all sites

for(i in 1:500){
  hold <- predicted_p[[1]][i,] # one chain, one iteration for all sites
  
  if(i == 1){
    plot(x = all_voc, y = hold, type = 'l', lty = 1, lwd = 2, col = alpha(randomcoloR::randomColor(count = 1), 0.2),
         xlab = 'Visual Obstruction', ylab = 'Detection Probability')
  } else{
    lines(all_voc, hold, lty = 1, lwd = 2, col = alpha(randomcoloR::randomColor(count = 1), 0.2))
  }
}

# add in mean
mean_p <- apply(predicted_p[[1]],2,mean)
lines(all_voc, mean_p, lty = 1, lwd = 6, col = 'black')
