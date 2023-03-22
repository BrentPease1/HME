# Linear Regression Review
library(abd)
library(ggplot2)
library(dplyr)
library(nimble)

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# SIMPLE LINEAR REGRESSION                                                  ####
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

data("LionNoses")
#These data come from a paper by Whitman, Starfield, Quadling, & Packer (2004), 
# in which the authors address the impact of trophy hunting on lion population dynamics. 
# The authors note that removing male lions can lead to increases in infanticide, 
# but the authors’ simulation models suggest that removal of only older males 
# (e.g., greater than 6 years of age) could minimize these impacts.3

# How could a researcher (or hunter/guide) tell the age of a lion from afar, though?
# It turns out that it is possible to get a rough idea of how old a male lion is 
# from the amount of black pigmentation on its nose

head(LionNoses)

# plot the relationship
ggplot(LionNoses, aes(proportion.black, age)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)+ xlab("Proportion Black") +
  ylab("Age")

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## FREQ
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

lm.noses <- lm(age ~ proportion.black, data = LionNoses)
lm.noses

# easier to interpret unit change with percentages
LionNoses <- LionNoses %>% 
  mutate(percentage.black = 100*proportion.black)
lm.noses2<-lm(age ~ percentage.black, data=LionNoses)
lm.noses2

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## NIMBLE
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

m1 <- nimbleCode({
  
  # priors
  B0 ~ dnorm(mean = 0, sd = 5)
  B1 ~ dnorm(mean = 0, sd = 5)
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0 + B1*percentage.black[i]
  }
  
  # SSE
  for(i in 1:nObs){
    sq.res[i] <- pow(y[i] - mu[i], 2)
  }
  
  mse <- sum(sq.res[1:nObs]) / (nObs - 2)
  manual.sig <- sqrt(mse)
  
})

nimData <- list(y = LionNoses$age)
nimConsts <- list(nObs = length(LionNoses$age),
                  percentage.black = LionNoses$percentage.black)

nimInits <- list(B0 = rnorm(1,0,5),
                 B1 = rnorm(1,0,5),
                 tau = rgamma(1,1,1))

keepers <- c('B0', 'B1', 'sig', 'manual.sig')

nim.noses <- nimbleMCMC(code = m1,
                        constants = nimConsts,
                        data = nimData,
                        inits = nimInits,
                        monitors = keepers,
                        niter = 7000,
                        nburnin = 2000,
                        thin = 1,
                        nchains = 3,
                        summary = T)

# Don't be bayesic, check your traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.noses$samples, coda::mcmc))

par(mfrow=c(1,4))
coda::traceplot(samples_mcmc)

# Check Rhat
coda::gelman.diag(samples_mcmc)

# Checks out, check parameter estimates
nim.noses$summary$all.chains

# compare to freq analysis
summary(lm.noses2)
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## MULTIPLE LINEAR REGRESSION                                               ####
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

## TWO CONTINUOUS ####
# fake to discuss parameter interp with more than 1 continuous predictor
elev <- runif(100,-1,1)
forest <- runif(100,-1,1)
b0 <- 2
b1 <- -1.1
b2 <- 1.2
m <- b0 + b1*elev + b2*forest
y <- rnorm(100, mean = m, sd = 1)
lm.fake <- lm(y ~ elev + forest)
summary(lm.fake)


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# GOLDEN JACKEL MANDIBLE - one categorical ####
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
males<-c(120, 107, 110, 116, 114, 111, 113, 117, 114, 112)
females<-c(110, 111, 107, 108, 110, 105, 107, 106, 111, 111)
jack <- data.frame(jaws = c(males, females),
                   sex = c(rep('male',10), rep('female', 10)))
head(jack)


## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## FREQ
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
lm.cat <- lm(jaws ~ sex, data = jack)
summary(lm.cat)

mean(females)
mean(males)
mean(females) + lm.cat$coefficients[2]
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## NIM
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# effects parameterization
# the model coefficients represent the effect of each level of the predictor variable, 
# relative to the overall mean of the response variable. 

m2 <- nimbleCode({
  
  # priors
  B0 ~ dnorm(mean = 0, sd = 5)
  B1[1] <- 0 # 'corner constraints' for parameter identifiability
  B1[2] ~ dnorm(mean = 0, sd = 5)


  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0 + B1[sex[i]]
  }
  
})

nimData <- list(y = jack$jaws)
nimConsts <- list(nObs = length(jack$jaws),
                  int_mean = mean(jack[jack$sex == 'female', 1]),
                  sex = ifelse(jack$sex == 'female', 1,2))

nimInits <- list(B0 = rnorm(1,0,5),
                 B1 = c(NA, rnorm(1,0,5)),
                 tau = rgamma(1,1,1))

keepers <- c('B0', 'B1', 'sig')

nim.jaws <- nimbleMCMC(code = m2,
                        constants = nimConsts,
                        data = nimData,
                        inits = nimInits,
                        monitors = keepers,
                        niter = 7000,
                        nburnin = 2000,
                        thin = 1,
                        nchains = 3,
                        summary = T)

# traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.jaws$samples, coda::mcmc))

par(mfrow=c(2,2))
coda::traceplot(samples_mcmc[, 1:4])

# Check Rhat
coda::gelman.diag(samples_mcmc[, c(1,3:4)]) #check column index

# Checks out, check parameter estimates
nim.jaws$summary$all.chains

# check freq
summary(lm.cat)



## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# means parameterization
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# the model coefficients represent the expected value of the response variable 
# for each level of the predictor variable


lm.cat2 <- lm(jaws ~ sex-1, data = jack)
summary(lm.cat2)

m3 <- nimbleCode({
  
  # priors
  for(k in 1:2){
    B0[k] ~ dnorm(mean = 100, sd = 50)
  }
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0[sex[i]]
  }
  
})

nimData <- list(y = jack$jaws)
nimConsts <- list(nObs = length(jack$jaws),
                  sex = ifelse(jack$sex == 'female', 1,2))

nimInits <- list(B0 = rnorm(2,100,50),
                 tau = rgamma(1,1,1))

keepers <- c('B0','sig')

nim.jaws.means <- nimbleMCMC(code = m3,
                        constants = nimConsts,
                        data = nimData,
                        inits = nimInits,
                        monitors = keepers,
                        niter = 6000,
                        nburnin = 1000,
                        thin = 1,
                        nchains = 3,
                        summary = T)

# check your traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.jaws.means$samples, coda::mcmc))

par(mfrow=c(1,3))
coda::traceplot(samples_mcmc)

# Check Rhat
coda::gelman.diag(samples_mcmc) #check column index

# Checks out, check parameter estimates
nim.jaws.means$summary$all.chains
summary(lm.cat2)

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# CHICKADEES - one continuous and one categorical ####
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Templeton, C.N., E. Greene, and K. Davis. 2005. 
# Allometry of alarm calls: Black-capped Chickadees encode information about predator size. Science 308: 1934-1937.


spp <- factor(sample(c("BCCH", "CACH"), 100, replace = T))
#spp <- ifelse(spp == "BCCH", 1, 0)
pred_mass <- seq(0,2, l = 100)

b0 <- 3
b1 <- -1.1
b2 <- 2

m <- b0 + b1*pred_mass + ifelse(spp == 'BCCH', 0, 2)

dees <- rnorm(100, mean = m, sd = 1)

# Create scatter plot
plot(pred_mass, dees, col = as.numeric(spp), pch = 16, xlab = "x", ylab = "y")

# Add legend
legend("topleft", legend = levels(spp), col = 1:length(levels(spp)), pch = 16)


mr <- lm(dees ~ pred_mass + spp)
mr2 <- lm(dees ~ pred_mass + spp - 1)

# CHICKADEES - effects parameterization ####

m4 <- nimbleCode({
  
  # priors
  B0 ~ dnorm(mean = 0, sd = 10) 
  B1[1] <- 0 # 'corner constraints' for parameter identifiability
  B1[2] ~ dnorm(mean = 0, sd = 5)
  B2 ~ dnorm(mean = 0, sd = 10)
  
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0 + B1[spp[i]] + B2*pred_mass[i]
  }
  
})

# CHICKADEES - means parameterization ####

m5 <- nimbleCode({
  
  # priors
  for(k in 1:2){
    B0[k] ~ dnorm(mean = 0, sd = 10)
  }
  B1 ~ dnorm(mean = 0, sd = 10)
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0[spp[i]] + B1*pred_mass[i]
  }
  
})


nimData <- list(y = dees)
nimConsts <- list(nObs = length(dees),
                  spp = ifelse(spp == 'BCCH',1,2),
                  pred_mass = pred_mass)

nimInits <- list(B0 = rnorm(2,0,10),
                 B1 = rnorm(1, 0, 10),
                 tau = rgamma(1,1,1))

keepers <- c('B0','B1', 'sig')

nim.dees <- nimbleMCMC(code = m5,
                             constants = nimConsts,
                             data = nimData,
                             inits = nimInits,
                             monitors = keepers,
                             niter = 6000,
                             nburnin = 1000,
                             thin = 1,
                             nchains = 3,
                             summary = T)

nim.dees$summary$all.chains
summary(mr2)

# NOTE: HOW COULD WE CODE UP CATEGORICAL * CONTINUOUS INTERACTION?
# means param
like_this <- nimbleCode({
    
    # priors
    for(k in 1:2){
      B0[k] ~ dnorm(mean = 0, sd = 10)
      B1[k] ~ dnorm(mean = 0, sd = 10)
    }
    
    tau ~ dgamma(1,1)
    sig <- 1/sqrt(tau)
    
    # likelihood
    for(i in 1:nObs){
      y[i] ~ dnorm(mean = mu[i], sd = sig)
      mu[i] <- B0[spp[i]] + B1[spp[i]]*pred_mass[i]
    }
    
  })

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## TWO CATEGORICAL
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Create the first categorical predictor
veg <- sample(c("forest", "grassland"), 100, replace = T)

# Create the second categorical predictor
crp <- sample(c("enrolled", "not"), 100, replace = T)


nests <- rep(NA, 100)
for (i in 1:100) {
  if (veg[i] == "forest" && crp[i] == "enrolled") {
    nests[i] <- rnorm(1, mean = 20, sd = 2)
  } else {
    nests[i] <- rnorm(1, mean = 10, sd = 3)
  }
}

# EFFECTS PARAMETERIZATION
cat.mr <- lm(nests ~ veg + crp)
summary(cat.mr)
predict(cat.mr, data.frame(veg = 'forest', crp = 'enrolled'))
predict(cat.mr, data.frame(veg = 'forest', crp = 'not'))
predict(cat.mr, data.frame(veg = 'grassland', crp = 'enrolled'))
predict(cat.mr, data.frame(veg = 'grassland', crp = 'not'))




m6 <- nimbleCode({
  
  # priors
  B0 ~ dnorm(0, sd = 10) # intercept
  B1[1] <-  0 # corner constraints
  B1[2] ~ dnorm(0, sd = 10) # effect of vegetation type
  B2[1] <- 0
  B2[2] ~ dnorm(0, sd = 10) # effect of crp
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0 + B1[veg[i]] + B2[crp[i]]
  }
  
})

nimData <- list(y = nests)
nimConsts <- list(nObs = length(nests),
                  veg = ifelse(veg == 'forest',1,2),
                  crp = ifelse(crp == 'enrolled',1,2))

nimInits <- list(B0 = rnorm(1,0,10),
                 B1 = c(NA, rnorm(1, 0, 10)),
                 B2 = c(NA, rnorm(1, 0, 10)),
                 tau = rgamma(1,1,1))

keepers <- c('B0','B1', 'B2', 'sig')

nim.nests <- nimbleMCMC(code = m6,
                       constants = nimConsts,
                       data = nimData,
                       inits = nimInits,
                       monitors = keepers,
                       niter = 6000,
                       nburnin = 1000,
                       thin = 1,
                       nchains = 3,
                       summary = T)


# check your traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.nests$samples, coda::mcmc))

par(mfrow=c(2,3))
coda::traceplot(samples_mcmc)

# Check Rhat
coda::gelman.diag(samples_mcmc[, c(1,3,5:6)]) #check column index


nim.nests$summary$all.chains
summary(cat.mr)


# MEANS PARAMETERIZATION

cat.mr2 <- lm(nests ~ veg + crp-1)
summary(cat.mr2)
predict(cat.mr2, data.frame(veg = 'forest', crp = 'enrolled'))
predict(cat.mr2, data.frame(veg = 'grassland', crp = 'enrolled'))
predict(cat.mr2, data.frame(veg = 'forest', crp = 'not'))
predict(cat.mr2, data.frame(veg = 'grassland', crp = 'not'))



m7 <- nimbleCode({
  
  # priors
  for(k in 1:2){
    B1[k] ~ dnorm(0, sd = 10) # effect of vegetation type
    B2[k] ~ dnorm(0, sd = 10) # effect of crp
  }

  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B1[veg[i]] + B2[crp[i]]
  }
  
  # derived values
  forest_crp <- B1[1] + B2[1] 
  forest_no <- B1[1] + B2[2] 
  grass_crp <- B1[2] + B2[1] 
  grass_no <- B1[2] + B2[2] 
  
})

nimData <- list(y = nests)
nimConsts <- list(nObs = length(nests),
                  veg = ifelse(veg == 'forest',1,2),
                  crp = ifelse(crp == 'enrolled',1,2))

nimInits <- list(B1 = rnorm(2, 0, 10),
                 B2 = rnorm(2, 0, 10),
                 tau = rgamma(1,1,1))

#keepers <- c('B1', 'B2', 'sig', 'forest_crp', 'forest_no', 'grass_crp', 'grass_no')
keepers <- c('forest_crp', 'forest_no', 'grass_crp', 'grass_no','sig')

nim.nests.m <- nimbleMCMC(code = m7,
                        constants = nimConsts,
                        data = nimData,
                        inits = nimInits,
                        monitors = keepers,
                        niter = 6000,
                        nburnin = 1000,
                        thin = 1,
                        nchains = 3,
                        summary = T)

samples_mcmc <- coda::as.mcmc.list(lapply(nim.nests.m$samples, coda::mcmc))
dev.off()
plot(samples_mcmc[[1]][,1])

nim.nests.m$summary$all.chains
summary(cat.mr2)



# END ####

# MODEL ASSESSMENT in BAYES ANALYSES #####

m5 <- nimbleCode({
  
  # priors
  for(k in 1:2){
    B0[k] ~ dnorm(mean = 0, sd = 10)
  }
  B1 ~ dnorm(mean = 0, sd = 10)
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  # likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0[spp[i]] + B1*pred_mass[i]
  }
  
  
  # assess model using a sums-of-squares style discrepancy
  for(i in 1:nObs){
    residual[i] <- y[i] - mu[i]             # Residuals for observed data
    predicted[i] <- mu[i]                   # Predicted values are in the mean of the distribution
    sq.res[i] <- pow(residual[i], 2)        # Squared residuals for observed data
  
  
  # generate replicate data and compute fit statistics
    y.new[i] ~ dnorm(mean = mu[i], sd = sig) # New dataset perfectly adhering to assumptions of model
    sq.new[i] <- pow(y.new[i] - predicted[i], 2) # squared residuals for new data
  }
  
  fit.orig <- sum(sq.res[1:nObs]) # SSE for original dataset
  fit.new <- sum(sq.new[1:nObs])  # SSE for new dataset
  
  # test whether new data are more extreme than observed
  test <- step(fit.new - fit.orig)
  
})


nimData <- list(y = dees)
nimConsts <- list(nObs = length(dees),
                  spp = ifelse(spp == 'BCCH',1,2),
                  pred_mass = pred_mass)

nimInits <- list(B0 = rnorm(2,0,10),
                 B1 = rnorm(1, 0, 10),
                 tau = rgamma(1,1,1))

keepers <- c('B0','B1', 'sig', 'fit.orig', 'fit.new', 'test', 'residual', 'predicted')

nim.dees <- nimbleMCMC(code = m5,
                       constants = nimConsts,
                       data = nimData,
                       inits = nimInits,
                       monitors = keepers,
                       niter = 6000,
                       nburnin = 1000,
                       thin = 1,
                       WAIC = T,
                       nchains = 3,
                       summary = T)

# get samples usable
samples_mcmc <- coda::as.mcmc.list(lapply(nim.dees$samples, coda::mcmc))

samples <- do.call(rbind, samples_mcmc)

# BEFORE LOOKING AT ESTIMATES, two things should be done:
# 1. Assess model convergence
# 2. Assess whether fitted model is adequate for dataset

par(mfrow=c(2,2))
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0\\['))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])
coda::traceplot(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'sig'))])

# Check Rhat
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B0\\['))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'B1'))])
coda::gelman.diag(samples_mcmc[, which(stringr::str_detect(string = colnames(samples), pattern = 'sig'))])


# 2a. Plot residuals
dev.off()
plot(apply(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'predicted'))], 2, mean),
     apply(samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'residual'))], 2, mean), xlab = 'predicted values',
     ylab = 'residuals')
abline(h = 0)


# 2b. Posterior predictive check (graphical and numerical summary)
# pull out our fit measures
fit.new <- samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.new'))]
fit.orig <- samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'fit.orig'))]

# visualize relationship - looking for 1:1 here
# if model fits the data, then about half of the points should be lie above and half 
# below a 1:1 line
plot(fit.orig, fit.new, main = "Graphical posterior predictive check",
xlab = "SSQ for actual data set", ylab = "SSQ for ideal (new) data sets")
abline(0, 1)

# bayesian p-value (numerical summary)
mean(fit.new > fit.orig)

# same as the derived value that we monitored in model
test <- samples[, which(stringr::str_detect(string = colnames(samples), pattern = 'test'))]
mean(test)
