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

data("LionNoses") #abd
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
  B0 ~ dnorm(mean = 0, sd = 10)
  B1 ~ dnorm(mean = 0, sd = 10)
  
  tau ~ dgamma(1,1)
  sig <- 1/sqrt(tau)
  
  #likelihood
  for(i in 1:nObs){
    y[i] ~ dnorm(mean = mu[i], sd = sig)
    mu[i] <- B0 + B1 * pct.black[i]
  }
  
})



nimData <- list(y = LionNoses$age)
nimConsts <- list(nObs = length(LionNoses$age),
                  pct.black = LionNoses$percentage.black)

nimInits <- list(B0 = rnorm(n = 1, mean = 0, sd = 10),
                 B1 = rnorm(n = 1, mean = 0, sd = 10),
                 tau = rgamma(n = 1, 1, 1))

keepers <-  c('B0', 'B1', 'sig')

nim.noses <- nimbleMCMC(code = m1, 
                        data = nimData,
                        constants = nimConsts,
                        inits = nimInits,
                        monitors = keepers,
                        niter = 6000,
                        nburnin = 1000,
                        thin = 1,
                        nchains = 3,
                        summary = T)
                        

# Don't be bayesic, check your traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.noses$samples, coda::mcmc))

par(mfrow=c(1,3))
coda::traceplot(samples_mcmc[, 1:3])

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


# traceplots
samples_mcmc <- coda::as.mcmc.list(lapply(nim.jaws$samples, coda::mcmc))

par(mfrow=c(2,2))
coda::traceplot(samples_mcmc[, 1:4])

# Check Rhat
coda::gelman.diag(samples_mcmc[, c(1,3:4)]) #check column index

# Checks out, check parameter estimates
nim.jaws$summary$all.chains

# compare to freq
summary(lm.cat)



## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# means parameterization
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# the model coefficients represent the expected value of the response variable 
# for each level of the predictor variable

lm.cat2 <- lm(jaws ~ sex-1, data = jack)
summary(lm.cat2)

# NIM


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



# CHICKADEES - means parameterization ####


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


# traceplots
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

