library(unmarked)
library(nimble)
library(here)
library(coda)
library(mcmcplots)

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# simulate simplest dataset 

set.seed(277)
M <- 100
J <- 2

# true parameter values
lambda <- 3 # expected abundance
det.prob <- 0.4 #prob of detection per individual

# matrix to hold count data
C <- matrix(data = NA, nrow = M, ncol = J)

# simulate true abundance at each site
# this is the thing we care about, but do not perfectly observe
# the "ecological process" here, is lambda - lambda controls nature
N <- rpois(n = M, lambda = lambda)

# simulate observation process (counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = det.prob)
}

# Take a look at truth
table(N)

# total number of individuals in population
sum(N) #why this number?

# Number of occupied sites
sum(N > 0)

# mean abundance
mean(N)

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Here is our observation of truth
table(apply(C,1,max))

# observed population size
sum(apply(C,1,max))

# Observed number of occupied sites
sum(apply(C,1,max) >0)

# observed mean abundance
mean(apply(C,1,max))

# comparison of truth vs observations
head(cbind(N=N, count1=C[,1], count2=C[,2]))

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# FIT MODELS

# unmarked is workhorse for a range of models
# https://cran.r-project.org/web/packages/unmarked/vignettes/unmarked.html

# fit n-mixture models using pcount function

# create data frame that library can ingest
umf <- unmarkedFramePCount(y = C)
summary(umf)

# fit model, get estimates on link scale
# specify linear model of each submodel
# observation is first, ecological is second
m1 <- pcount(~1 ~1, data = umf)

# get estimates on natural scale
backTransform(m1, "state")
backTransform(m1, "det")

# Do estimates reflect truth?

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


# N-mix in Nimble

nmix1 <- nimbleCode({
  
  # Priors
  lambda ~ dgamma(0.01,0.01)
  p ~ dbeta(1,1)
  
  # likelihood
  for(i in 1:M){
    N[i] ~ dpois(lambda)           # state model
    for(j in 1:J){
      C[i,j] ~ dbin(p, N[i])       # observation model
    }
  }
  
  # derived values of interest
  Ntotal <- sum(N[1:M])
})


nimData <- list(C = C)
nimConsts <- list(M = M,
                  J = J) 

# THIS IS NEW
# N is a latent state variable - there is some true but unknown N at each site
# for any variable like this, we need to supply starting values
# based on observed data

Nst <- apply(C,1,max)

nimInits <- list(N = Nst,
                 lambda = rgamma(1, 0.01, 0.01),
                 p = rbeta(1,1,1))

keepers <- c('lambda', 'p', 'Ntotal')

m2 <- nimbleMCMC(code = nmix1,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 inits = nimInits,
                 niter = 50000,
                 nburnin = 10000,
                 thin = 20,
                 nchains = 3,
                 summary = T,
                 samplesAsCodaMCMC = TRUE) #this is also new - can skip some data mgmt steps with this


params_of_interest <- c("lambda", "p")
# inspect posterior distributions
mcmcplot(m2$samples[,params_of_interest])

# inspect output
m2$summary



