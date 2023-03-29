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

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



nmix1 <- nimbleCode({

  # Priors
  lambda ~ dgamma(0.01,0.01)
  p ~ dbeta(1,1)
  
  # Likelihood
  for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
      C[i,j] ~ dbin(p, N[i]) # Observation model
    }
  }
  
  # derived values of interest
  Ntotal <- sum(N[1:M])
  
  # Posterior predictive distributions of chi2 discrepancy
  for (i in 1:M) {
    for (j in 1:J) {
      C.sim[i,j] ~ dbin(p, N[i]) # Create new data set under model
      e.count[i,j] <- N[i] * p   # Expected Count
      
      # Add small value e to denominator to avoid division by zero (e = 0.001)
      # Chi-square discrepancy for the actual data
      chi2.actual[i,j] <- pow((C[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      # Chi-square discrepancy for the simulated ('perfect') data
      chi2.sim[i,j] <- pow((C.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
      
    }
  }
  # Add up individual chi2 values for overall fit statistic
  fit.actual <- sum(chi2.actual[1:M,1:J])  # Fit statistic for actual data set
  fit.sim <- sum(chi2.sim[1:M,1:J])        # Fit statistic for a fitting model
  c.hat <- fit.actual / fit.sim      # Over-dispersion measure (c-hat) - happens when variance is >> than mean 
  bpv <- step(fit.sim-fit.actual)    # Bayesian p-value

  
})


nimData <- list(C = C)
nimConsts <- list(M = M,
                  J = J,
                  e = 0.001) 

# THIS IS NEW
# N is a latent state variable - there is some true but unknown N at each site
# for any variable like this, we need to supply starting values
# based on observed data

Nst <- apply(C,1,max)

nimInits <- list(N = Nst,
                 lambda = rgamma(1,0.1, 0.1),
                 p = rbeta(1,1,1))

keepers <- c('lambda', 'p', 'Ntotal', 'c.hat', 'bpv')

m2 <- nimbleMCMC(code = nmix1,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 inits = nimInits,
                 niter = 25000,
                 nburnin = 5000,
                 thin = 20,
                 nchains = 3,
                 summary = T,
                 samplesAsCodaMCMC = TRUE)


params_of_interest <- c("lambda", "p")
# inspect posterior distributions
mcmcplot(m2$samples[,params_of_interest])

# inspect output
m2$summary$all.chains
