library(unmarked)
library(nimble)
library(here)
library(coda)
library(mcmcplots)
library(ggplot2)


# Choose sample sizes and prepare observed data array y
set.seed(24)                  # So we all get same data set
M <- 100                      # Number of sites
J <- 2                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8                    # Probability of occupancy or presence
p <- 0.5                      # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi)  # R has no Bernoulli

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}

# Look at data
sum(z)                        # True number of occupied sites
sum(apply(y, 1, max))         # Observed number of occupied sites
head(cbind(z=z, y))           # Truth and measurements for first 6 sites


umf <- unmarkedFrameOccu(y = y)  # Create unmarked data frame
summary(umf)                     # Summarize data frame
(fm1 <- occu(~1 ~1, data = umf)) # Fit model

backTransform(fm1, "state")      # Get estimates on probability scale
backTransform(fm1, "det")



# Specify model in BUGS language
m1 <- nimbleCode({
  
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  
  # Likelihood
  for (i in 1:M) {    # Loop over sites
    z[i] ~ dbern(psi) # State model
    for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(mu[i])
    }
     mu[i] <- z[i]*p  
  }
})


# Initial values
zst <- apply(y, 1, max)       # Avoid data/model/inits conflict
inits <- list(z = zst,
              psi = runif(1),
              p = runif(1))

nimData <- list(y = y)
nimConsts <- list(M = M,
                  J = J) 

# Parameters monitored
keepers <- c('psi', 'p')

m1 <- nimbleMCMC(code = m1,
                 data = nimData,
                 constants = nimConsts,
                 monitors = keepers,
                 inits = inits,
                 niter = 25000,
                 nburnin = 5000,
                 thin = 20,
                 nchains = 3,
                 summary = T,
                 samplesAsCodaMCMC = TRUE) #this is also new - can skip some data mgmt steps with this


params_of_interest <- c("psi", "p")
# inspect posterior distributions
mcmcplot(m1$samples[,params_of_interest])

# inspect output
m1$summary

