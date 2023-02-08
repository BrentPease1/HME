## ----model01, echo = T-------------------------------------------------------------------------------------------------
library(nimble)

tree_diameter <- c(42,43,58,70,47,51,85,63,58,46)


tree_model01 <- nimbleCode({
  
  
  # Note: historically had to specify "precision" in BUGS/JAGS
  # precision = 1/pop variance
  # pop variance = pop sd*pop sd
  # In Nimble, we can just specify SD, but need to be explicit in likelihood distribution
  
  # Priors
  pop.mean ~ dnorm(53, sd = 5) # need the "sd =" or have to provide precision (1/(5*5) = 0.04)
  
  # need to specify a prior for the standard deviation of this sample.
  # we could approach it several ways, but first we'll use the SD of the actual observations (we'll check out other options in subsequent models)
  pop.sd <- tree_sd # we'll pass `tree_sd` in as data
  
  # likelihood
  for(i in 1:nObs){
    tree[i] ~ dnorm(pop.mean, sd = pop.sd) 
  }
  
})


## ----bundle, echo = T--------------------------------------------------------------------------------------------------
tree_data <- list(tree = tree_diameter,
                  tree_sd = sd(tree_diameter))
tree_constants <- list(nObs = length(tree_diameter))


## ----inits, echo = T---------------------------------------------------------------------------------------------------
inits <- list(pop.mean = rnorm(n = 1, 53, 5))


## ----mcmc, echo = T----------------------------------------------------------------------------------------------------

# things we want `NIMBLE` to keep track of:
# (very useful with complexity)
keepers <- c('pop.mean') # do we really need to monitor pop.sd?

# MCMC settings
nc = 3 # why chains
nb = 1000 # Why burn-ins
ni = nb + 2000 # why inits
nt = 1 # why thinning




## ----samples01, echo = T-----------------------------------------------------------------------------------------------

# one call
samples <- nimbleMCMC(
  code = tree_model01,
  constants = tree_constants,
  data = tree_data,
  inits = inits,
  monitors = keepers,
  niter = ni,
  nburnin = nb,
  thin = nt,
  nchains = nc,
  summary = T) # get jags-style summary of posterior




## ----inspect, echo = T-------------------------------------------------------------------------------------------------

# .......................................................................
# .......................................................................

# First, "Summary" gives us some simple stuff
samples$summary$all.chains