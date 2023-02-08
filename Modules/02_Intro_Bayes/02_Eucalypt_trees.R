## ----setup, include=FALSE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)


## ----setup2, out.width='50%', fig.align='center', fig.cap='', warning = F, message = F---------------------------------
knitr::include_graphics(here::here('Modules/02_Intro_Bayes/slide_supps/tree.jpeg'))


## ----data, echo = T----------------------------------------------------------------------------------------------------
# here is our dataset

tree_diameter <- c(42,43,58,70,47,51,85,63,58,46)


## ----mean/sd, echo = T-------------------------------------------------------------------------------------------------

(data.frame(summary = c('mean', 'SD', 'SE'),
           value = c(round(mean(tree_diameter),1), 
                     round(sd(tree_diameter),1), 
                     round(sd(tree_diameter)/sqrt(length(tree_diameter)),1))))



## ----tree lm, echo = T-------------------------------------------------------------------------------------------------

summary(lm(tree_diameter ~ 1))


## ----model01, echo = T-------------------------------------------------------------------------------------------------
library(nimble)

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
keepers <- c('pop.mean', 'pop.sd') # do we really need to monitor pop.sd?

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
# function to summarize samples
getValues <- function(x){
  mean <- mean(x)
  sd <- sd(x)
  quants <- quantile(x, probs = c(0.025, 0.5, 0.975))
  out <- matrix(data = c(mean, sd, quants), nrow = 1, dimnames = list(1,c('mean',
                                                                          'sd',
                                                                          'lower95',
                                                                          'median',
                                                                          'upper95')))
  return(out)
}
# .......................................................................
# .......................................................................

# First, "Summary" gives us some simple stuff
samples$summary$all.chains

# this is important - let's talk through it.
# prior knowledge said 43-63
# we updated prior knowledge with new knowledge and reduced 95% CI. 
# pop sd is exactly the same value as we passed it - does not get sampled - just a constant


# (HOW DOES 95% CRED INT differ from 95% CONF INT?)
# but also,

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples$samples, coda::mcmc))

# Look at traceplots of the parameters
par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, 1:2])

# calculate Rhat convergence diagnostic of parameters
# "Gelman-Rubin Statitsic" - compares ratio of the variation of the samples within a chain and the variation of samples when the chains are pooled; variation within and between chains should stabilize with convergence (i.e., go to 1)
# rule of thumb is anything <1.1
coda::gelman.diag(samples_mcmc[,1]) # just look at pop.mean = all good

# extract mean and SD

samplesdf <- do.call(rbind, samples_mcmc)
pop.mean <- samplesdf[, 1]
pop.sd <- samplesdf[, 2]
getValues(pop.mean)
getValues(pop.sd)

# MCMCplots is nice too
library(mcmcplots)
mcmcplot(samples$samples, dir = here::here('Modules/02_Intro_Bayes/output'), filename = "tree_model01")


## ----analytically, echo = T--------------------------------------------------------------------------------------------
tree_diameter <- c(42,43,58,70,47,51,85,63,58,46)

posterior_mean <- function(prior_mean, prior_var, data_mean, data_var, n){
  ((prior_mean / prior_var) + (data_mean*(n / data_var))) / ((1/prior_var) + (n/data_var))
}

posterior_var <- function(prior_var, data_var, n){
  (prior_var * (data_var / n)) / ((data_var / n) + prior_var)
}

# posterior mean
print(paste("posterior mean:", round(posterior_mean(prior_mean = 53, prior_var = 5^2, data_mean = mean(tree_diameter), data_var = var(tree_diameter), n = length(tree_diameter)),2)))

# posterior variance
print(paste("posterior variance:", round(posterior_var(prior_var = 5^2, data_var = var(tree_diameter), n = length(tree_diameter)),2)))

# posterior SD, since it is reported
print(paste("posterior SD:", round(sqrt(posterior_var(prior_var = 5^2, data_var = var(tree_diameter), n = length(tree_diameter))),2)))

# posterior 95% CI is same as nimble
print(paste("posterior 95% CI:", 
            round(posterior_mean(prior_mean = 53, prior_var = 5^2, data_mean = mean(tree_diameter), data_var = var(tree_diameter), n = length(tree_diameter)),1) - 
              1.96 * round(sqrt(posterior_var(prior_var = 5^2, data_var = var(tree_diameter), n = length(tree_diameter))),1), "-",
            round(posterior_mean(prior_mean = 53, prior_var = 5^2, data_mean = mean(tree_diameter), data_var = var(tree_diameter), n = length(tree_diameter)),1) + 
              1.96 * round(sqrt(posterior_var(prior_var = 5^2, data_var = var(tree_diameter), n = length(tree_diameter))),1)
              ))


## ----uniform, echo = T-------------------------------------------------------------------------------------------------
tree_model02 <- nimbleCode({
  
  ## Priors ##
  pop.mean ~ dunif(0,200) # the population mean has an equal probability of being any number between 0 and 200. 200 is arbitrary - what if we put 500? 1000?

  pop.sd ~ dunif(0, 100) # the pop.sd has an equal probability of being any number between 0 and 100
  
  # NOTE: before nimble, had to do something like this:
  # pop.sd ~ dunif(0, 100)
  # pop.var <- pop.sd * pop.sd
  # pop.prec <- 1/pop.var # pass pop.prec to dnorm below
  
  # likelihood
  for(i in 1:nObs){
    tree[i] ~ dnorm(pop.mean, sd = pop.sd) 
  }

})


## ----uniform inits, echo = T-------------------------------------------------------------------------------------------
inits <- list(pop.mean = runif(n = 1, min = 0, max = 200),
              pop.sd = runif(n = 1, min = 0, max = 100))


## ----echo = T----------------------------------------------------------------------------------------------------------
tree_data <- list(tree = tree_diameter)
tree_constants <- list(nObs = length(tree_diameter))


## ----echo = T----------------------------------------------------------------------------------------------------------

# things we want `NIMBLE` to keep track of:
# (very useful with complexity)
keepers <- c('pop.mean', 'pop.sd')

# MCMC settings
nc = 3
nb = 1000
ni = nb + 2000
nt = 1




## ----echo = T----------------------------------------------------------------------------------------------------------

# one call
samples02 <- nimbleMCMC(
    code = tree_model02, # changed model name
    constants = tree_constants,
    data = tree_data,
    inits = inits,
    monitors = keepers,
    niter = ni,
    nburnin = nb,
    thin = nt,
    nchains = nc,
    summary = T) # get jags-style summary of posterior




## ----uniform output, echo = T------------------------------------------------------------------------------------------
# First, "Summary" gives us some simple stuff
samples02$summary$all.chains

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples02$samples, coda::mcmc))

# Look at traceplots of the parameters
par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, 1:2])

# calculate Rhat convergence diagnostic of parameters
# "Gelman-Rubin Statitsic" - compares ratio of the variation of the samples within a chain and the variation of samples when the chains are pooled; variation within and between chains should stabilize with convergence (i.e., go to 1)
# rule of thumb is anything <1.1
coda::gelman.diag(samples_mcmc) # we can now look at both mean and sd

# extract mean and SD lambda of each grid cell

samplesdf <- do.call(rbind, samples_mcmc)
pop.mean <- samplesdf[, 1]
pop.sd <- samplesdf[, 2]
getValues(pop.mean)
getValues(pop.sd)

# MCMCplots is nice too
library(mcmcplots)
mcmcplot(samples02$samples, dir = here::here('Modules/02_Intro_Bayes/output'), filename = "tree_model02")


## ----uniform2, echo = T------------------------------------------------------------------------------------------------

tree_model03 <- nimbleCode({
  
  ## Priors ##
  pop.mean ~ dunif(0,1000) # the population mean has an equal probability of being any number between 0 and 200. 200 is arbitrary - what if we put 500? 1000?

  pop.sd ~ dunif(0, 500) # the pop.sd has an equal probability of being any number between 0 and 100
  
  # NOTE: before nimble, had to do something like this:
  # pop.sd ~ dunif(0, 100)
  # pop.var <- pop.sd * pop.sd
  # pop.prec <- 1/pop.var # pass pop.prec to dnorm below
  
  # likelihood
  for(i in 1:nObs){
    tree[i] ~ dnorm(pop.mean, sd = pop.sd) 
  }

})

# don't really need to change initial values - they will fall within acceptable values and will still be radndom
inits <- list(pop.mean = runif(n = 1, min = 0, max = 200),
              pop.sd = runif(n = 1, min = 0, max = 100))

# data and constants
tree_data <- list(tree = tree_diameter)
tree_constants <- list(nObs = length(tree_diameter))


# params to monitor
keepers <- c('pop.mean', 'pop.sd')

# MCMC settings
nc = 3
nb = 1000
ni = nb + 2000
nt = 1


# one call
samples03 <- nimbleMCMC(
    code = tree_model03, # changed model name
    constants = tree_constants,
    data = tree_data,
    inits = inits,
    monitors = keepers,
    niter = ni,
    nburnin = nb,
    thin = nt,
    nchains = nc,
    summary = T) # get jags-style summary of posterior

# First, "Summary" gives us some simple stuff
samples03$summary$all.chains

# this is important - let's talk through it.
# prior knowledge said 43-63
# we updated prior knowledge with new knowledge and reduced 95% CI. 
# pop sd is exactly the same value as we passed it - does not get sampled - just a constant


# (HOW DOES 95% CRED INT differ from 95% CONF INT?)
# but also,

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples03$samples, coda::mcmc))

# Look at traceplots of the parameters
par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, 1:2])

# calculate Rhat convergence diagnostic of parameters
# "Gelman-Rubin Statitsic" - compares ratio of the variation of the samples within a chain and the variation of samples when the chains are pooled; variation within and between chains should stabilize with convergence (i.e., go to 1)
# rule of thumb is anything <1.1
coda::gelman.diag(samples_mcmc) # we can now look at both mean and sd

# extract mean and SD lambda of each grid cell

samplesdf <- do.call(rbind, samples_mcmc)
pop.mean <- samplesdf[, 1]
pop.sd <- samplesdf[, 2]
getValues(pop.mean)
getValues(pop.sd)

# MCMCplots is nice too
library(mcmcplots)
mcmcplot(samples03$samples, dir = here::here('Modules/02_Intro_Bayes/output'), filename = "tree_model03")



## ----normal prior, echo = T--------------------------------------------------------------------------------------------
# model
tree_model04 <- nimbleCode({
  
  ## Priors ##
  pop.mean ~ dnorm(0, sd = 100) # still flat 'uninformative', but draws from a random normal
  
  # to set up the standard deviation, we need to make sure that the values are not negative.
  # If you rememeber, a gamma distribution has support from 0 -> Inf. This is a good choice.
  # we'll set the gamma on the precision term and then transform that to something more natural to us - the std. dev. 
  prec ~ dgamma(0.1, 0.1)
  pop.sd <- 1/sqrt(prec)
  

  # likelihood
  for(i in 1:nObs){
    tree[i] ~ dnorm(pop.mean, sd = pop.sd) 
  }

})

# inits
inits <- list(pop.mean = rnorm(n = 1, 100, 10), # now rnorm
              prec = rgamma(n = 1, 0.1,0.1))

# gather data and constants
tree_data <- list(tree = tree_diameter)
tree_constants <- list(nObs = length(tree_diameter))

# parameters to monitor
keepers <- c('pop.mean', 'pop.sd')

# MCMC settings
nc = 3
nb = 5000
ni = nb + 5000
nt = 1

# get samples
samples04 <- nimbleMCMC(
    code = tree_model04, # changed model name
    constants = tree_constants,
    data = tree_data,
    inits = inits,
    monitors = keepers,
    niter = ni,
    nburnin = nb,
    thin = nt,
    nchains = nc,
    summary = T)


# check it out
samples04$summary$all.chains

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples04$samples, coda::mcmc))

# Look at traceplots of the parameters
par(mfrow=c(1,2))
coda::traceplot(samples_mcmc[, 1:2])

# calculate Rhat convergence diagnostic of parameters
# "Gelman-Rubin Statitsic" - compares ratio of the variation of the samples within a chain and the variation of samples when the chains are pooled; variation within and between chains should stabilize with convergence (i.e., go to 1)
# rule of thumb is anything <1.1
coda::gelman.diag(samples_mcmc) # we can now look at both mean and sd

library(mcmcplots)
mcmcplot(samples04$samples, dir = here::here('Modules/02_Intro_Bayes/output'), filename = "tree_model04")



## ----falcon, out.width='50%', fig.align='center', fig.cap='', warning = F, message = F---------------------------------
knitr::include_graphics(here::here('Modules/02_Intro_Bayes/slide_supps/falcon.jpg'))


## ----peregrines, echo = T----------------------------------------------------------------------------------------------
peregrines <- c(616, 653, 658, 608, 575, 621, 583, 602, 581, 604, 584, 604)


## ----scipt, echo = F---------------------------------------------------------------------------------------------------
# knitr::purl(input = here::here('Modules/02_Intro_Bayes/02_Eucalypt_trees.Rmd'),
#             output = here::here('Modules/02_Intro_Bayes/02_Eucalypt_trees.R'))

