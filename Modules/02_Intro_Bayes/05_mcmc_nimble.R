library(nimble)

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

## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
## -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --


# Three ways to run an MCMC in nimble

# Each approach gives you different levels of control over how the MCMC algorithm runs. 

# ranked least to most control: 1. nimbleMCMC(), 2. runMCMC(), and 3. the mcmc$run() MCMC object function.

# I think nimble's recommendation is the runMCMC() approach, but nimbleMCMC() is great for straightforward models

# Demonstration of each approach using the peregrine falcon data.

peregrines <- c(616, 653, 658, 608, 575, 621, 583, 602, 581, 604, 584, 604)

# model code
falcon <- nimbleCode({
  
  ## Priors ##
  pop.mean ~ dnorm(611, sd = 50) # generous sd to get weakly informative prior
  
  pop.sd ~ dgamma(10, 0.4) # explored properties of gamma dist to get something that reflects what I know about pop
  
  
  # likelihood
  for(i in 1:nObs){
    falcon[i] ~ dnorm(pop.mean, sd = pop.sd) 
  }
  
})

# Create some constants, data, and initial values to pass to the model builder
constants <- list(nObs = length(peregrines))

data <- list(falcon = peregrines)

inits <- list(pop.mean = rnorm(1, mean = 611, sd = 50),
              pop.sd = rgamma(1, 10, 0.4))

## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##

# approach 1 - nimbleMCMC()

## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##


# "The nimbleMCMC() function handles everything from model building to 
# compilation to running the MCMC. It takes model code and model inputs 
# (constants, data, and initial values) and returns MCMC samples. 
# The nimbleMCMC() function is easy and quick to use but provides little 
# control over the MCMC algorithm. 
# Also, since the objects (the uncompiled and compiled model and MCMC) are 
# used only in the function, they aren’t saved and can’t be manipulated for debugging."

nimbleMCMC_samples <- nimbleMCMC(code = falcon, 
                                 constants = constants, 
                                 data = data, 
                                 inits = inits,
                                 nburnin = 2000, 
                                 niter = 8000,
                                 nchains = 3,
                                 summary = T)

head(nimbleMCMC_samples)

# We can look at the results:

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(nimbleMCMC_samples, coda::mcmc))

# Look at traceplots of the parameters
mcmcplots::mcmcplot(lapply(nimbleMCMC_samples, coda::mcmc))

# check Gelman-Rubin Statistic for convergence
coda::gelman.diag(samples_mcmc[,1:2])

# extract mean and SD
samplesdf <- do.call(rbind, samples_mcmc)
popMean <- getValues(samplesdf[,1])
popSD <- getValues(samplesdf[,2])


## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##

# approach 2 - runMCMC()

## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##

# The runMCMC() function provides more control 
# and a closer look at what NIMBLE is doing. 
# runMCMC() takes a NIMBLE MCMC object, so we need to back up and build these things.

# First, build the NIMBLE model.

falconModel <- nimbleModel(code = falcon,
                           constants = constants,
                           data = data, 
                           inits = inits)


# Next, build an MCMC object for this model with specified settings
# in configureMCMC() using buildMCMC().
conf <- configureMCMC(falconModel, 
                      monitors = c("pop.mean", "pop.sd"))

# note that since we are not really changing any MCMC settings in "conf", 
# it is unncessary, but good to get in habit

falconMCMC <- buildMCMC(falconModel, conf = conf)

# this tells us some stuff that the MCMC is planning on doing

# ===== Monitors =====
#   thin = 1: pop.mean, pop.sd
# ===== Samplers =====
#   RW sampler (1)
# - pop.sd
# conjugate sampler (1)
# - pop.mean

# from here, We can use runMCMC() on this uncompiled model if we want to, 
# but we’ll get a warning reminding us that compiled MCMC will be much faster.
# Uncompiled MCMC will be really slow, but running in R allows easy testing and debugging. 
# Here is how you would run it for 5 iterations:

# runMCMC(falconMCMC, niter = 5)

# To compile the model and MCMC object, use compileNimble(). 
# Note that the model object must be compiled first for the MCMC to compile.

CfalconModel <- compileNimble(falconModel)

CfalconMCMC <- compileNimble(falconMCMC, project = falconModel)

runMCMC_samples <- runMCMC(CfalconMCMC,
                           nburnin = 2000, 
                           niter = 8000,
                           nchains = 3)

# We can look at the results:

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(runMCMC_samples, coda::mcmc))

# Look at traceplots of the parameters
mcmcplots::mcmcplot(lapply(runMCMC_samples, coda::mcmc))

# check Gelman-Rubin Statistic for convergence
coda::gelman.diag(samples_mcmc[,1:2])

# extract mean and SD
samplesdf <- do.call(rbind, samples_mcmc)
(popMean <- getValues(samplesdf[,1]))
(popSD <- getValues(samplesdf[,2]))


## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##

# approach 3 - mcmc$run()

## -- -- -- -- -- -- -- -- ##
## -- -- -- -- -- -- -- -- ##

# You can directly call the $run() method of an MCMC object. 
# Similarly to runMCMC(), the MCMC object and model object must be built and compiled.
# Once this is done, you can use it as:

CfalconMCMC$run(10000)
run_samples <- as.matrix(CfalconMCMC$mvSamples)
plot(run_samples[ , 'pop.mean'], type = 'l', xlab = 'iteration',  ylab = "Population Mean")
plot(run_samples[ , 'pop.sd'], type = 'l', xlab = 'iteration',  ylab = "Population SD")


# There are several advantages to using runMCMC() over $run(). 
# runMCMC() can handle multiple chains and burn-in. 
# It also returns a friendlier output, whereas $run() requires retrieving samples
# from the MCMC object and converting to a matrix.


# using approach 2, we can also have complete control over our MCMC.
# for example, there might be a point when you want to modify the samplers being used in a model
# leaving off at line 110,

falconMCMCconfig <- configureMCMC(falconModel) #note that buildMCMC did this in one step - make default MCMC config and build mcmc

falconMCMCconfig$printSamplers()

# we can see that our pop mean parameter has been assigned a conjugate sampler - this is good - it will be quick and clean
# the pop.sd parameter has been assigned an adaptive random walk Metropolis-Hasting sampler.
# This is the type of sampler we just explored outside of nimble.

# we can do some cool stuff here.
# we can completely change the sampler that is being used - for more on this,
# check out nimble documentation, section 7.2

# we can also just modify basic mcmc settings such as iterations, thinning, etc.,
# *FOR EACH individual parameter*. -> lots of control

# looking back at traceplots:
# Look at traceplots of the parameters
mcmcplots::mcmcplot(lapply(runMCMC_samples, coda::mcmc))

# we had correlation in the pop.sd posterior that reduced our effective sample size for that parameter:
coda::effectiveSize(lapply(runMCMC_samples, coda::mcmc))

# we got back every darn sample from the conjugate sampler but not for our RW sampler on the SD term. 
# correlation had *something* (not everything) to do with this.
# thinning the chain can reduce correlation.

# BUT, no need to thin pop.mean - it is fine as it is. So let's leave it as it is

falconMCMCconfig$monitors <- 'pop.mean'
falconMCMCconfig$monitors2 <- 'pop.sd'

falconMCMCconfig$thin <- 1
falconMCMCconfig$thin2 <- 10 #drawback is you need more iterations to get same chain length

falconMCMC <- buildMCMC(falconMCMCconfig)
CfalconMCMC <- compileNimble(falconMCMC, project = falconModel, resetFunctions = T)

thinned_samples <- runMCMC(CfalconMCMC, niter = 8000, nburnin = 2000, nchains = 3)

# result is that thin/thin2 separates into samples / samples 2
mcmcplots::mcmcplot(lapply(thinned_samples$samples, coda::mcmc))
mcmcplots::mcmcplot(lapply(thinned_samples$samples2, coda::mcmc))

# check ESS
coda::effectiveSize(lapply(thinned_samples$samples, coda::mcmc))  # same
coda::effectiveSize(lapply(thinned_samples$samples2, coda::mcmc)) # fewer, but accepting samples at higher rate

maxS <- 18000/falconMCMCconfig$thin2 #(6000 iterations * 3 chains = 18000 / thinning rate of 10)
# with thinning
coda::effectiveSize(lapply(thinned_samples$samples2, coda::mcmc)) / maxS

# without thinning
coda::effectiveSize(lapply(runMCMC_samples, coda::mcmc))[2] / 18000

# without thinning, we were retaining about 2 out of every 10 samples

# with thinning, retaining 9 out of every 10 samples