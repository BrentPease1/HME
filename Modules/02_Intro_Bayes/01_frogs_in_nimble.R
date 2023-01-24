library(nimble)

frog_code01 <- nimbleCode({
  
  ## ## ## ## ## ## ## ## 
  #  Prior Probability  # 
  ## ## ## ## ## ## ## ##
  
  # prior probability of presence given what we know about frog
  prior.knowledge <- 0.5        
  
  # probability of detection is relatively high but depends on presence/absence 
  # incorporating prior knowledge of species dectection
  detect.prob <- 0.8 * presence 
  
  ## ## ## ## ## ## ## ## 
  #     Likelihood      # 
  ## ## ## ## ## ## ## ##
  
  # actual presence drawn from a Bernoulli distribution
  presence ~ dbern(prior.knowledge) 
  
  # actual detection occurs with random variation that follows a Bernoulli distribution
  detected ~ dbern(detect.prob)

  
}
)

## frog_code01 <- nimbleCode({
## 
##   ## ## ## ## ## ## ## ##
##   #  Prior Probability  #
##   ## ## ## ## ## ## ## ##
## 
##   # prior probability of presence given what we know about frog
##   prior.knowledge <- 0.5
## 
##   # probability of detection is relatively high but depends on presence/absence
##   detect.prob <- 0.8
## 
##   ## ## ## ## ## ## ## ##
##   #     Likelihood      #
##   ## ## ## ## ## ## ## ##
## 
##   # actual presence drawn from a Bernoulli distribution
##   presence ~ dbern(prior.knowledge)
## 
##   # actual detection occurs with random variation that follows a Bernoulli distribution
##   detected ~ dbern(detect.prob*presence) #conditional on presence of species
## 
## 
##   }
## )

frog_data <- list(detected = 0)

frog_samples01 <- nimbleMCMC(frog_code01,
                            data = frog_data,
                            niter = 10000)

# what did we just create?
str(frog_samples01) #10,000 numbers with a column header of "presence"

summary(frog_samples01) # point estimate nearly identical to solution found analytically 


# with just one chain and a quite simple model, our "model diagnostics" are a bit limited, but we'll do a full inspection once we have a bigger model

frog_code02 <- nimbleCode({
  
  ## ## ## ## ## ## ## ## 
  #  Prior Probability  # 
  ## ## ## ## ## ## ## ##
  
  # prior probability of presence given what we know about the pond
  prior.knowledge <- 0.75        
  
  # probability of detection is relatively high but depends on presence/absence 
  # incorporating prior knowledge of species dectection
  detect.prob <- 0.8 * presence 
  
  ## ## ## ## ## ## ## ## 
  #     Likelihood      # 
  ## ## ## ## ## ## ## ##
  
  # actual presence drawn from a Bernoulli distribution
  presence ~ dbern(prior.knowledge) 
  
  # actual detection occurs with random variation that follows a Bernoulli distribution
  detected ~ dbern(detect.prob)

  }
)

frog_samples02 <- nimbleMCMC(frog_code02,
                            data = frog_data,
                            niter = 10000)

summary(frog_samples02) # point estimate nearly identical to solution found analytically 
