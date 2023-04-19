library(MASS)
library(raster)
library(nimble)
library(here)
library(spdep)
library(fasterize)
library(verification)
library(unmarked)
library(ggplot2)
library(dplyr)

# we're all the same
set.seed(73273)

# set sample sizes
M <- 20^2 # sites
J <- 5    # surveys


# Set up a square lattice region
simgrid <- expand.grid(seq(-1, 1, length.out = sqrt(M)), seq(-1, 1, length.out = sqrt(M)))
simgrid <- simgrid[order(-simgrid[, 2], simgrid[, 1]), ]
simgrid <- cbind(simgrid, 1:nrow(simgrid))

# Set up distance matrix
distance <- as.matrix(dist(simgrid[, 1:2]))

#simulate the covariates
index <- 1:M
phi  <-  runif(1, 0, 0.5)
forest <- MASS::mvrnorm(1, rep(0, M), exp(-phi*distance))

# create raster to interact with covariates and hold other simulated data
ras <- raster(nrow = sqrt(M), ncol = sqrt(M), xmn = -1, xmx = 1, ymn = -1, ymx = 1)
ras$index <- index
ras$forest <- forest

# covariate effects
B0 <- -1 # logit-scale intercept - what is baseline probability of occupancy?
B1 <- 0.7 # logit-scale effect of forest cover

#truth
ras$psi <- plogis(B0 + B1 * forest)

# generate true occurrence
ras$Z <- rbinom(n = M, size = 1, prob = raster::values(ras$psi))

# generate imperfect counts of the individuals
# Create a detection covariate called distance from road and site elev (line of sight)
set.seed(277) #unnessary, just messing around
dist <- array(runif(M, -1, 1), dim = c(M,J))
set.seed(0817)
elev <- array(runif(M, -1, 1), dim = c(M,J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- 0.5 # Logit-scale intercept
alpha1 <- -3 # Logit-scale slope for wind
alpha2 <- -1 # logit-scale slope for slope
p <- plogis(alpha0 + alpha1 * dist + alpha2 * elev) # Detection probability

# # Take J [ 3 presence/absence measurements at each site
y <- matrix(NA, nrow = M, ncol = J)
for(j in 1:J) {
  y[,j] <- rbinom(nrow(simgrid), 1, p[,j]*raster::values(ras$Z))
}

# check raster
# add observations then plot
ras$y <- apply(y,1,max)
plot(ras)
mapview::mapview(ras)

# Learn about dataset
# TRUTH
# total number of occupied sites
sum(raster::values(ras$Z))

# observed occupied sites
sum(apply(y,1,max))


# bring covariates and observations together 
dat <- data.frame(coords = raster::coordinates(ras), raster::values(ras), y, dist, elev)
names(dat)[8:22] <- c(paste0('visit_', 1:5), paste0('dist_',1:5), paste0('elev_', 1:5))
head(dat)

# --- --- --- --- end simulation data prep --- --- --- --- 


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# --- check truth inside of unmarked --- 


# also shows how to bring in covariates into unmarked
umf <- unmarkedFrameOccu(y = y,
                         siteCovs = data.frame(forest = dat[, c('forest')]), # talk about this
                         obsCovs = list(dist = dist, elev = elev))

summary(umf)                     # Summarize UMF
(fm1 <- occu(~1 ~1, data = umf)) # Fit model
(fm2 <- occu(~dist + elev ~ forest, data = umf)) # Fit model

# reminder
c(B0, B1, alpha0, alpha1, alpha2)

# get on prob scale
backTransform(fm2, 'state') # doesn't work
lc <- linearComb(fm2, c(1, 0), type="state") # Estimate occupancy on the log scale when forest=0
backTransform(lc)                           
# truth
plogis(B0)


lc <- linearComb(fm2, c(1, max(forest)), type="state") # Estimate occupancy on the logit scale when forest= max(forest)
(lc_fm2 <- backTransform(lc))

# get 95% CI and can do this for a range of predictor values
preds <- predict(fm2, 'state', newdata = data.frame(forest = max(forest)))
preds$model <- "full_data"
# truth
plogis(B0 + B1*max(forest))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# OKAY, let's explore what happens when we induce reality


# Visit X percent of the grid
# 20% = 400*.2 = 80 sites
sites <- sample(1:M, M * 0.20)

#visualize
plot(ras$Z, main = "True Occurrence with sampling sites")
points(dat[sites,1:2], pch = 19, cex=1)


umf_sampled <- unmarkedFrameOccu(y = y[sites,],
                         siteCovs = data.frame(forest = dat[sites, c('forest')]),
                         obsCovs = list(dist = dist[sites,], elev = elev[sites,]))
summary(umf)
(fm3 <- occu(~dist + elev ~ forest, data = umf_sampled)) # Fit model
lc <- linearComb(fm3, c(1, max(forest)), type="state") # Estimate occupancy on the log scale when forest=0
(lc_fm3 <- backTransform(lc))

# get 95% CI and can do this for a range of predictor values
fm3_pred <- predict(fm3, 'state', newdata = data.frame(forest = max(forest)))
fm3_pred$model <- 'sampled'
preds <- rbind(preds, fm3_pred)

# Create plot
ggplot(preds, aes(x = model, y = Predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = plogis(B0 + B1*max(forest)), linewidth = 1, col = 'red') +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Model", y = "Predicted") +
  ggtitle("Predicted Values with 95% Confidence Intervals") +
  theme_minimal()

#

# Was that a fluke? 
holder <- data.frame(n_sites = length(sites), prop_sites = length(sites)/M, n_surveys = ncol(y), rep = 0,
                     forest_est = 0, forest_low = 0, forest_high = 0,
                     psi_mean = 0, psi_low = 0, psi_high = 0,
                     psi_max_mean = 0, psi_max_low = 0, psi_max_high = 0)

nsim <- 100
for(i in 1:nsim){
  set.seed(73273*i)
  sites <- sample(1:nrow(simgrid), nrow(simgrid) * 0.20)
  umf_sampled <- unmarkedFrameOccu(y = y[sites,],
                                   siteCovs = data.frame(forest = dat[sites, c('forest')]),
                                   obsCovs = list(dist = dist[sites,], elev = elev[sites,]))
  
  fmInf <- occu(~dist + elev ~ forest, data = umf_sampled)
  con <- predict(fmInf, data.frame(forest = 0), type = 'state')
  con_max <- predict(fmInf, data.frame(forest = max(forest)), type = 'state')
  
  holder[i, 'n_sites'] <- length(sites)
  holder[i, 'prop_sites'] <- length(sites)/M
  holder[i, 'n_surveys'] <- ncol(y)
  
  holder[i, 'forest_est'] <- coef(fmInf)[2]
  holder[i, 'forest_low'] <- confint(fmInf, type = 'state')[2,1]
  holder[i, 'forest_high'] <- confint(fmInf, type = 'state')[2,2]
  
  holder[i, 'psi_mean'] <- con[1]
  holder[i, 'psi_low'] <- con[3]
  holder[i, 'psi_high'] <- con[4]
  
  holder[i, 'psi_max_mean'] <- con_max[1]
  holder[i, 'psi_max_low'] <- con_max[3]
  holder[i, 'psi_max_high'] <- con_max[4]
  
  holder[i, 'rep'] <- i
  
}


ggplot(holder, aes(rep, forest_est)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = forest_low, ymax = forest_high)) +
  geom_hline(yintercept = B1, linewidth = 1, col = 'red') +
  labs(x = "Replicate", y = "Regression Coefficient") +
  ggtitle("Estimated Effect of Forest Cover") +
  theme_minimal()

# how many times did 95% CI cover true effect of forest cover
length(which(holder$forest_low <= B1 & B1 <= holder$forest_high)) / nsim


# est of occupancy at forest=0
ggplot(holder, aes(rep, psi_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_low, ymax = psi_high)) +
  geom_hline(yintercept = plogis(B0), linewidth = 1, col = 'red')+
  labs(x = "Replicate", y = "Predicted Occurence") +
  ggtitle("Estimated Occupancy Prob when forest = 0") +
  theme_minimal()

# est of occupancy at forest=max(forest)
ggplot(holder, aes(rep, psi_max_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_max_low, ymax = psi_max_high)) +
  geom_hline(yintercept = plogis(B0 + B1*max(forest)), linewidth = 1, col = 'red')+
  labs(x = "Replicate", y = "Predicted Occurence") +
  ggtitle("Estimated Occupancy Prob when forest = max(forest") +
  theme_minimal()

# how many times did 95% CI cover true effect of forest cover
length(which(holder$psi_max_low <= plogis(B0 + B1*max(forest)) & plogis(B0 + B1*max(forest)) <= holder$psi_max_high)) / nsim

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# varying sampling effort ####

# OKAY, let's try different sampling efforts
effort <- c(0.2, 0.5, 0.75)


holder <- data.frame(n_sites = length(sites), prop_sites = length(sites)/M, n_surveys = ncol(y), rep = 0,
                     forest_est = 0, forest_low = 0, forest_high = 0,
                     psi_mean = 0, psi_low = 0, psi_high = 0,
                     psi_max_mean = 0, psi_max_low = 0, psi_max_high = 0)

nsim <- 100
counter = 0
for(e in 1:length(effort)){
  for(i in 1:nsim){
    counter = counter + 1
    set.seed(73273*i*e)
    sites <- sample(1:nrow(simgrid), nrow(simgrid) * effort[e])
    umf_sampled <- unmarkedFrameOccu(y = y[sites,],
                                     siteCovs = data.frame(forest = dat[sites, c('forest')]),
                                     obsCovs = list(dist = dist[sites,], elev = elev[sites,]))
    
    fmInf <- occu(~dist + elev ~ forest, data = umf_sampled)
    con <- predict(fmInf, data.frame(forest = 0), type = 'state')
    con_max <- predict(fmInf, data.frame(forest = max(forest)), type = 'state')
    
    holder[counter, 'n_sites'] <- length(sites)
    holder[counter, 'prop_sites'] <- length(sites)/M
    holder[counter, 'n_surveys'] <- ncol(y)
    
    holder[counter, 'forest_est'] <- coef(fmInf)[2]
    holder[counter, 'forest_low'] <- confint(fmInf, type = 'state')[2,1]
    holder[counter, 'forest_high'] <- confint(fmInf, type = 'state')[2,2]
    
    holder[counter, 'psi_mean'] <- con[1]
    holder[counter, 'psi_low'] <- con[3]
    holder[counter, 'psi_high'] <- con[4]
    
    holder[counter, 'psi_max_mean'] <- con_max[1]
    holder[counter, 'psi_max_low'] <- con_max[3]
    holder[counter, 'psi_max_high'] <- con_max[4]
    
    holder[counter, 'rep'] <- i
    
  }
  cat('completed effort', effort[e], '\n')
}

# est of effect of forest
ggplot(holder, aes(rep, forest_est)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = forest_low, ymax = forest_high)) +
  geom_hline(yintercept = B1, linewidth = 1, col = 'red') +
  facet_wrap(~n_sites)
# est of occupancy at forest=0
ggplot(holder, aes(rep, psi_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_low, ymax = psi_high)) +
  geom_hline(yintercept = plogis(B0), linewidth = 1, col = 'red') +
  facet_wrap(~n_sites)
# est of occupancy at forest=max(forest)
ggplot(holder, aes(rep, psi_max_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_max_low, ymax = psi_max_high)) +
  geom_hline(yintercept = plogis(B0 + B1*max(forest)), linewidth = 1, col = 'red') +
  facet_wrap(~n_sites)



# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# varying occu and det probs
# What about letting other things vary, too?
B0_values <- c(0:2) #expected occu on logit-scale with forest value of 0
alpha0_values <- c(-1:1) #expected det prob on logit-scale with wind and site elev at 0
effort <- c(0.2, 0.5)


holder <- data.frame(n_sites = length(sites), prop_sites = length(sites)/M, n_surveys = ncol(y), rep = 0,
                     b0 = 0, alpha = 0,         # this is new                                                          
                     forest_est = 0, forest_low = 0, forest_high = 0,
                     psi_mean = 0, psi_low = 0, psi_high = 0,
                     psi_max_mean = 0, psi_max_low = 0, psi_max_high = 0)

nsim <- 100
counter = 0
for(b in 1:length(B0_values)){
  set.seed(277*b)
  # covariate effects
  B0 <- B0_values[b]
  B1 <- 0.7 # logit-scale effect of forest cover
  
  #truth
  psi <- plogis(B0 + B1 * forest)
  
  # generate true occurrence
  Z <- rbinom(n = M, size = 1, prob = psi)
  
  for(a in 1:length(alpha0_values)){
    set.seed(0817*a)
    
    alpha0 <- alpha0_values[a]
    alpha1 <- -3
    alpha2 <- -1 
    p <- plogis(alpha0 + alpha1 * dist + alpha2 * elev) 
    
    # # Take J [ 3 presence/absence measurements at each site
    y <- matrix(NA, nrow = M, ncol = J)
    for(j in 1:J) {
      y[,j] <- rbinom(nrow(simgrid), 1, p[,j]*Z)
    }

    for(e in 1:length(effort)){
      for(i in 1:nsim){
        counter = counter + 1
        set.seed(73273*i*e)
        sites <- sample(1:nrow(simgrid), nrow(simgrid) * effort[e])
        umf_sampled <- unmarkedFrameOccu(y = y[sites,],
                                         siteCovs = data.frame(forest = forest[sites]),
                                         obsCovs = list(dist = dist[sites,], elev = elev[sites,]))
        
        fmInf <- occu(~dist + elev ~ forest, data = umf_sampled)
        con <- predict(fmInf, data.frame(forest = 0), type = 'state')
        con_max <- predict(fmInf, data.frame(forest = max(forest)), type = 'state')
        
        holder[counter, 'n_sites'] <- length(sites)
        holder[counter, 'prop_sites'] <- length(sites)/M
        holder[counter, 'n_surveys'] <- ncol(y)
        holder[counter, 'b0'] <- b
        holder[counter, 'alpha'] <- a
        
        
        
        holder[counter, 'forest_est'] <- coef(fmInf)[2]
        holder[counter, 'forest_low'] <- confint(fmInf, type = 'state')[2,1]
        holder[counter, 'forest_high'] <- confint(fmInf, type = 'state')[2,2]
        
        holder[counter, 'psi_mean'] <- con[1]
        holder[counter, 'psi_low'] <- con[3]
        holder[counter, 'psi_high'] <- con[4]
        
        holder[counter, 'psi_max_mean'] <- con_max[1]
        holder[counter, 'psi_max_low'] <- con_max[3]
        holder[counter, 'psi_max_high'] <- con_max[4]
        
        holder[counter, 'rep'] <- i
        
      }
      cat('completed effort', effort[e], 'with baseline occu prob at', round(plogis(B0_values[b]),2),
          'and baseline detection prob at', round(plogis(alpha0_values[a]),2), '\n')
    }
  }
}
# stuff for plotting
holder$group <- paste(holder$n_sites, holder$b0, holder$alpha, sep = '_')
intercepts <- holder %>%
  group_by(group) %>%
  summarize(plogis(b0)) %>%
  filter(!duplicated(group))

# talk about this
holder <- holder %>%
  filter(forest_high < 10)

# est of effect of forest
ggplot(holder, aes(rep, forest_est)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = forest_low, ymax = forest_high)) +
  geom_hline(yintercept = B1, linewidth = 1, col = 'red') +
  facet_wrap(~n_sites + b0 + alpha, scales = 'fixed')
# est of occupancy at forest=0
ggplot(holder, aes(rep, psi_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_low, ymax = psi_high)) +
  geom_hline(data = intercepts, aes(yintercept = `plogis(b0)`), linewidth = 1, col = 'red') +
  facet_wrap(~group)
# est of occupancy at forest=max(forest)
ggplot(holder, aes(rep, psi_max_mean)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = psi_max_low, ymax = psi_max_high)) +
  geom_hline(data = intercepts, aes(yintercept = `plogis(b0)`), linewidth = 1, col = 'red') +
  facet_wrap(~group)



