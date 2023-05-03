library(AHMbook)

# simulation of a metacommunity
set.seed(73273)
simComm(type="det/nondet", nsite=100, nrep=3, nspec=10,
        mean.psi=0.5, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,  
        mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)


(mean.psi <- rnorm(n = 10, #species
                  mean = plogis(0.5), 
                  sd = 1))
plogis(mean.psi)

# effect of, say, forest cover is positive but there is variation
simComm(type="det/nondet", nsite=100, nrep=3, nspec=10,
        mean.psi=0.5, sig.lpsi=1, mu.beta.lpsi=2, sig.beta.lpsi=1,  
        mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)

# effect of ecological and detection covariates
simComm(type="det/nondet", nsite=100, nrep=3, nspec=10,
        mean.psi=0.5, sig.lpsi=1, mu.beta.lpsi=2, sig.beta.lpsi=1,  
        mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=1, show.plot = TRUE)



# Some possibly interesting settings of the function
data <- simComm(nsite = 267, nspec = 190, mean.psi = 0.25, sig.lpsi = 2,
                mean.p = 0.12, sig.lp = 2) # similar to Swiss BBS


data <- simComm(mean.psi = 1) # all species occur at every site

data <- simComm(mean.p = 1) # no measurement error (perfect detection)


# Effect of spatial sample size (nsite) on species richness in sample (Ntotal.fs)
data <- simComm(nsite=50, nspec = 200) # 1-3 are usually missed in sample
data <- simComm(nsite=30, nspec = 200) # 4-6 usually missed
data <- simComm(nsite=10, nspec = 200) # around 30 typically missed
