rm(list = ls(all = TRUE))
graphics.off() 
gc() 

# require libraries -------------------------------------------------------
x = c("BRugs", "R2WinBUGS", "R2jags")
lapply(x, require, character.only = TRUE)


# Model in BUGS code ------------------------------------------------------
sink("surplusProduction.txt")
cat("
    model
    {
    # priors K
    K ~ dlnorm(15.92, 10)I(10000, 500000)
    #priors r
    r ~ dnorm(0.6, 500)I(0.01, 1.2)
    
    #prior q
    iq ~ dgamma(0.001, 0.001)#I(0.5, 200) #Non informative prior
    q <- 1/iq
    
    #priors isigma itau
    isigma2 ~ dunif(0.02, 400)       #According to Gelman  
    sigma2 <- 1/isigma2
    
    itau2 ~ dgamma(1.7, 0.01)
    tau2 <- 1/itau2
    
    #time step [1] conditions
    Pmed[1] <- 0
    P[1] ~ dlnorm(Pmed[1], isigma2)T(0.05, 1.6)
    
    #time steps of model
    for(t in 2:N){
    Pmed[t] <- log(max(P[t-1] + (r*P[t-1])*(1-P[t-1]) - C[t-1]/K, 0.001))
    P[t] ~ dlnorm(Pmed[t], isigma2)T(0.05, 1.5)
    }
    
    #Sampling Distribution
    for (t in 1:N)                              
    {
    Imed[t] <- log(q*K*P[t])
    I[t] ~ dlnorm(Imed[t], itau2)
    
    #posterior predictions
    index[t] <- log(q*K*P[t])
    I.new[t] ~ dlnorm(index[t], itau2)
    }
    
    #additional parameters and preditions
    MSP <-  r*K/4
    EMSP <-  r/(2*q)
    BMSY <- K/2
    P2015 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
    B2015 <-  P2015*K
    }
    ",fill=TRUE)
sink()
# End of BUGS model code --------------------------------------------------


# Data --------------------------------------------------------------------

C = c(28025, 29787, 35651, 31456, 37078, 33755, 35333, 49473, 57153, 53359, 43688, 58961, 55830, 50000) #Official catch data FAO
I = c(2.2587839, 1.6809916, 1.1716501, 0.8808423, 0.8579199, 0.8192171, 0.7820640, 0.8267496, 0.7489146,
      0.6116128, 0.6667181, 0.8404752, 0.7942241, 0.9448188) #Standardized catch rate (using hold capacity)
N = length(C)

data.list = list("N" = N, "C" = C, "I" = I)
inits.list = function(){
  list(K = runif(1, 10000, 500000),
       r = runif(1, 0.01, 1.2), 
       P = runif(14, 0.5, 1), 
       iq = runif(1, 0.001, 0.001), 
       isigma2 = 100, 
       itau2 = 100)
}



# Specify which stochastic quantities ("parameters") to parametros --------

parametros = c("r", "K",  "q", "BMSY", "MSP", "sigma2", "tau2", "I.new")

out.jags = jags(data.list,  inits.list, parametros, "surplusProduction.txt", 
                n.chains = 3, n.thin = 10, n.iter = 500000, n.burnin = 1000)

out.jags
