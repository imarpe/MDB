rm(list = ls(all = TRUE))
graphics.off() 
gc() 

library(BRugs)
library(R2WinBUGS)
library(R2jags)

#model in BUGS code
sink("surplusProduction.txt")
cat("
    model
    {
    # priors K
    K ~ dlnorm(11.92, 10)I(10000, 200000)
    
    #priors r
    r ~ dnorm(0.84, 0.13)I(0.01, 1.2)
    
    #prior q
    iq ~ dlnorm(-9, 16)#I(0.0001, 0.01)
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
    P2015 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
    B2015 <-  P2015*K
    }
    ",fill=TRUE)
sink()
#end of BUGS model code

#No Oficial data
C = c()
I = c()
N = 15

data.list = list("N" = N, "C" = C, "I" = I)
inits.list = function(){
  list(K = runif(1, 10000, 200000),
       r = runif(1, 0.01, 1.2), 
       P = runif(15, 0.5, 1), 
       iq = runif(1, 0.1, 0.1), 
       isigma2 = runif(1, 80, 150), 
       itau2 = runif(1, 90, 120))
}


# specify which stochastic quantities ("parameters") to parametros
parametros = c("r", "K", "MSP", "EMSP", "B2015", "sigma2", "tau2", "q")

out.jags = jags(data.list,  inits.list, parametros, "surplusProduction.txt", 
                n.chains = 3, n.thin = 10, n.iter = 100000, n.burnin = 1000)
