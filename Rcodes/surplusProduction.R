#Rcode for a Schaefer Production Model

#Data_section
perico.data = read.table("perico.data", header = TRUE)
yr          = perico.data$year
ct          = perico.data$catch
yt1         = perico.data$catch/perico.data$nTrips 
yt2         = perico.data$catch/perico.data$nHooks 
yt3         = perico.data$catch/perico.data$hours 
yt4         = perico.data$catch/perico.data$Ehours 
yt5         = perico.data$catch/perico.data$nHooks_Ehours*1000 

##Parameter_section
theta = c(k=log(71000), r=log(0.9), q = log(0.00000026), sigma=log(0.31))  #Review de inputs

##Procedure_section
schaefer = function(theta)
{
  bt=vector()
  with(as.list(theta),{
    r=exp(r); k=exp(k); q=exp(q); sigma=exp(sigma)

    bt[1]=k
    for(i in 1:length(yr))
    {
      bt[i+1]=bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i]
    }
    epsilon = log(yt)-log(q*bt[1:length(yt)])
    nloglike = -sum(dnorm(epsilon,0,sigma,log=T))
    
    return(list(nloglike=nloglike,bt=bt,epsilon=epsilon,it=yt/exp(epsilon)))
  })
}

solver = function(pars,fn,hess=FALSE)
{
  fit = optim(pars,fn,method="BFGS",hessian=hess)
  if(hess){	fit$V = solve(fit$hessian)  			#Variance-Covariance
  fit$S = sqrt(diag(fit$V))				  #Standard deviations
  fit$R = fit$V/(fit$S %o% fit$S) }	#Parameter correlation
  return(fit)
}


## Main_section
fit = solver(theta,fn=function(theta) schaefer(theta)$nloglike, hess=T)
A = schaefer(fit$par)  #CPUE, Biomass, etc

MSY = as.numeric((exp(fit$par)[1]*exp(fit$par)[2])/4) #Maximum sustainable yield
MSY
Fmsy = as.numeric(exp(fit$par)[2]/2)                  #Fishing mortality rate at MSY
Fmsy
Bmsy = as.numeric(exp(fit$par)[1]/2)                  #Biomass giving maximum sustainable yield 
Bmsy
Emsy = as.numeric(exp(fit$par)[2]/2*exp(fit$par)[3])  #Effort that should lead to MSY
