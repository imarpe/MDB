set.seed(19880703)

filename = "MSYdorado.csv"
outfile  = "CatchMSY_OutputDorado.csv"
cdat = read.csv(filename, header=T)
id = "dorado" 

## Asumiendo que hay un solo stock (Posibilidad de correr con varios stocks o unir los stocks)
for(stock in id) {
  yr   = cdat$yr[as.character(cdat$stock) == id]
  ct   = as.numeric(cdat$ct[as.character(cdat$stock) == id])  
  res  = unique(as.character(cdat$res[as.character(cdat$stock) == id])) 
  nyr  = length(yr)
  
  cat("\n","Stock",stock,"\n")
  flush.console()
  
  ## PARAMETER SECTION
  
  ## If resilience is to be used, delete ## in rows 1-4 below and set ## in row 5	below
  start_r  = if(res == "Very low"){c(0.015, 0.1)}
  else if(res == "Low") {c(0.05, 0.5)}
  else if(res == "High") {c(0.4, 1.5)}
  else {c(0.2,1)} ## Medium, o por defecto si "res" no está en la base 
  ## start_r     = c(0.5,1.5)  ## no usar esta línea si hay datos de resilencia
  start_k     = c(max(ct),20*max(ct)) 
  ## startbio 	= c(0.8,1)  
  startbio    = if(ct[1]/max(ct) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} 
  interyr 	= yr[2]   ## 
  interbio 	= c(0, 1) ## 
  ## finalbio 	= c(0.8, 0.9) ## biomass range after last catches, as fraction of k
  finalbio    = if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)}  
  n           = 100000     ## número de iteraciones
  sigR        = 0      ## process error; 0 = modelo determinístico; 0.05? 0.2 valor muy alto.
  
  startbt     = seq(startbio[1], startbio[2], by = 0.05) 	
  parbound = list(r = start_r, k = start_k, lambda = finalbio, sigR)
  
  cat("Last year =",max(yr),", last catch =",ct[nyr],"\n")
  cat("Resilience =",res,"\n")
  cat("Process error =", sigR,"\n")
  cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
  cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
  cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
  cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")
  cat("Initial bounds for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
  
  flush.console()
  
  ## FUNCTIONS
  .schaefer	= function(theta)
  {
    with(as.list(theta), {
      bt=vector()
      ell = 0 
      for (j in startbt)
      {
        if(ell == 0) 
        {
          bt[1]=j*k*exp(rnorm(1,0, sigR))  
          for(i in 1:nyr) 
          {
            xt=rnorm(1,0, sigR)
            bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) 
          }
          
          #Bernoulli likelihood, asigna 0 o 1 para cada combinación de r y k
          ell = 0
          if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
            ell = 1
        }	
      }
      return(list(ell=ell))
      
    })
  }
  
  sraMSY	=function(theta, N)
  {
    with(as.list(theta), 
         {
           ri = exp(runif(N, log(r[1]), log(r[2])))  
           ki = exp(runif(N, log(k[1]), log(k[2])))  
           itheta = cbind(r = ri, k = ki, lam1 = lambda[1], lam2 = lambda[2], sigR = sigR) 
           M = apply(itheta, 1, FUN = .schaefer)
           i = 1:N
           ## función objetivo
           get.ell = function(i) M[[i]]$ell
           ell = sapply(i, get.ell) 
           return(list(r = ri,k = ki, ell = ell))	
         })
  }
  
  ## 
  R1 = sraMSY(parbound, n)  
  
  ## Estadísticos  de r, k, MSY e intervalos
  r1 	= R1$r[R1$ell==1]
  k1 	= R1$k[R1$ell==1]
  msy1  = r1*k1/4
  mean_msy1 = exp(mean(log(msy1))) 
  max_k1a  = min(k1[r1<1.1*parbound$r[1]]) 
  max_k1b  = max(k1[r1*k1/4<mean_msy1]) 
  max_k1 = if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
  
  if(length(r1)<10) {
    cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
    flush.console()
  }
  
  if(length(r1)>=10) {
    
    parbound$r[2] = 1.2*max(r1)
    parbound$k 	  = c(0.9 * min(k1), max_k1)
    
    cat("First MSY =", format(mean_msy1, digits=3),"\n")
    cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
    cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")	
    cat("New range for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
    
    ## Repetir el análisis con nuenos límites de r-k 
    R1 = sraMSY(parbound, n)
    
    ## Estadísticos de r, k y MSY
    r = R1$r[R1$ell==1]
    k = R1$k[R1$ell==1]
    msy = r * k / 4
    mean_ln_msy = mean(log(msy))
    bmsy = k/2
    mean_ln_bmsy = mean(log(bmsy))
    
    ## Plot MSY vs Captura
    par(mfcol=c(2,3))
    plot(yr, ct, type="l", ylim = c(0, 1.1*max(ct)), xlab = "Year", ylab = "Catch (t)", main = stock, lwd = 2)
    abline(h = exp(mean(log(msy))), col = "red", lwd = 2)
    abline(h = exp(mean_ln_msy - 2 * sd(log(msy))), col="red")
    abline(h = exp(mean_ln_msy + 2 * sd(log(msy))), col="red")
    
    hist(r, freq = F, xlim=c(0, 1.2 * max(r)), main = "", col = "gray62")
    abline(v=exp(mean(log(r))),col="red",lwd=2)
    abline(v=exp(mean(log(r))-2*sd(log(r))),col="red")
    abline(v=exp(mean(log(r))+2*sd(log(r))),col="red")
    
    plot(r1, k1/1000, xlim = c(0.2, 1.6), ylim = c(150, 700), xlab="r", ylab="k (1000 t)", col = "gray62")
    
    hist(k/1000, freq=F, xlim=c(0, 1.5 * max(k/1000)), xlab="k (1000 t)", main = "", col = "gray62")
    abline(v=exp(mean(log(k))),col="red", lwd=2)	
    abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
    abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")
    
    plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)" , col = "grey62")
    abline(v=mean(log(r)))
    abline(h=mean(log(k)))
    abline(mean(log(msy))+log(4),-1, col="red",lwd=2)
    abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
    abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
    
    hist(msy/1000, freq=F, xlim=c(0, 1.2 * max(msy/1000)), xlab="MSY (1000 t)",main = "", col = "gray62")
    abline(v=exp(mean(log(msy))),col="red", lwd=2)
    abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
    abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
    
    cat("Possible combinations = ", length(r),"\n")
    cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
    cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
    cat("geom. mean k =", format(exp(mean(log(k))),digits=3), "\n")
    cat("k +/- 2 SD =", format(exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
    cat("geom. mean MSY =", format(exp(mean(log(msy))),digits=3),"\n")
    cat("MSY +/- 2 SD =", format(exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
    cat("geom. mean BMSY =", format(exp(mean(log(bmsy))),digits=3),"\n")
    cat("BMSY +/- 2 SD =", format(exp(mean_ln_bmsy - 2 * sd(log(bmsy))),digits=3), "-", format(exp(mean_ln_bmsy + 2 * sd(log(bmsy))),digits=3), "\n")
    
    output = data.frame(stock, sigR, startbio[1], startbio[2], interbio[1], interbio[2], finalbio[1], finalbio[2], 
                        min(yr), max(yr), res, max(ct), ct[1], ct[nyr], length(r), exp(mean(log(r))), sd(log(r)), min(r), 
                        quantile(r,0.05), quantile(r,0.25), median(r), quantile(r,0.75), quantile(r,0.95), max(r), exp(mean(log(k))), 
                        sd(log(k)), min(k), quantile(k, 0.05), quantile(k, 0.25), median(k), quantile(k, 0.75), quantile(k, 0.95), 
                        max(k), exp(mean(log(msy))), sd(log(msy)), min(msy), quantile(msy, 0.05), quantile(msy, 0.25), median(msy), 
                        quantile(msy, 0.75), quantile(msy, 0.95), max(msy),
                        exp(mean(log(bmsy))), sd(log(bmsy)), min(bmsy), quantile(bmsy, 0.05), quantile(bmsy, 0.25), median(bmsy), 
                        quantile(bmsy, 0.75), quantile(bmsy, 0.95), max(bmsy)) 
    
    write.csv(output, file = outfile)
    
  }
} 
