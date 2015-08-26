BHcpue = function(data, stdVessel, plot = FALSE, col, ... ){
  baseGroup            = aggregate(cbind(catch, effort) ~ year + holdCapacityGroup, FUN = sum, data = data, na.action=na.pass)
  catchByHoldcapacity  = tapply(baseGroup$catch, INDEX = baseGroup$holdCapacityGroup, sum)
  effortByHoldcapacity = tapply(baseGroup$effort, INDEX = baseGroup$holdCapacityGroup, sum)
  std1                 = catchByHoldcapacity/effortByHoldcapacity
  std                  = data.frame(holdCapacityGroup = seq_along(unique(baseGroup$holdCapacityGroup)), RFP = as.vector(std1/std1[stdVessel]))
  baseGroup            = merge(baseGroup, std)
  baseGroup$stdEffort  = baseGroup$effort*baseGroup$RFP
  stdCPUE              = tapply(baseGroup$catch, INDEX = baseGroup$year, sum)/tapply(baseGroup$stdEffort, INDEX = baseGroup$year, sum)
  
  if(plot==TRUE){
    plot(x = stdCPUE/mean(stdCPUE,na.rm = TRUE), axes=FALSE, type = "b", pch=19, lwd = 2,
         ylab = "Estandardized catch rate",xlab = "year", col = col)
    axis(side = 1, at = 1:length(unique(data$year)),
         labels = seq(sort(unique(data$year))[1], rev(sort(unique(data$year)))[1]))
    axis(2,las=1)
    box()
  }
  return(stdCPUE/mean(stdCPUE, na.rm=TRUE))
}


lag = function(x, lag=0, freq=12) {
  lag = lag %% freq
  if(lag==0) return(x)
  yini = x[1] - freq/12
  output = c(rep(yini, lag), x)
  length(output) = length(x)
  return(output)
}





