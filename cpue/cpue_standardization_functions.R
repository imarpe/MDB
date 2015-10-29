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
    axis(side = 1, at = 1:length(unique(data$year)), las = 2,
         labels = seq(sort(unique(data$year))[1], rev(sort(unique(data$year)))[1]))
    axis(2, las=1)
    box()
  }
  return(stdCPUE/mean(stdCPUE, na.rm=TRUE))
}

BHcpue2 = function(data, stdVessel, plot = FALSE, col, ...){
  catchByHoldcapacity  = tapply(data$catch, INDEX = data$holdCapacityGroup, sum)
  effortByHoldcapacity = tapply(data$effort, INDEX = data$holdCapacityGroup, sum)
  std1                 = catchByHoldcapacity/effortByHoldcapacity
  std                  = data.frame(holdCapacityGroup = seq_along(unique(data$holdCapacityGroup)),
                                    RFP = as.vector(std1/std1[stdVessel]))
  data                 = merge(data, std)
  data$stdEffort       = data$effort*data$RFP
  CPUE                 = tapply(data$catch, INDEX = data$year, sum)/tapply(data$stdEffort, INDEX = data$year, sum)
  stdCPUE              = CPUE/mean(CPUE,na.rm = TRUE)
  if(plot==TRUE){
    plot(x = stdCPUE, axes=FALSE, type = "b", pch=19, lwd = 2, ylim = c(0.9*min(stdCPUE), 1.1*max(stdCPUE)),
         ylab = "Estandardized catch rate",xlab = "year", col = col)
    axis(side = 1, at = 1:length(unique(data$year)), las = 2,
         labels = seq(sort(unique(data$year))[1], rev(sort(unique(data$year)))[1]))
    axis(2, las = 1)
    box()
  }
  return(stdCPUE)
}


lag = function(x, lag=0, freq=12) {
  lag = lag %% freq
  if(lag==0) return(x)
  yini = x[1] - freq/12
  output = c(rep(yini, lag), x)
  length(output) = length(x)
  return(output)
}

daysTrip = function(date_ini, date_end, format = "%d/%m/%Y"){
  date_ini = as.character(date_ini)
  date_end = as.character(date_end)
  output     = as.numeric(as.Date(date_end, format = format) - as.Date(date_ini, format = format))
  if(length(date_ini) != length(date_end)){
    stop("date_ini and date_end do not have teh same length")
  }
  for(i in seq_along(output))
    if(output[i] < 0){
      stop("date_end must be greater than date_ini")
    }
  return(output)
}


