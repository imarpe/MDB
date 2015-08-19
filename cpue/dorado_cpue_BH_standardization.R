# Clean:
rm(list = ls()) #clear all
graphics.off() # close all
gc() # Clear memmory 

source("cpue_standardization_functions.R")

# Data pre-procesamiento
dorado = read.csv("dorado_data.csv")
catch = (dorado$doradokg)/1000 #in tons
effort = dorado$nHook   #effort

#groups hold capacity
breaks = c(5, 10, 20)
dorado$holdCapacityGroup = cut(dorado$holdCapacity, breaks=c(0, breaks, Inf), labels = seq(length(breaks)+1))

#CPUE Standardization (Beverton y Holt)
BHcpue(dorado, 2, plot = TRUE, col = "red")
