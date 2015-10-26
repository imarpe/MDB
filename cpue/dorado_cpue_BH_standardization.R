# Clean:
rm(list = ls()) #clear all
graphics.off() # close all
gc() # Clear memmory 

source("cpue_standardization_functions.R")

# Data pre-procesamiento
perico = read.csv("perico.csv")
catch = (perico$pericokg)/1000 #expresada en toneladas
effort = perico$nHook   #elecci?n del esfuerzo

#Definiendo los grupos de embarcaciones por capacidad de bogeda
breaks = c(5,10,20)
perico$holdCapacityGroup = cut(perico$holdCapacity, breaks=c(0, breaks, Inf), labels = seq(length(breaks)+1))

#Estandarizaci?n CPUE (Beverton y Holt)
BHcpue(perico, 2, plot = TRUE, col = "red")
title("Artisanal 2000 - 2014")

