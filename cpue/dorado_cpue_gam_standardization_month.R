# Clean:
rm(list = ls(all = TRUE)) #clear all;
graphics.off() # close all;
gc() # Clear memmory (residuals of operations?, cache? Not sure)

require(gam)
require(pgirmess)
source("pota_functions.R")

# Data pre-processing -----------------------------------------------------

perico = read.csv("perico.csv")
perico = perico[perico$year < 2015, ]
names(perico)
perico$nTrips = rep(1, length=nrow(perico))
perico$daysTrip = daysTrip(perico$date_ini, perico$date_end)
environment = read.table("http://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices",
                   header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
environment = environment[ ,c("YR","MON", "NINO1.2", "ANOM")]
colnames(environment) = c("year", "month", "nino12", "anom")

# enviroment = read.csv("enviroment.csv")
# mei = mei[, c("year", "mei")]
# mei$month = rep(1:12, length(nrow(mei)))


factor = 1e-3
names(perico)
pericoByMonth = aggregate(cbind(total, pericokg, holdCapacity, nHook, daysTrip, nHours, nTrips) ~ year + month, FUN = sum, 
                         data = perico, na.action=na.pass)
newBase = expand.grid(year = unique(pericoByMonth$year), month = c(1:12))
pericoByMonth = merge(newBase, pericoByMonth, all = TRUE)
pericoByMonth = merge(pericoByMonth, environment, all.x = TRUE)
pericoByMonth$time = pericoByMonth$year + ((pericoByMonth$month)/12)
pericoByMonth = pericoByMonth[order(pericoByMonth$time), ]
pericoByMonth = pericoByMonth[pericoByMonth$year >= 1999 & pericoByMonth$year <= 2015, ]
pericoByMonth$semester = ifelse(pericoByMonth$month<=6, pericoByMonth$year+0, pericoByMonth$year+0.5)
pericoByMonth$quarter  = pericoByMonth$year + rep(c(0, 0.25, 0.5, 0.75), each=3)

season = rep(c("summer", "fall", "winter", "spring"), each = 3, len = nrow(pericoByMonth) + 1)

# year effect
for(i in 0:11) pericoByMonth[, paste0("year", i)] = as.factor(lag(pericoByMonth$year, i))

# semester effect
for(i in 0:5) pericoByMonth[, paste0("semester", i)] = as.factor(lag(pericoByMonth$semester, i, freq=6))

# quarter effect
for(i in 0:2) pericoByMonth[, paste0("quarter", i)] = as.factor(lag(pericoByMonth$quarter, i, freq=3))

pericoByMonth$month = as.factor(pericoByMonth$month)
pericoByMonth$year = as.factor(pericoByMonth$year)
pericoByMonth$season1 = as.factor(season[-length(season)])
pericoByMonth$season2 = as.factor(season[-1])
pericoByMonth$yearSeason = pericoByMonth$year3

# CPUE:  ----------------------------------------------------

pericoByMonth$cpue = (factor*pericoByMonth$pericokg)/pericoByMonth$nHook
pairsrp(pericoByMonth[, c("cpue", "holdCapacity", "nHook", "nHours", "nTrips")], meth = "pearson",
        cex = 1.5, col = "grey50")

# year effect
tot.year = NULL
for(i in 0:11) {
  fmla = sprintf("log(cpue) ~  year%d", i)
  tot.year[[paste0("year", i)]] = gam(as.formula(fmla), data=pericoByMonth)
  print(c(round(AIC(tot.year[[i+1]]),2), names(tot.year)[i+1]))
}

  
# semester effect
tot.sem = NULL
for(i in 0:5) {
  fmla = sprintf("log(cpue) ~ semester%d", i)
  tot.sem[[paste0("semester", i)]] = gam(as.formula(fmla), data=pericoByMonth)
  print(c(round(AIC(tot.sem[[i+1]]),2), names(tot.sem)[i+1]))
}

# quarter effect
tot.quart = NULL
for(i in 0:2) {
  fmla = sprintf("log(cpue) ~  quarter%d + month", i)
  tot.quart[[paste0("quarter", i)]] = gam(as.formula(fmla), data=pericoByMonth)
  print(c(round(AIC(tot.quart[[i+1]]),2), names(tot.quart)[i+1]))
}


#Modelos
model = NULL 

model[[1]] = gam (log(cpue) ~ year10, data=pericoByMonth)
model[[2]] = gam (log(cpue) ~ semester4, data=pericoByMonth)
model[[3]] = gam (log(cpue) ~ quarter0, data=pericoByMonth)

model[[4]] = gam (log(cpue) ~ year10 + season2, data=pericoByMonth)
model[[5]] = gam (log(cpue) ~ semester4 + season2, data=pericoByMonth)
model[[6]] = gam (log(cpue) ~ quarter0 + season2, data=pericoByMonth)

model[[7]] = gam (log(cpue) ~ year10 + season2 + s(nino12), data=pericoByMonth)
model[[8]] = gam (log(cpue) ~ semester4 + season2 + s(nino12), data=pericoByMonth)
model[[9]] = gam (log(cpue) ~ quarter0 + season2 + s(nino12), data=pericoByMonth)

model[[10]] = gam (log(cpue) ~ year10 + season2 + s(anom), data=pericoByMonth)
model[[11]] = gam (log(cpue) ~ semester4 + season2 + s(anom) , data=pericoByMonth)
model[[12]] = gam (log(cpue) ~ quarter0 + season2 + s(anom) , data=pericoByMonth)

model[[13]] = gam (log(cpue) ~ year10 + season2 + s(nHours) , data=pericoByMonth)
model[[14]] = gam (log(cpue) ~ semester4 + season2 + s(nHours) , data=pericoByMonth)
model[[15]] = gam (log(cpue) ~ quarter0 + season2 + s(nHours) , data=pericoByMonth)

for(i in seq_along(model)){
  print(AIC(model[[i]]))
}

for(i in seq_along(model)){
  print(summary(model[[i]]))
}

pericoByMonth$cpueP1[!is.na(pericoByMonth$cpue)] = exp(predict(model[[1]]))
pericoByMonth$cpueP2[!is.na(pericoByMonth$cpue)] = exp(predict(model[[2]]))
pericoByMonth$cpueP3[!is.na(pericoByMonth$cpue)] = exp(predict(model[[3]]))
pericoByMonth$cpueP4[!is.na(pericoByMonth$cpue)] = exp(predict(model[[4]]))
pericoByMonth$cpueP5[!is.na(pericoByMonth$cpue)] = exp(predict(model[[5]]))
pericoByMonth$cpueP6[!is.na(pericoByMonth$cpue)] = exp(predict(model[[6]]))
pericoByMonth$cpueP7[!is.na(pericoByMonth$cpue)] = exp(predict(model[[7]]))
pericoByMonth$cpueP8[!is.na(pericoByMonth$cpue)] = exp(predict(model[[8]]))
pericoByMonth$cpueP9[!is.na(pericoByMonth$cpue)] = exp(predict(model[[9]]))
pericoByMonth$cpueP10[!is.na(pericoByMonth$cpue)] = exp(predict(model[[10]]))
pericoByMonth$cpueP11[!is.na(pericoByMonth$cpue)] = exp(predict(model[[11]]))
pericoByMonth$cpueP12[!is.na(pericoByMonth$cpue)] = exp(predict(model[[12]]))
pericoByMonth$cpueP13[!is.na(pericoByMonth$cpue)] = exp(predict(model[[13]]))
pericoByMonth$cpueP14[!is.na(pericoByMonth$cpue)] = exp(predict(model[[14]]))
pericoByMonth$cpueP15[!is.na(pericoByMonth$cpue)] = exp(predict(model[[15]]))

plot(pericoByMonth$time, pericoByMonth$cpue, type="b", col="gray", xlim=c(1999, 2015),
     ylab = "CPUE", xlab = "time")
lines(pericoByMonth$time, pericoByMonth$cpueP9, col="blue", lwd=2)
lines(pericoByMonth$time, pericoByMonth$cpueP15, col="red", lwd=2)

print(cor(pericoByMonth[, grep(pat="cpue", names(pericoByMonth),value=TRUE)], use="complete")[-1,1])



##############
ref = "2000"
MODELyear = 7 

#Asumiendo year
x = model[[MODELyear]]
xperico = complete.cases(pericoByMonth[, c("nHook")])
pT = predict(x, type="terms", se=TRUE)
pT[["index"]] = pericoByMonth$year10[xperico]
pT[["time"]] = pericoByMonth$year10[xperico]

split(pT$fit[, "year10"], f = pT$index)
lyear = tapply(pT$fit[, "year10"], INDEX=pT$index, FUN=unique)
se    = tapply(pT$se.fit[, "year10"], INDEX=pT$index, FUN=unique)
year = exp(lyear + 0.5*se^2)
year.perico = year/year[ref]

year.perico = data.frame(time = as.numeric(names(lyear)), 
                         lyear = lyear, se = se, year = year, ind = year.perico)

#par(mfrow=c(1,2), mar=c(3,3,1,1))
#par(mfrow=c(1,2))
plot(year.perico$time[-c(1:3)], year.perico$year[-c(1:3)], type = "b", lwd = 2, pch = 19, 
     xlab = "Year", ylab = "Standardized catch rate", axes = FALSE, ylim = c(0, 2))
axis(1, las = 2, seq(1998, 2014, 1))
axis(2, las = 1)
box()
title("year")


#Asumiendo quarter
MODELquarter = 15

y = model[[MODELquarter]]
yperico = complete.cases(pericoByMonth[, c("nHook")])
pT = predict(y, type="terms", se=TRUE)
pT[["index"]] = pericoByMonth$quarter0[yperico]
pT[["time"]] = pericoByMonth$quarter0[yperico]

split(pT$fit[, "quarter0"], f = pT$index)
lyear = tapply(pT$fit[, "quarter0"], INDEX=pT$index, FUN=unique)
se    = tapply(pT$se.fit[, "quarter0"], INDEX=pT$index, FUN=unique)
year = exp(lyear + 0.5*se^2)
quarter.perico = year/year[ref]

quarter.perico = data.frame(time=as.numeric(names(lyear)), 
                         lyear=lyear, se=se, year=year, ind=quarter.perico)

plot(quarter.perico$time, quarter.perico$year, type = "b", lwd = 2, pch = 19, 
     xlab = "Year", ylab = "Standardized catch rate", col = "red", axes = FALSE)
axis(1, las = 2, seq(1999, 2015, 1))
axis(2, las = 1)
box()
title("quarter")
