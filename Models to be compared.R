library(forecast)
library(xts)
library(caret)
library(WaveletANN)
library(dplyr)
library(WaveletArima)
library(zoo)
library(changepoint)
library(strucchange)

testingmodel.holt<-function(train, test){
  x<-train; n.forecast<-nrow(test)
  x.exp<-holt(x, h = n.forecast)
  x.resid<-residuals(x.exp)
  hyb<-x.exp$fitted
  hyb.ext<-x.exp$mean
  k1<-xyplot.ts(as.xts(cbind.data.frame("True Value" = x, "Fitted Value" = hyb)), 
                superpose = T)
  k2<-xyplot.ts(as.xts(cbind.data.frame("True Value" = test, "Fitted Value" = hyb.ext)), 
                superpose = T)
  mylist<-list(k1, k2, "Fitted.Values" = hyb, "Forecasted.Values" = hyb.ext, 
               "In sample RMSE" = RMSE(hyb, x), "In sample MAE" = MAE(hyb, x),
               "Out of sample RMSE" = RMSE(hyb.ext, test), 
               "Out of sample MAE" = MAE(hyb.ext, test))
  return(mylist)
}


testingmodel.holt.wbann<-function(train, test){
  x<-train; n.forecast<-nrow(test)
  x.exp<-holt(x, h = n.forecast)
  x.resid<-residuals(x.exp)
  wv<-WaveletFittingann(x.resid, Waveletlevels = floor(log(length(x.resid))), 
                        boundary = "periodic", FastFlag = T, nonseaslag = 5, seaslag = 5, 
                        NForecast = n.forecast)
  hyb<-wv$FinalPrediction+x.exp$fitted
  hyb.ext<-wv$Finalforecast+x.exp$mean
  k1<-xyplot.ts(as.xts(cbind.data.frame("True Value" = x, "Fitted Value" = hyb)), 
                superpose = T)
  k2<-xyplot.ts(as.xts(cbind.data.frame("True Value" = test, "Fitted Value" = hyb.ext)), 
                superpose = T)
  mylist<-list(k1, k2, "Fitted.Values" = hyb, "Forecasted.Values" = hyb.ext, 
               "In sample RMSE" = RMSE(hyb, x), "In sample MAE" = MAE(hyb, x),
               "Out of sample RMSE" = RMSE(hyb.ext, test), 
               "Out of sample MAE" = MAE(hyb.ext, test))
  return(mylist)
}



testingmodel.arima.wbf<-function(train, test){
  x<-train; n.forecast<-nrow(test)
  x.arima<-auto.arima(x)
  arima.forecast<-forecast(x.arima, h = n.forecast)
  x.resid<-residuals(x.arima)
  wv<-WaveletFittingarma(x.resid, Waveletlevels = floor(log(length(x.resid))), 
                         boundary = "periodic", FastFlag = T, MaxARParam = 5, 
                         MaxMAParam = 5, NForecast = n.forecast)
  hyb<-wv$FinalPrediction+x.arima$fitted
  hyb.ext<-wv$Finalforecast+arima.forecast$mean
  k1<-xyplot.ts(as.xts(cbind.data.frame("True Value" = x, "Fitted Value" = hyb)), 
                superpose = T)
  k2<-xyplot.ts(as.xts(cbind.data.frame("True Value" = test, "Fitted Value" = hyb.ext)), 
                superpose = T)
  mylist<-list(k1, k2, "ARIMA Order (p, d, q)" = x.arima$arma[c(1,6,2)], "Fitted.Values" = hyb, 
               "Forecasted.Values" = hyb.ext, 
               "In sample RMSE" = RMSE(hyb, x), "In sample MAE" = MAE(hyb, x),
               "Out of sample RMSE" = RMSE(hyb.ext, test), 
               "Out of sample MAE" = MAE(hyb.ext, test))
  return(mylist)
}

testingmodel.arima<-function(train, test){
  x<-train; n.forecast<-nrow(test)
  x.arima<-auto.arima(x)
  arima.forecast<-forecast(x.arima, h = n.forecast)
  k1<-xyplot.ts(as.xts(cbind.data.frame("True Value" = x, "Fitted Value" = x.arima$fitted)), 
                superpose = T)
  k2<-xyplot.ts(as.xts(cbind.data.frame("True Value" = test, "Fitted Value" = arima.forecast$mean)), 
                superpose = T)
  mylist<-list(k1, k2, "ARIMA Order (p, d, q)" = x.arima$arma[c(1,6,2)], 
               "Fitted.Values" = x.arima$fitted, "Forecasted.Values" = arima.forecast$mean, 
               "In sample RMSE" = RMSE(x.arima$fitted, x), 
               "In sample MAE" = MAE(x.arima$fitted, x),
               "Out of sample RMSE" = RMSE(arima.forecast$mean, test), 
               "Out of sample MAE" = MAE(arima.forecast$mean, test))
  
  return(mylist)
}

ensemble<-function(data, percent.train){
  k<-as.integer((percent.train/100)*nrow(data))
  train<-data[1:k,]
  test<-data[(k+1):nrow(data),]
  arima1<-testingmodel.arima(train, test)
  arima.wbf<-testingmodel.arima.wbf(train, test)
  holt1<-testingmodel.holt(train, test)
  holt.ann<-testingmodel.holt.wbann(train, test)
  c1<-c("ARIMA", "ARIMA w/ WBF", "Holt", "Holt's Model w/ WBANN")
  c2<-c(arima1$`In sample RMSE`, arima.wbf$`In sample RMSE`, holt1$`In sample RMSE`, holt.ann$`In sample RMSE`)
  c3<-c(arima1$`In sample MAE`, arima.wbf$`In sample MAE`, holt1$`In sample MAE`, holt.ann$`In sample MAE`)
  c4<-c(arima1$`Out of sample RMSE`, arima.wbf$`Out of sample RMSE`, holt1$`Out of sample RMSE`, holt.ann$`Out of sample RMSE`)
  c5<-c(arima1$`Out of sample MAE`, arima.wbf$`Out of sample MAE`, holt1$`Out of sample MAE`, holt.ann$`Out of sample MAE`)
  tab<-cbind.data.frame("Model Name" = c1, "In sample RMSE" = c2, "In sample MAE" = c3,
                        "Out of sample RMSE" = c4, "Out of sample MAE" = c5)
  return(tab)
}


model.holt.wbann<-function(data, n.forecast){
  x<-data
  x.exp<-holt(x, h = n.forecast)
  x.resid<-residuals(x.exp)
  wv<-WaveletFittingann(x.resid, Waveletlevels = floor(log(length(x.resid))), 
                        boundary = "periodic", FastFlag = T, nonseaslag = 5, seaslag = 5, 
                        NForecast = n.forecast)
  hyb<-wv$FinalPrediction+x.exp$fitted
  hyb.ext<-wv$Finalforecast+x.exp$mean
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(x)[length(x)]+i, after = length(dates))
  hyb.ext<-as.xts(cbind.data.frame("Forecasts" = hyb.ext), order.by = dates)
  
  op1<-vector("integer", length = n.forecast)
  op1[1]<-sum(hyb.ext[1],x)
  if(n.forecast > 1){
  for (i in 2:nrow(hyb.ext))
    op1[i]<-op1[i-1]+hyb.ext[i]
  }
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(data)[length(data)]+i, after = length(dates))
  ext<-as.xts(cbind.data.frame(hyb.ext, op1), order.by = dates)
  colnames(ext)<-c("Daily", "Cumulative")

  return(ext)
  
}

model.holt<-function(data, n.forecast){
  x<-data
  x.exp<-holt(x, h = n.forecast)
  x.resid<-residuals(x.exp)
  hyb<-x.exp$fitted
  hyb.ext<-x.exp$mean
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(x)[length(x)]+i, after = length(dates))
  hyb.ext<-as.xts(cbind.data.frame("Forecasts" = hyb.ext), order.by = dates)
  
  op1<-vector("integer", length = n.forecast)
  op1[1]<-sum(hyb.ext[1],x)
  for (i in 2:nrow(hyb.ext))
    op1[i]<-op1[i-1]+hyb.ext[i]
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(data)[length(data)]+i, after = length(dates))
  ext<-as.xts(cbind.data.frame(hyb.ext, op1), order.by = dates)
  colnames(ext)<-c("Daily", "Cumulative")
  
  return(ext)
  
}


model.arima.wbf<-function(data, n.forecast){
  x<-data
  x.arima<-auto.arima(x)
  arima.forecast<-forecast(x.arima, h = n.forecast)
  x.resid<-residuals(x.arima)
  wv<-WaveletFittingarma(x.resid, Waveletlevels = floor(log(length(x.resid))), 
                         boundary = "periodic", FastFlag = T, MaxARParam = 5, 
                         MaxMAParam = 5, NForecast = n.forecast)
  hyb<-wv$FinalPrediction+x.arima$fitted
  hyb.ext<-wv$Finalforecast+arima.forecast$mean
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(x)[length(x)]+i, after = length(dates))
  hyb.ext<-as.xts(cbind.data.frame("Forecasts" = hyb.ext), order.by = dates)
  op1<-vector("integer", length = n.forecast)
  op1[1]<-sum(hyb.ext[1],x)
  for (i in 2:nrow(hyb.ext))
    op1[i]<-op1[i-1]+hyb.ext[i]
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(data)[length(data)]+i, after = length(dates))
  ext<-as.xts(cbind.data.frame(hyb.ext, op1), order.by = dates)
  colnames(ext)<-c("Daily", "Cumulative")
  
  return(ext)
}

model.arima<-function(data, n.forecast){
  x<-data
  x.arima<-auto.arima(x)
  arima.forecast<-forecast(x.arima, h = n.forecast)
  hyb.ext<-arima.forecast$mean
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(x)[length(x)]+i, after = length(dates))
  hyb.ext<-as.xts(cbind.data.frame("Forecasts" = hyb.ext), order.by = dates)
  op1<-vector("integer", length = n.forecast)
  op1[1]<-sum(hyb.ext[1],x)
  for (i in 2:nrow(hyb.ext))
    op1[i]<-op1[i-1]+hyb.ext[i]
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(data)[length(data)]+i, after = length(dates))
  ext<-as.xts(cbind.data.frame(hyb.ext, op1), order.by = dates)
  colnames(ext)<-c("Daily", "Cumulative")
  
  return(ext)
}

ensemble.final<-function(data, n.forecast, tuner.percent){
  k<-ensemble(data, tuner.percent)
  acc<-(k$`Out of sample RMSE`+k$`Out of sample MAE`)/2
  model<-k[(acc == min(acc)),1]
  arima1<-model.arima(data, n.forecast)
  arimawbf1<-model.arima.wbf(data, n.forecast)
  holt1<-model.holt(data, n.forecast)
  holtwbann1<-model.holt.wbann(data, n.forecast)
  if(model==levels(k$`Model Name`)[1])
    {op<-arima1; mdl<-"ARIMA"}
  if(model==levels(k$`Model Name`)[2])
    {op<-arimawbf1; mdl<-"ARIMA W/ WBF"}
  if(model==levels(k$`Model Name`)[3])
    {op<-holt1; mdl<-"HOLT"}
  if(model==levels(k$`Model Name`)[4])
  {op<-holtwbann1; mdl<-"HOLT W/ WBANN"}
  
  op1<-vector("numeric", length = n.forecast)
  op1[1]<-sum(op[1],data)
  for (i in 2:nrow(op))
    op1[i]<-op1[i-1]+op[i]
  dates<-NULL
  for (i in 1:n.forecast)
    dates<-append(dates, index(data)[length(data)]+i, after = length(dates))
  ext<-as.xts(cbind.data.frame(op, op1), order.by = dates)
  colnames(ext)<-c("Daily", "Cumulative")
  ext1<-list("Forecasts" = ext, "Model" = mdl)
  return(ext1)
}


adj.mod.1<-function(States, Country){
  st.fore<-NULL
  for (i in 1:ncol(States)){
    k<-model.holt.wbann(States[,i], 1)
    st.fore<-append(st.fore, as.numeric(k[,1]), after = length(st.fore))
  }
  st.fore<-append(st.fore, sum(st.fore), after = length(st.fore))
  k1<-as.character(colnames(States))
  k1<-append(k1, "Ind", after = length(k1))
  zz<-cbind.data.frame("States" = k1, "Forecast" = st.fore)
  return(zz)
}

adj.mod.2<-function(States, Country, w){
  st.fore<-NULL
  for (i in 1:ncol(States)){
    k<-model.holt.wbann(States[,i], 1)
    st.fore<-append(st.fore, as.numeric(k[,1]), after = length(st.fore))
  }
  kk<-model.holt.wbann(Country,1)
  d1<-as.numeric(kk[,1])-sum(st.fore)
  d<-rep(d, times = length(w))
  st.fore.1<-st.fore+w*d
  st.fore.1<-append(st.fore.1, as.numeric(kk[,1]), after = length(st.fore.1))
  k1<-as.character(colnames(States))
  k1<-append(k1, "Ind", after = length(k1))
  zz<-cbind.data.frame("States" = k1, "Forecast" = st.fore.1)
  return(zz)
}




adj.holt.wbann<-function(States, Country){
  z<-NULL
  z.s<-NULL
  for (i in 1:ncol(States)) {
    forecast1<-model.holt.wbann(States[1:(nrow(States)-1),i],1)
    z<-append(z, as.numeric(forecast1[,1]-States[nrow(States),i]), after = length(z))
    z.s<-append(z, as.numeric(forecast1[,1]), after = length(z))
  }
  
  z1<-model.holt.wbann(Country[1:(nrow(Country)-1),1], 1)
  z.country<-as.numeric(z1[,1]-Country[nrow(Country)])
  
  adj.States<-NULL
  if(abs(sum(z)) <= abs(z.country)){
    adj.States<-adj.mod.1(States, Country)
   }
  else{
    w<-z^2/sum(z^2)
    adj.States<-adj.mod.2(States, Country, w)
  }
  return(adj.States)
}

#adj.out<-adj.holt.wbann(States, Country)
#adj.out[,2]<-round(adj.out[,2],0)
#adj.out

#z12<-NULL
#for (i in 1:ncol(States)) {
#  kk<-model.holt.wbann(States[,i],1)
#  z12<-append(z12, as.numeric(kk[,1]), after = length(z12))
#}


#kk2<-model.holt.wbann(Country, 1)
#z12<-append(z12, kk2[,1], after = ncol(States))
#View(as.integer(z12))

m<-150
dff<-Country
kk<-model.holt.wbann(dff[1:(nrow(dff)-m)], m)
k1<-Country[(nrow(dff)-m):(nrow(dff))]

ape<-abs((kk$Daily-k1)/k1)*100

plot(ape, col = "red", ylab = "Absolute Percent Error", 
     main = "Plot of APE vs Date", xlab = "Date", lwd = 2.5)

ape1<-cbind.data.frame(ape, "t" = c(1:nrow(ape)))
lm1<-lm(Daily~t, data = ape1)
plot(ape1$Daily, col = "red", ylab = "Absolute Percent Error",
     main = "Plot of APE w/ Regression Line", xlab = "Size of Testing Data (in Days)", lwd = 2.5)
abline(lm1, col = "darkgreen", lwd = 3)
ape_5perc <- as.numeric((5-lm1$coefficients[1])/lm1$coefficients[2])

print(round(ape_5perc, 0))


