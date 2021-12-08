ensemble.mw<-function(data, train.size, window.size){
  k<-train.size
  train<-data[1:k,]
  test<-data[(k+1):(k+window.size),]
  arima<-testingmodel.arima(train, test)
  arima.wbf<-testingmodel.arima.wbf(train, test)
  holt<-testingmodel.holt(train, test)
  holt.ann<-testingmodel.holt.wbann(train, test)
  c1<-c("ARIMA", "ARIMA w/ WBF", "Holt's Model", "Holt's Model w/ WBANN")
  c2<-c(arima$`In sample RMSE`, arima.wbf$`In sample RMSE`, holt$`In sample RMSE`, holt.ann$`In sample RMSE`)
  c3<-c(arima$`In sample MAE`, arima.wbf$`In sample MAE`, holt$`In sample MAE`, holt.ann$`In sample MAE`)
  c4<-c(arima$`Out of sample RMSE`, arima.wbf$`Out of sample RMSE`, holt$`Out of sample RMSE`, holt.ann$`Out of sample RMSE`)
  c5<-c(arima$`Out of sample MAE`, arima.wbf$`Out of sample MAE`, holt$`Out of sample MAE`, holt.ann$`Out of sample MAE`)
  tab<-cbind.data.frame("Model Name" = c1, "In sample RMSE" = c2, "In sample MAE" = c3,
                        "Out of sample RMSE" = c4, "Out of sample MAE" = c5)
  return(tab)
}


model.perform<-function(state, window.size){

  arima.out<-NULL
  arima.wbf.out<-NULL
  holt.out<-NULL
  holt.wbann.out<-NULL
  init<-nrow(state)/2
  term<-nrow(state)-window.size

  for (i in init:term) {
    df<-ensemble.mw(state, i, window.size)
    out.temp<-(df$`Out of sample RMSE`+ df$`Out of sample MAE`)/2
    arima.out<-append(arima.out, out.temp[1], after = length(arima.out))
    arima.wbf.out<-append(arima.wbf.out, out.temp[2], after = length(arima.wbf.out))
    holt.out<-append(holt.out, out.temp[3], after = length(holt.out))
    holt.wbann.out<-append(holt.wbann.out, out.temp[4], after = length(holt.wbann.out))
  }

  k1<-cbind.data.frame("ARIMA" = arima.out, "ARIMA-WBF" = arima.wbf.out, "Holt" = holt.out, "Holt-WBANN" = holt.wbann.out)
  k2<-vector("numeric", length = nrow(k1))
  k3<-vector("numeric", length = nrow(k1))
  k4<-vector("numeric", length = nrow(k1))
  k5<-vector("numeric", length = nrow(k1))
  
  for(i in 1:nrow(k1)){
    k2[i]<-ifelse(k1[i,1]==min(k1[i,]), 1, 0)
    k3[i]<-ifelse(k1[i,2]==min(k1[i,]), 1, 0)
    k4[i]<-ifelse(k1[i,3]==min(k1[i,]), 1, 0)
    k5[i]<-ifelse(k1[i,4]==min(k1[i,]), 1, 0)
  }
  
  h1<-100*mean(k2)
  h2<-100*mean(k3)
  h3<-100*mean(k4)
  h4<-100*mean(k5)
  ll<-c("ARIMA", "ARIMA-WBF", "Holt", "Holt-WBANN")
  Timeline<-cbind.data.frame("Model" = ll, 
                             "Timeline %" = c(h1, h2, h3, h4))
  dates<-NULL
  for (i in 1:(as.integer(term-init+1)))
    dates<-append(dates, index(state)[init]+i, after = length(dates))
  
  Out.perform<-as.xts(as.data.frame(k1), order.by = dates)
  k<-xyplot.ts(Out.perform, superpose = T, main = "Model Performace over Moving Window", ylab = "Model Performance Metric")
  return(list(k, "Accuracy Timeline" = Timeline, "Performance" = Out.perform))
}

k<-model.perform(hotspots[,1], 4)
