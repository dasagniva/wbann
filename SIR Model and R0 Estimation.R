require("deSolve")

derivative_calc_func=function(t, x, vparameters){
  S = x[1]  
  I = x[2]  
  R = x[3]  
  
  with(as.list(vparameters),{
    
    npop = S+I+R   
    
    dS = -beta*S*I/npop            
    dI = +beta*S*I/npop - gamma*I  
    dR = +gamma*I                  
    
    vout = c(dS,dI,dR)
    list(vout)
  })
}

library("R0")
R0_and_mean_shape_estimation<-function(incidence
                                       ,dates
                                       ,N=100
                                       ,mean_min=0.1
                                       ,mean_max=10
                                       ,shape_min=0.1
                                       ,shape_max=10
                                       ,npop = 100000
                                       ,I_0 = 1
                                       ,R_0 = 0
                                       ,gamma = 1/14 ## 14 days to recover        
                                       
                                       
){
  mean.star<-seq(mean_min,mean_max,length.out = N)
  shape.star<-seq(shape_min,shape_max,length.out = N)
  R0.star<-matrix(NA,nrow=N^2,ncol=6)
  colnames(R0.star)<-c("mean","shape","R0","R0.lower","R0.upper","MSE")
  S_0 = npop-I_0-R_0
  
  tbegin = 0
  tend   = length(incidence)
  vt = seq(tbegin,tend,1)  
  
  k<-1
  for(i in 1:N){
    for(j in 1:N){
      tryCatch(
        {
          R0.star[k,"mean"]<-mean.star[i]
          R0.star[k,"shape"]<-shape.star[j]
          mGT<-generation.time("gamma"
                               ,c(mean.star[i]
                                  ,mean.star[j]))
          EG.R<-est.R0.EG(epid = incidence,t = dates ,GT = mGT)
          R0.star[k,"R0"]<-EG.R$R
          R0.star[k,"R0.lower"]<-EG.R$conf.int[1]
          R0.star[k,"R0.upper"]<-EG.R$conf.int[2]
          R0    = EG.R$R       
          beta  = R0*gamma  ## This is beta parameter of SIR model 
          
          vparameters = c(gamma=gamma,beta=beta)
          inits = c(S=S_0,I=I_0,R=R_0)
          
          solved_model = as.data.frame(lsoda(inits, vt, derivative_calc_func, vparameters))
          
          vI = solved_model$I
          R0.star[k,"MSE"]<-mean((vI-incidence)^2)
        }, error=function(e){}
      )
      rm(list = c("EG.R")) 
      if(k%%1000==0)cat("k = ",k,"\n")
      k<-k+1
    }
  }
  min.mse.indx<-which.min(R0.star[,"MSE"])
  result<-R0.star[min.mse.indx,]
  return(result)
  
}


x<-R0_and_mean_shape_estimation(incidence = as.vector(mh.conf), 
                             dates = as.Date(index(mh.conf)),
                             N=100, 
                             mean_min=0.1, 
                             mean_max=10, 
                             shape_min=0.1, 
                             shape_max=10, 
                             npop = 100000,
                             I_0 = 1, 
                             R_0 = 0, 
                             gamma = 1/14)
x
R0.mh<-x


head(hotspots)

R_naught<-list(NULL)
for (i in 1:(ncol(hotspots))) {
  x<-R0_and_mean_shape_estimation(incidence = as.vector(hotspots[,i]), 
                                  dates = as.Date(index(hotspots)),
                                  N=100, 
                                  mean_min=0.1, 
                                  mean_max=10, 
                                  shape_min=0.1, 
                                  shape_max=10, 
                                  npop = 100000,
                                  I_0 = 1, 
                                  R_0 = 0, 
                                  gamma = 1/14)
  R_naught[[i]]<-as.data.frame(x)
}
R_naught
final_R0<-as.data.frame(R_naught)
colnames(final_R0)<-colnames(hotspots)
View(final_R0)
write.csv(final_R0, "C:\\Users\\HP\\OneDrive\\Desktop\\Estimated_R0.csv")

cc<-R0_and_mean_shape_estimation(incidence = as.vector(Country), 
                             dates = as.Date(index(hotspots)),
                             N=100, 
                             mean_min=0.1, 
                             mean_max=10, 
                             shape_min=0.1, 
                             shape_max=10, 
                             npop = 100000,
                             I_0 = 1, 
                             R_0 = 0, 
                             gamma = 1/14)
