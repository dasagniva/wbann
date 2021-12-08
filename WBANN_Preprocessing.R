library(dplyr)
library(xts)
library(httr)
library(jsonlite)
res<-GET("https://api.covid19india.org/states_daily.json")
dataa<-fromJSON(rawToChar(res$content))
df<-dataa$states_daily
df<-df[,-9]
colnames(df)<-toupper(colnames(df))
colnames(df)[8]<-"Date"
colnames(df)[33]<-"Status"
k<-c(1:7, 9:32, 34:41)
for (i in k) 
  df[,i]<-as.numeric(df[,i])
df$Date<-as.Date(as.character(df$Date), format = "%d-%B-%y")
#View(df.conf)
df.conf<-df[df$Status=="Confirmed",]
df.conf<-na.omit(df.conf)
df.dec<-df[df$Status=="Deceased",]
df.rec<-df[df$Status=="Recovered",]
ind.conf<-as.xts(df.conf[,-c(8,9,33)], order.by = df.conf$Date)

#nrow(ind.conf)
Country<-ind.conf[1:272, 34]
States<-ind.conf[1:272,-34]

Country[174]
k<-c(20,2,16,32,7,17)
hotspots<-States[,k]


#Hotspots
mh.conf<-hotspots[,1] 
ap.conf<-hotspots[,2]
ka.conf<-hotspots[,3]
tn.conf<-hotspots[,4]
ct.conf<-hotspots[,5]
kl.conf<-hotspots[,6]

#Others
ot.conf<-States[,-k]
