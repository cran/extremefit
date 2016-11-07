\name{dataWind}
\alias{dataWind}
\docType{data}
\title{Wind speed for Brest (France)}
\description{The data frame provides the wind speed of Brest from 1976 to 2005.}
\usage{data("dataWind")}
\format{The data is the wind speed in meters per second (m/s) every day from 1976 to 2005.
\describe{
\item{\code{Year}}{The year of the measure.}
\item{\code{Month}}{The month of the measure.}
\item{\code{Day}}{the day of the measure.}
\item{\code{Speed}}{The wind speed in meters per second}
}}

\examples{
data("dataWind")




\dontrun{


Xwind<-dataWind$Speed
# As the speed is in meter/sec, we transform the speed into a more common scale km/h

Xwind<-Xwind*3600/1000

n<-length(Xwind)
twind<-1:n

plot(twind,Xwind)



#Choice of a period : 1990 to 2000 for example

indices<-which((dataWind$Year>=1990)&(dataWind$Year<=2000))
X<-Xwind[indices]
t<-twind[indices]

#we choose a bandwidth equal to 60 (two months of data)

h<-60
seqTgrid<-seq(1,length(t),by=10)
Tgrid<-t[seqTgrid]
TS.W<-hill.ts(X,t=t,Tgrid,h=h,Gaussian.kernel,
               CritVal=2.7,gridlen=100,initprop = 1/4, r1 = 1/4, r2 = 1/20)

pgrid<-c(0.95,0.99) # can be changed to a vector of
# probabilities without affecting the computing time
pred.quant.W<-predict(TS.W,X,t=t,pgrid=pgrid,type="quantile")


plot(t,X,type="l",ylim=c(0,100))

quant<-as.numeric(pred.quant.W$quantile[1,])
lines(Tgrid,quant,col="red")

quant<-as.numeric(pred.quant.W$quantile[2,])
lines(Tgrid,quant,col="blue")


legend("bottomright",c("0.95-quantile","0.99-quantile"),col=c("red","blue"),pch="-")

# with the time on the axis

# date<-paste(dataWind$Year,"-",dataWind$Month,"-",dataWind$Day)
# dateWind<-as.Date(date,format="%Y - %m - %d")


# plot(dateWind,Xwind)

# Date<-dateWind[indices]

# dateGRid<-Date[seqTgrid]

# plot(Date,X,type="l",ylim=c(0,100))

# quant<-as.numeric(pred.quant.W$quantile[1,])
# lines(dateGRid,quant,col="red")

# quant<-as.numeric(pred.quant.W$quantile[2,])
# lines(dateGRid,quant,col="blue")


# legend("bottomright",c("0.95-quantile","0.99-quantile"),col=c("red","blue"),pch="-")
}

}
