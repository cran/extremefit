\name{dataLoadCurve}
\alias{dataLoadCurve}
\docType{data}
\title{Load curves from ERDF customers}
\description{The data frame provides electric consumption from ERDF customers in France.}
\usage{data("dataLoadCurve")}
\format{The data is the average electric consumption of the customers in Kilovolt-amps (kVA) every 30 minutes during 3 years. The averages are calculated according to the following four patterns:
\describe{
\item{\code{RES1}}{the average customers consumption with a 3 or 6 kVA contract.}
\item{\code{RES11}}{the average customers consumption with a 9 or more kVA contract.}
\item{\code{RES2HP}}{the average consumption of all customers during peak time.}
\item{\code{RES2HC}}{the average consumption of all customers during off-peak time.}
}}
\source{Electricite Reseau Distribution France}
\examples{
data("dataLoadCurve")
t<-dataLoadCurve$Time
X<-dataLoadCurve$RES1

# uncomment one of the following line to change the RES1 pattern to an other one.
#X<-dataLoadCurve$RES11
#X<-dataLoadCurve$RES2HP
#X<-dataLoadCurve$RES2HC

n<-length(X)

#choice of the bandwidth
\dontrun{ # For computation time purpose

# due to the periodicity in the data, we select the grid on one year.

L<-(max(t)-min(t))/3

ind<-which((t>=0)&(t<=L))
Xselect<-X[ind]
tselect<-t[ind]

Tgrid<-seq(min(tselect),max(tselect),length=30)

#choice of the bandwidth grid between hmin and hmax
hmin<-L/12/4 #One week
hmax<-L/4  #four weeks
hgrid<-bandwidth.grid(hmin,hmax,20,type="geometric")


#performs the cross validation function
H<-bandwidth.CV(Xselect,tselect,Tgrid,hgrid,
              Gaussian.kernel,pcv=0.99,CritVal=2.7,plot=TRUE)

hcv<-H$h.cv


#Calculate the extremes quantiles on one year

#choice of the evaluation grid: one measurement per day
# (this corresponds to the choice with
# step = 30 minutes * 2 * 24h
# since the data are collected every 30 minutes)

TMax<-(max(t)-min(t))/3
Tgrid<-seq(min(t),TMax,by=30*2*24)
ind<-which((t>=min(t))&(t<=TMax))
Xyear<-X[ind]
tyear<-t[ind]

TS.LC<-hill.ts(Xyear,t=tyear,Tgrid,h=hcv,Gaussian.kernel,
              CritVal=2.7,gridlen=100,initprop = 1/10, r1 = 1/4, r2 = 1/20)
pgrid<-c(0.95,0.99,0.999) # can be changed to a vector
# of probabilities without affecting the computing time
pred.quant.LC<-predict(TS.LC,X,t=t,pgrid=pgrid,type="quantile")


plot(tyear,Xyear,type="l",xlab="t (minutes)",col="lightgrey",
  ylab="Electric consumption (kVA)",main="Extreme quantiles",
  ylim=c(0.5,2.3))

quant<-as.numeric(pred.quant.LC$quantile[1,])
lines(Tgrid,quant,col="red")

quant<-as.numeric(pred.quant.LC$quantile[2,])
lines(Tgrid,quant,col="blue")

quant<-as.numeric(pred.quant.LC$quantile[3,])
lines(Tgrid,quant,col="magenta")


legend("bottomright",c("0.95-quantile","0.99-quantile","0.999-quantile"),
      col=c("red","blue","magenta"),pch="-")




#Calculate the extremes quantiles on three years

#choice of the evaluation grid: one measurement per day
# (this corresponds to the choice with step = 30 minutes * 2 * 24h
# since the data are collected every 30 minutes)

Tgrid<-seq(min(t),max(t),by=30*2*24)


TS.LC<-hill.ts(X,t=t,Tgrid,h=hcv,Gaussian.kernel,
                     CritVal=2.7,gridlen=100,initprop = 1/10, r1 = 1/4, r2 = 1/20)
pgrid<-c(0.95,0.99,0.999) # can be changed to a vector of
# probabilities without affecting the computing time
pred.quant.LC<-predict(TS.LC,X,t=t,pgrid=pgrid,type="quantile")


plot(t,X,type="l",xlab="t (minutes)",col="lightgrey",
    ylab="Electric consumption (kVA)",main="Extreme quantiles",
    ylim=c(0.5,2.3))

quant<-as.numeric(pred.quant.LC$quantile[1,])
lines(Tgrid,quant,col="red")

quant<-as.numeric(pred.quant.LC$quantile[2,])
lines(Tgrid,quant,col="blue")

quant<-as.numeric(pred.quant.LC$quantile[3,])
lines(Tgrid,quant,col="magenta")

legend("bottomright",c("0.95-quantile","0.99-quantile","0.999-quantile"),
      col=c("red","blue","magenta"),pch="-")
}

}
