\name{dataOyster}
\alias{dataOyster}
\docType{data}
\title{
High-frequency noninvasive valvometry data
}
\description{
 The data frame provides the opening amplitude of one oyster's shells (in mm) with respect to the time (in hours). The opening velocity of the oyster's shells is also given.
}
\usage{data("dataOyster")}
\format{
  A data frame with 18000 observations for 3 variables.
  \describe{
    \item{\code{time}}{Time of measurement (in hours).}
    \item{\code{opening}}{opening amplitude between the two shells (in mm).}
    \item{\code{velocity}}{a numeric vector (in mm/s).}
  }
}
\references{
Durrieu, G., Grama, I., Pham, Q. & Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.

Azais, R., Coudret R. & Durrieu G. (2014). A hidden renewal model for monitoring aquatic systems biosensors. Environmetrics, 25.3, 189-199.

Schmitt, F. G., De Rosa, M., Durrieu, G., Sow, M., Ciret, P., Tran, D., & Massabuau, J. C. (2011). Statistical study of bivalve high frequency microclosing behavior: Scaling properties and shot noise analysis. International Journal of Bifurcation and Chaos, 21(12), 3565-3576.

Sow, M., Durrieu, G., Briollais, L., Ciret, P., & Massabuau, J. C. (2011). Water quality assessment by means of HFNI valvometry and high-frequency data modeling. Environmental monitoring and assessment, 182(1-4), 155-170.

website : http://molluscan-eye.epoc.u-bordeaux1.fr/
}
\examples{

data("dataOyster")
Velocity<-dataOyster$data[,3]
time<-dataOyster$data[,1]
plot(time,Velocity,type="l",xlab="time (hour)",ylab="Velocity (mm/s)")

Tgrid<-seq(0,24,0.1) 
#Grid with positive velocity
new.Tgrid<-dataOyster$Tgrid

X<-Velocity+(-min(Velocity)) #We adjust the data to be positive

\dontrun{ #For computing time purpose
#We find the h by minimizing the cross validation function 
hgrid<-bandwidth.grid(0.05,0.8,30,type="geometric")
indices<-round(seq(1,length(new.Tgrid),length=30))
TgridCV<-new.Tgrid[indices]
#the following function takes 15.87 secs on 2,5 GHz Intel Core i7, 16 GB 1600 MHz DDR3
H<-bandwidth.CV(X,time,TgridCV,hgrid,
                  TruncGauss.kernel,kpar=c(sigma=1),pcv=0.99,CritVal=3.6,plot=TRUE)
hcv<-H$h.cv

#we use our method with the h found previously
TS.Oyster<-hill.ts(X,t=time,new.Tgrid,h=hcv,TruncGauss.kernel,kpar=c(sigma=1),
          CritVal=3.6,gridlen=100,initprop = 1/10, r1 = 1/4, r2 = 1/20)
          
plot(time,Velocity,type="l",ylim=c(-0.6,1.6),
    main="Extreme quantiles estimator",
    xlab="Time (hour)",ylab="Velocity (mm/s)")
pgrid<-c(0.999)
pred.quant.Oyster<-predict(TS.Oyster,X,t=time,pgrid=pgrid,type="quantile")

quant0.999<-rep(0,length(Tgrid))
for(ii in 1:length(new.Tgrid)){
  quant0.999[which(Tgrid==new.Tgrid[ii])]<-as.numeric(pred.quant.Oyster$quantile)[ii]-
  (-min(Velocity))
  }
lines(Tgrid,quant0.999,col="magenta")   

legend("topright",c("0.999-quantile"),col=c("magenta"),pch="-")

}
}
