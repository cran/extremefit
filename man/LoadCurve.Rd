\name{LoadCurve}
\alias{LoadCurve}
\docType{data}
\title{Load curve of an habitation}
\description{The data frame provides electric consumption of an habitation in France over one month.}
\usage{data("LoadCurve")}
\format{The data is the electric consumption of an habitation in Kilovolt-amps (kVA) every 10 minutes during one month. The habitation has a contract that allows a maximum power of 6 kVA.A list of 2 elements.
  \describe{
    \item{\code{$data : a data frame with 24126 observations for 2 variables}}{\describe{
\describe{
\item{\code{Time}}{the number of day since the 1st of January, 1970.}
\item{\code{Value}}{the value of the electric consumtion in kVA.}
}}}
\item{\code{$Tgrid : A grid of time to perform the procedure.}}{}
    } 
}

\source{Electricite Reseau Distribution France}
\examples{
data("LoadCurve")

X<-LoadCurve$data$Value
days<-LoadCurve$data$Time
new.Tgrid <- LoadCurve$Tgrid
\dontrun{ #For computing time purpose
# Choice of the bandwidth by cross validation. 
# We choose the Gaussian kernel and the critical value 2.7. 
# As the computing time is high, we give the value of the bandwidth.

#hgrid <- bandwidth.grid(0.8, 5, 60)
#hcv<-bandwidth.CV(X, days, new.Tgrid, hgrid, pcv = 0.99, 
# kernel = Gaussian.kernel, #CritVal = 2.7, plot = FALSE)
#h.cv <- hcv$h.cv

h.cv <- 3.444261
HH<-hill.ts(X, days, new.Tgrid, h=h.cv, kernel = Gaussian.kernel,CritVal = 2.7)

Tgrid<-seq(min(days),max(days),length=400)
Quant<-rep(NA,length(Tgrid))
Quant[match(new.Tgrid, Tgrid)]<-as.numeric(predict(HH, 
            X, days, pgrid = 0.99, type = "quantile")$quantiles)

plot(Date, X/1000, ylim = c(0, 8),
      type = "l", ylab = "Electric consumption (kVA)", xlab = "Time")

lines(as.POSIXlt((Tgrid)*86400, origin = "1970-01-01",
                 tz = "Europe/Paris"), Quant/1000, col = "red")
}
}
