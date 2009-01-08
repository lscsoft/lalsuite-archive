
#
# This is a simulation of fourth moment statistics produced in linear version of statistics functions
#

sum4<-function(nbins=501, half_window=10) {
	Y<-rnorm(nbins)
	max_dx_bin<-which.max(Y)
	Idx<-1:nbins
	F<- abs(Idx-max_dx_bin)>half_window
	Idx<-Idx[F]

	Y<-Y[Idx]

	M<-mean(Y)
	sum4<- mean((Y-M)^4)/sd(Y)^4
	return(sum4)
	}

Run<- unlist(lapply(1:1000, function(x)sum4()))
