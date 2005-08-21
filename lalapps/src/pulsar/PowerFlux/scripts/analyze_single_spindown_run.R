source("params.R")


vn<-function(field, i) {
	return(get(paste(var_prefix, field, i, sep="")))
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}

#
# Load large data sets first
#
ksVeto<- NULL
for(i in Bands) {
	#
	# Read ks_hist here as it takes too much RAM otherwise
	#
	fn<-paste(prefix, "ks_hist", ".", i, suffix, sep="")
	cat("Loading data from ", fn, "\n", sep="")
	data<-read.table(pipe(paste("cut -f 1,2,3,4 -d \\  < ", fn, sep="")), header=FALSE)
	if(is.null(ksVeto)) {
		ksVeto <- data[,4]<ksLarge
		} else {
		ksVeto <- ksVeto + (data[,4]<ksLarge)
		}
	#
	# Free RAM
	#
	rm(data)
	#gc(verbose=FALSE, reset=TRUE)
	}


fn<-paste(prefix, "hist_residuals", suffix, sep="")
cat("Loading data from ", fn, "\n", sep="")
data<-read.table(pipe(paste("cut -f 1,2,3,4 -d \\  < ", fn, sep="")), header=FALSE)
highRes<-data[,4]
highResMax<-highRes>ResLarge
#
# Free RAM
#
rm(data)
#gc(verbose=FALSE, reset=TRUE)

for(field in c("band", "cputime", "max_dx", "masked", "npoints", 
		"max_high_ul", "max_circ_ul", "median", "TMedian",
		"spindown")) {
	fn<-paste(prefix, field, suffix, sep="")
	cat("Loading data from ", fn, "\n", sep="")
	assign(paste(var_prefix, field, sep=""), read.table(fn))
	}

for(field in c("max_band", "masked_max_band", "max_high_ul_band",
		"max_circ_ul_band")) {
	for(i in Bands) {
		fn<-paste(prefix, field, ".", i, suffix, sep="")
		cat("Loading data from ", fn, "\n", sep="")
		assign(paste(var_prefix,field,i,sep=""), read.table(fn))
		}
	}

plus<-max_dx[,4]=="plus"
cross<-max_dx[,4]=="cross"
pol1<-max_dx[,4]=="pi_1_8"
pol2<-max_dx[,4]=="pi_3_8"
circ<-max_dx[,4]=="circular"

nomasked<-masked[,3]==0
# Assume 1 Hz within 60 Hz is 
no60hz <- (((band[,4]+1)%% 60)>2) & (((band[,4]+1.25)%% 60)>2)

ksVetoP <- pmax(ksVeto[plus], ksVeto[cross], ksVeto[pol1], ksVeto[pol2], ksVeto[circ])==6

highulfactor <- 1.0/cos(pi/8)
circulfactor <- 1.0

# All sky limits
pulsarUL<-max_high_ul[,4]*highulfactor
circularUL<-max_circ_ul[,4]*circulfactor
# Good bands only
pulsarUL <- NULL
circularUL <- NULL
DxMaxBand <- NULL
for(i in Bands) {
	if(is.null(pulsarUL)) {
		pulsarUL<-vn("max_high_ul_band", i)[,3]
		circularUL<-vn("max_band", i)[circ, 7]
		DxMaxBand <- vn("max_band", i)[,6]
		} else {
		pulsarUL<-pmax(pulsarUL, vn("max_high_ul_band", i)[,3])
		circularUL<-pmax(circularUL, vn("max_band", i)[circ, 7])
		DxMaxBand <- pmax(DxMaxBand, vn("max_band", i)[,6])
		}
	}

pulsarUL <- highulfactor*pulsarUL
circularUL <- circulfactor*circularUL

DxMaxBandP <- pmax(DxMaxBand[plus], DxMaxBand[cross], 
	DxMaxBand[pol1], DxMaxBand[pol2], DxMaxBand[circ]);
DxMaxBandPCapped <- pmin(DxMaxBandP, DxCap)
highDxAny <- DxMaxBandP > DxLarge

dir.create(output_dir)	
PrevDevice <- getOption("device")
#
# Compute time estimates
#
#plot(band[,4],cputime[,3], xlab="Frequency", ylab="seconds", main="Compute time")
start_plot("timing")
plot(band[,4],cputime[,3]/3600.0, xlab="Frequency", ylab="hours", main="Compute time")
dev.off()

CPUDays<-sum(cputime[,3]/(3600.0*5.0*24.0))
#
#plot(band[,4],masked[,3], xlab="Frequency", ylab="Count", main="Number of masked points")


#
# Histogram of detected strength
#
Factor<-no60hz[plus] & ksVetoP & !highResMax[plus]
start_plot("detection_strength_hist")
plot(hist(pmin(DxMaxBandPCapped, 15)[Factor], 200), main="Histogram of capped detection strength", xlab="Detection strength capped at 15", ylab="Count")
dev.off()


FancyPlot <- function(UL, f0=NULL, f1=NULL, ylab="y", title=ylab, median_shift=2, median_scale=1) {
	if(is.null(f0)) {
		f0<- min(band[plus, 4])
		}
	if(is.null(f1)) {
		f1<- max(band[plus, 4])
		}
		
	plot(band[plus,4], UL, type="n", main=title, xlim=c(f0, f1), xlab="Frequency", ylab=ylab)
	grid()

	points(band[plus,4], median_scale*(log10(2.0*median[plus,3])*0.5-log10(1800*16384))+median_shift, col="pink", type="s")

	points(band[plus & !no60hz,4], UL[!no60hz[plus]], col="light blue", pch=19)
	factor=no60hz[plus]

	points(band[plus,4][factor & !ksVetoP], UL[factor & !ksVetoP], col="blue", pch=20)
	factor=factor & ksVetoP

	points(band[plus,4][factor & highResMax[plus]], UL[factor & highResMax[plus]], col="red", pch=23)
	factor=factor & !highResMax[plus]

	points(band[plus,4][factor & !highDxAny], UL[factor & !highDxAny], col="dark green", type="s")

	points(band[plus,4][factor & highDxAny], UL[factor & highDxAny], col="red", pch=20)
	}

# Plot of general limits on h0
start_plot("h0UL")
FancyPlot(log10(pulsarUL[plus]), title=h0ULtitle, ylab="Log10 Strain")
dev.off()
#
# Plot of most sensitive polarization - circular
start_plot("circUL")
FancyPlot(log10(circularUL), title=circULtitle, ylab="Log10 Strain")
dev.off()
#
# Detection strength plot, over all polarizations sampled
start_plot("max_dx")
FancyPlot(DxMaxBandPCapped, title=DXtitle, ylab="Capped detection strength", median_shift=370, median_scale=10)
dev.off()

#
# Write out special points data
#

p_max_dx<-pmax(max_dx[cross,5], 
		max_dx[plus,5], 
		max_dx[pol1,5], 
		max_dx[pol2,5], 
		max_dx[circ,5])
		
p_max_dx_pol<-max_dx[cross,]
p_max_dx_pol[p_max_dx==max_dx[plus,5],]<-max_dx[plus,][p_max_dx==max_dx[plus,5],]
p_max_dx_pol[p_max_dx==max_dx[pol1,5],]<-max_dx[pol1,][p_max_dx==max_dx[pol1,5],]
p_max_dx_pol[p_max_dx==max_dx[pol2,5],]<-max_dx[pol2,][p_max_dx==max_dx[pol2,5],]
p_max_dx_pol[p_max_dx==max_dx[circ,5],]<-max_dx[circ,][p_max_dx==max_dx[circ,5],]

signalX<-cos(p_max_dx_pol[,2])*cos(p_max_dx_pol[,3])
signalY<-sin(p_max_dx_pol[,2])*cos(p_max_dx_pol[,3])
signalZ<-sin(p_max_dx_pol[,3])

signalCos<-BandAxis[1]*signalX+BandAxis[2]*signalY+BandAxis[3]*signalZ

Factor <- no60hz[plus] & !ksVetoP
BluePoints<-data.frame(Band=band[plus,4][Factor],
	pulsarUL=pulsarUL[plus][Factor],
	circularUL=circularUL[Factor],
	resMax=highRes[plus][Factor],
	fdot=spindown[plus,3][Factor],
	maxDx=p_max_dx[Factor],
	RA=p_max_dx_pol[Factor, 2],
	DEC=p_max_dx_pol[Factor, 3],
	f0=p_max_dx_pol[Factor, 8],
	ul=p_max_dx_pol[Factor, 6],
	pol=p_max_dx_pol[Factor,4],
	thetaCos=signalCos[Factor]
	)

file<-ofn("BluePoints.txt")
write.table(BluePoints, file)
close(file)



Factor<-no60hz[plus] & ksVetoP & highResMax[plus]
RedDiamonds<-data.frame(Band=band[plus,4][Factor],
	pulsarUL=pulsarUL[plus][Factor],
	circularUL=circularUL[Factor],
	resMax=highRes[plus][Factor],
	fdot=spindown[plus,3][Factor],
	maxDx=p_max_dx[Factor],
	RA=p_max_dx_pol[Factor, 2],
	DEC=p_max_dx_pol[Factor, 3],
	f0=p_max_dx_pol[Factor, 8],
	ul=p_max_dx_pol[Factor, 6],
	pol=p_max_dx_pol[Factor,4],
	thetaCos=signalCos[Factor]
	)

file<-ofn("RedDiamonds.txt")
write.table(RedDiamonds, file)
close(file)

Factor<-no60hz[plus] & ksVetoP & !highResMax[plus] & highDxAny

RedPoints<-data.frame(Band=band[plus,4][Factor],
	pulsarUL=pulsarUL[plus][Factor],
	circularUL=circularUL[Factor],
	resMax=highRes[plus][Factor],
	fdot=spindown[plus,3][Factor],
	maxDx=p_max_dx[Factor],
	RA=p_max_dx_pol[Factor, 2],
	DEC=p_max_dx_pol[Factor, 3],
	f0=p_max_dx_pol[Factor, 8],
	ul=p_max_dx_pol[Factor, 6],
	pol=p_max_dx_pol[Factor,4],
	thetaCos=signalCos[Factor]
	)

file<-ofn("RedPoints.txt")
write.table(RedPoints, file)
close(file)
