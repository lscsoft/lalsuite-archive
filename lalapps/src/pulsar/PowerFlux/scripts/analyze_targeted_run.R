source("params.R")

require("lattice")
require("RMySQL")

p<-function(...) {
	return(paste(..., sep=""))
	}

fn<-p(prefix, "stat", suffix)

LoadData<- FALSE
try({ if(length(header)>0)LoadData<-FALSE}, silent=TRUE)

NAStrings<-c("NA", "NaN", "NAN")

vn<-function(field, i, position) {
	return(data[,p(var_prefix, field, "_", i, "_", position)])
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}


if(LoadData) {
	cat("Loading data from ", fn, "\n", sep="")
	# Cheat - get a small sample first and use it to find out types and such.
	header<-read.table(pipe(p("head --lines=1001 ", fn)), header=TRUE, na.strings=NAStrings, sep=" ")
	cat("Found header with", dim(header)[2], "columns\n")
	
	Types<-lapply(header, class)
	Types<-lapply(Types, function(x)gsub("logical", "numeric", x))

	BandNames<-grep("grid_points.*2$", colnames(header), value=TRUE)
	
	NBands<- length(BandNames)
	
	BandData<-data.frame(i=1:(NBands/2), idx=NA, Name=NA, ref_idx=NA, ref_Name=NA)
	k<-1
	for(i in 1:NBands) {
		if(regexpr("_ref$", as.character(header[1, p("grid_points.", i-1, ".2")]))>=0)next
		BandData[k, 'idx']<-i-1
		BandData[k, 'Name']<- as.character(header[1, p("grid_points.", i-1, ".2")])
		for(m in 1:NBands) {
			if(p(header[1, p("grid_points.", i-1, ".2")], "_ref")==header[1, p("grid_points.", m-1, ".2")]) {
				BandData[k, 'ref_idx']<-m-1
				BandData[k, 'ref_Name']<- as.character(header[1, p("grid_points.", m-1, ".2")])
				break
				}
			}
		k<-k+1
		}

	FieldsUsed<-c("band.3", "spindown.2", "hist_residuals.3", "hist_residuals.4", "median.2",
			"max_dx.1", "max_dx.2", "max_dx.3", "max_dx.4",
			"max_dx.5", "max_dx.7", "masked.2", "cputime.2")
	
	for(skyband in 0:(NBands-1)) {
		FieldsUsed<-c(FieldsUsed, 
			p("max_dx_band.", skyband, ".1"), 
			p("max_dx_band.", skyband, ".2"), 
			p("max_dx_band.", skyband, ".3"), 
			p("max_dx_band.", skyband, ".4"), 
			p("max_dx_band.", skyband, ".5"), 
			p("max_dx_band.", skyband, ".6"), 
			p("max_dx_band.", skyband, ".7"), 
			p("max_dx_band.", skyband, ".8"), 
			p("ks_hist.", skyband, ".3"), 
			p("max_band.", skyband, ".5"),
			p("max_band.", skyband, ".6"),
			p("max_high_ul_band.", skyband, ".2"),
			p("max_ratio.", skyband, ".3")
			)
		}
	
	
	#Types<-lapply(Types, function(x) {
	#	if(x=="factor")return("character")
	#	return(x)
	#	})
	FieldsUnused<-setdiff(names(Types), FieldsUsed)
	Types[FieldsUnused]<-NULL
	
	data<-read.table(fn, sep=" ", header=TRUE, na.strings=NAStrings, colClasses=Types, as.is=TRUE)
	cat("Loaded table with", dim(data)[2], "columns and", dim(data)[1], "rows\n")
	#data<-read.table(pipe(paste("cut -f 1,2,3,4 -d \\  < ", fn, sep="")), sep="\t", header=TRUE, na.strings=NAStrings)
	cat("Reducing to only columns we want\n")
	data<-data[, FieldsUsed]
	gc()
	gc(reset=TRUE)
	}

# # Prune out extra headers
# data<-data[data[,'LogFile']!='LogFile', ]
# 
# # Save and reload - this properly converts numbers
# cat("Saving\n")
# write.table(data, "data.csv", sep="\t", row.names=FALSE)
# cat("Reloading\n")
# data<-read.table("data.csv", sep="\t", header=TRUE)

# Sort - as log files can come in any order
# cat("Sorting\n")
# P<-order(data[,"band.3"], data[,"spindown.2"], data[,"max_dx.3"])
# data<-data[P,]
# # Explicitly delete old data..
# rm(P)
gc()
# gc(reset=TRUE)
#mem.limits(vsize=1024^3)

cat("Connecting to the database\n")
con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

cat("Loading table", DataSet, "\n")
header<-dbGetQuery(con, p("SELECT * FROM ", DataSet, " LIMIT 10"))

BandLine<-unique(dbGetQuery(con, p("SELECT ", p(grep("^max_dx_sky_band_name", names(header), value=TRUE), collapse=", "), " FROM ", DataSet, " LIMIT 10")))

NBands<- dim(BandLine)[2]

BandData<-data.frame(i=1:(NBands/2), idx=NA, Name=NA, ref_idx=NA, ref_Name=NA)
k<-1
for(i in 1:NBands) {
	if(regexpr("_ref$", as.character(BandLine[1, i]))>=0)next
	BandData[k, 'idx']<-i-1
	BandData[k, 'Name']<- as.character(BandLine[1, i])
	for(m in 1:NBands) {
		if(p(BandLine[1, i], "_ref")==BandLine[1, m]) {
			BandData[k, 'ref_idx']<-m-1
			BandData[k, 'ref_Name']<- as.character(BandLine[1, m])
			break
			}
		}
	k<-k+1
	}

Fields<-c("band", "spindown", "hist_residuals_max", "masked", "cputime")
for(i in 0: (NBands-1)) {
	Fields<-c(Fields, 
		p("max_dx_snr_", i),
		p("max_dx_freq_", i),
		p("max_dx_pol_", i),
		p("max_dx_ra_", i),
		p("max_dx_dec_", i),
		p("MAX(ks_hist_", i, "_3) as ks_hist_", i, "_3"),
		p("max_high_ul_band_", i, "_2"),
		p("max_ratio_", i, "_3")
		)
	}

data<-dbGetQuery(con, p("SELECT ", p(Fields, collapse=", ")," FROM ", DataSet, " WHERE Band>50 GROUP BY band, spindown ORDER BY band, spindown"))

highRes<-data[,'hist_residuals_max']
highResMax<-highRes>ResLarge
BandName<-""
dir.create(output_dir)



# plus<-data[,'max_dx.3']=="plus"
# cross<-data[,'max_dx.3']=="cross"
# pol1<-data[,'max_dx.3']=="pi_1_8"
# pol2<-data[,'max_dx.3']=="pi_3_8"
# circ<-data[,'max_dx.3']=="circular"

nomasked<-data[,'masked']==0
# Assume 1 Hz within 60 Hz is 
no60hz <- (((data[,'band']+1)%% 60)>2) & (((data[,'band']+1.25)%% 60)>2)

Freq<-data[, 'band']
no60hzReduced<-no60hz
highResMaxReduced<-highResMax

start_plot("HighResHist")
plot(hist(highRes, 50))
dev.off()

start_plot("HighRes")
plot(data[,'band'], highRes)
points(data[,'band'], rep(ResLarge, dim(data)[1]), col="red", type="l")
dev.off()

	
FancyPlot <- function(UL, f0=NULL, f1=NULL, ylab="y", title=ylab, median_shift=2, median_scale=1, skip.high.dx=FALSE, summary.curve=NULL) {
	cat("Plotting \"", title, "\"\n", sep="")

	if(is.null(f0)) {
		f0<- min(Freq, na.rm=TRUE)
		}
	if(is.null(f1)) {
		f1<- max(Freq, na.rm=TRUE)
		}

	plot(Freq, UL, type="n", main=title, xlim=c(f0, f1), xlab="Frequency", ylab=ylab)
	grid()

	#points(Freq, median_scale*(log10(2.0*data[,'median.2'])*0.5-log10(1800*16384))+median_shift, col="pink", type="s")

	points(Freq[!no60hzReduced], UL[!no60hzReduced], col="magenta", pch="*")
	Factor=no60hzReduced

	points(Freq[Factor & !ksVetoP], UL[Factor & !ksVetoP], col="blue", pch=20)
	Factor=Factor & ksVetoP

	points(Freq[Factor & highResMaxReduced], UL[Factor & highResMaxReduced], col="red", pch=23)
	Factor=Factor & !highResMaxReduced


	if(skip.high.dx) {
		points(Freq[Factor], UL[Factor], col="dark green", pch="+")
		} else {
		points(Freq[Factor & !highDxAny], UL[Factor & !highDxAny], col="dark green", pch="+")

		points(Freq[Factor & highDxAny], UL[Factor & highDxAny], col="red", pch=20)
		}

	if(!is.null(summary.curve)) {
		points(Freq, summary.curve, col="cyan", type="l")
		}

	}

FancyMap <- function(UL, title="", levels=NULL, f0=NULL, f1=NULL) {
	sp<-data[,'spindown']
	ba<-Freq

	XX<-sort(unique(sp))
	YY<-unique(ba)
	if(!is.null(f0)) {
		YY <- YY[YY>=f0]
		}
	if(!is.null(f1)) {
		YY <- YY[YY<=f1]
		}
	YY<-sort(YY)


	NSpindown<-length(XX)
	if(NSpindown<2) {
		cat("Skipping map - spindown count is too low\n")
		return(NULL)
		}
	NBand<-length(YY)
	cat("Mapping \"", title, "\" NSpindown=", NSpindown, " NBand=", NBand, "\n", sep="")
	
	#ZZ<-array(NA, c(NSpindown, NBand))

# 	for(i in 1:length(ba)){
# 		X<-which.min((XX-sp[i])^2)
# 		Y<-which.min((YY-ba[i])^2)
# 		if(YY[Y]!=ba[i]) next
# 	
# 		ZZ[X,Y]<-UL[i]
# 		}

	#filled.contour(x=XX, y=YY, z=ZZ, color.palette=topo.colors, levels=levels, xlab="Spindown", ylab="Frequency", main=title)
	#image(x=XX, y=YY, z=ZZ, col=topo.colors(length(levels)-1), breaks=levels, xlab="Spindown", ylab="Frequency", main=title)
	if(is.null(f0)) {
		f0<-min(ba, na.rm=TRUE)
		} 
	if(is.null(f1)) {
		f1<-max(ba, na.rm=TRUE)
		}
	F<-(ba>=f0) & (ba<=f1) & !is.na(UL)
	grid<-data.frame(x=sp, y=ba, z=UL)

	if(is.null(levels)){
		levels<-pretty(range(UL[F]), 20)
		}

	plot(levelplot(z~x*y, grid[F,], col.regions=topo.colors(length(levels)-1), at=levels, 
		xlab="Spindown", ylab="Frequency", main=title))
	}


#
# Compute time estimates
#
#plot(band[,4],cputime[,3], xlab="Frequency", ylab="seconds", main="Compute time")

skyband<-0
BandName<-""

if(sum(!is.na(data[,'cputime']))>0) {
	start_plot("timing")
	plot(Freq, data[,'cputime']/3600.0, xlab="Frequency", ylab="hours", main="Compute time")
	dev.off()
	}

#CPUDays<-sum(cputime[,3]/(3600.0*5.0*24.0))
#
#plot(band[,4],masked[,3], xlab="Frequency", ylab="Count", main="Number of masked points")

		

for(i in 1:(dim(BandData)[1])) {
	skyband<- BandData[i, 'idx']
	ref_skyband<- BandData[i, 'ref_idx']
	BandName<- BandData[i, 'Name']
	cat("Processing skyband", skyband, BandName, "\n")
	#
	# Load large data sets first
	#

	#
	# Read ks_hist here as it takes too much RAM otherwise
	#
	ksVeto <- data[,p("ks_hist_", skyband, "_3")]<ksLarge
	#data[,'max_dx.3']
	
	
	ksVetoP <- ksVeto
	
	highulFactor <- 2.0/cos(pi/8)
	circulFactor <- 1.0
	
	# All sky limits
	#pulsarUL<-max_high_ul[,4]*highulFactor
	#circularUL<-max_circ_ul[,4]*circulFactor

	pulsarUL<-vn("max_high_ul_band", skyband, 2)
	ref_pulsarUL<-vn("max_high_ul_band", ref_skyband, 2)
	#circularUL<-vn("max_band", skyband, 6)[circ]
	DxMaxBand <- data[, p("max_dx_snr_", skyband)]
	ref_DxMaxBand <- data[, p("max_dx_snr_", ref_skyband)]
	
	pulsarUL <- highulFactor*pulsarUL
	ref_pulsarUL <- highulFactor*ref_pulsarUL
	#circularUL <- circulFactor*circularUL

	Diff_pulsarUL<- pmax(0, pulsarUL-ref_pulsarUL)/pulsarUL

	Diff_DxMaxBand <- (pmax(2, DxMaxBand)-pmax(2, ref_DxMaxBand))*15/(10+pmax(2, ref_DxMaxBand))

	DxMaxBandCapped <- pmin(DxMaxBand, DxCap)
	highDxAny <- DxMaxBand > DxLarge
	Diff_DxMaxBandCapped <- pmax(0, Diff_DxMaxBand)
	
	if(sum(!is.na(pulsarUL))<1) {
		cat("no data for skyband", skyband, "\n");
		next
		}
	
	SummaryCurve<-SummaryFunc(Freq)

	#
	# Histogram of detected strength
	#
	#Factor<-no60hz & ksVetoP & !highResMax
	#start_plot("detection_strength_hist")
	#plot(hist(pmin(DxMaxBandPCapped, 15)[Factor], 200), main="Histogram of capped detection strength", xlab="Detection strength capped at 15", ylab="Count")
	#dev.off()
	
	
	# Plot of general limits on h0
	start_plot("h0UL")
	FancyPlot(log10(pulsarUL), title=h0ULtitle, ylab="Log10 Strain", summary.curve=log10(SummaryCurve))
	dev.off()

	start_plot("h0_summary_hist")
	F<- no60hzReduced & ksVetoP
	plot(hist(log10(pulsarUL[F])-log10(SummaryCurve[F]), 100), xlab="log10(h0/summary)", main="Summary curve deviation", ylab="Count")
	dev.off()


	start_plot("diff_h0UL")
	FancyPlot(Diff_pulsarUL, title=h0ULtitle, ylab="Excess strain")
	dev.off()
	#
	# Plot of most sensitive polarization - circular
# 	start_plot("circUL")
# 	FancyPlot(log10(circularUL), title=circULtitle, ylab="Log10 Strain")
# 	dev.off()
	#
	# Detection strength plot, over all polarizations sampled
	start_plot("max_dx")
	FancyPlot(DxMaxBandCapped, title=DXtitle, ylab="Capped detection strength", median_shift=370, median_scale=10)
	dev.off()
	
	start_plot("diff_max_dx")
	FancyPlot(Diff_DxMaxBandCapped, title=DXtitle, ylab="Excess SNR", median_shift=370, median_scale=10)
	dev.off()
	
	#
	# 200 Hz chunks
	#
	f0<-floor(min(Freq, na.rm=TRUE))
	f1<-max(Freq, na.rm=TRUE)
	while(f0<f1) {
		start_plot(paste("max_dx_",f0,sep=""))
		FancyPlot(DxMaxBandCapped, title=DXtitle, ylab="Capped detection strength", median_shift=370, median_scale=10, f0=f0,f1=f0+200.0)
		dev.off()
	
		start_plot(paste("max_dx_map_",f0,sep=""))
		FancyMap(DxMaxBandCapped, title=DXtitle, levels=c(min(DxMaxBandCapped, na.rm=TRUE),4,5,6,7,10,15,20,50.0), f0=f0,f1=f0+200.0)
		dev.off()

		start_plot(paste("diff_max_dx_map_",f0,sep=""))
		FancyMap(Diff_DxMaxBandCapped, title="Diff SNR", levels=c(min(Diff_DxMaxBandCapped, na.rm=TRUE), 0.2, 0.5, 0.7, 1,2,3,5,max(Diff_DxMaxBandCapped, na.rm=TRUE)), f0=f0,f1=f0+200.0)
		dev.off()
	
		start_plot(paste("max_ratio_map_",f0,sep=""))
		FancyMap(vn("max_ratio", skyband, 3),title="Maximum veto ratio (by weight)", f0=f0,f1=f0+200.0)
		dev.off()
	
		f0<-f0+100
		}
	
	#
	# Spindown versus frequency maps
	
	# start_plot("h0ULmap")
	# FancyMap(log10(pulsarUL[plus]), title="Log10 h0 UL")
	# dev.off()
	
	start_pdf_plot("h0ULmap", title="Log10 h0 UL")
	FancyMap(log10(pulsarUL), title="Log10 h0 UL")
	dev.off()

	start_pdf_plot("diff_h0ULmap", title="Excess h0 UL")
	FancyMap(Diff_pulsarUL, title="Excess h0 UL")
	dev.off()
	
	# start_plot("max_dx_map")
	# FancyMap(DxMaxBandPCapped, title="Detection strength", levels=c(min(DxMaxBandPCapped),7,10,15,20,30,50))
	# dev.off()
	
	start_pdf_plot("max_dx_map", title="Detection strength")
	FancyMap(DxMaxBandCapped, title="Detection strength", levels=c(min(DxMaxBandCapped, na.rm=TRUE),4,5,6,7,10,15,20,50))
	dev.off()

	start_pdf_plot("diff_max_dx_map", title="Diff SNR")
	FancyMap(Diff_DxMaxBandCapped, title="Diff SNR", levels=c(min(Diff_DxMaxBandCapped, na.rm=TRUE), 0.2, 0.5, 0.7, 1,2,3,5,max(Diff_DxMaxBandCapped, na.rm=TRUE)))
	dev.off()
	
	start_pdf_plot("max_ratio_map")
	FancyMap(vn("max_ratio", skyband,3), title="Maximum veto ratio (by weight)")
	dev.off()

	F<- Diff_DxMaxBandCapped>DiffDxLarge
	F[is.na(F)]<-FALSE
	
	file<-ofn(paste("Outliers_", BandName, ".txt", sep=""))
	write.table(data[F,], file)
	close(file)

	next
	#
	# Write out special points data
	#
	
	p_max_dx<-pmax(data[cross, 'max_dx.4'], 
			data[plus,'max_dx.4'], 
			data[pol1,'max_dx.4'], 
			data[pol2,'max_dx.4'], 
			data[circ,'max_dx.4'])
			
	Fields<-c('max_dx.1','max_dx.2','max_dx.3', 'max_dx.5', 'max_dx.7')
	p_max_dx_pol<-data[cross, Fields]

	F<-p_max_dx==data[plus,'max_dx.4']
	F[is.na(F)]<-FALSE
	p_max_dx_pol[F,]<-data[plus, Fields][F,]

	F<-p_max_dx==data[cross,'max_dx.4']
	F[is.na(F)]<-FALSE
	p_max_dx_pol[F,]<-data[cross, Fields][F,]

	F<-p_max_dx==data[pol1,'max_dx.4']
	F[is.na(F)]<-FALSE
	p_max_dx_pol[F,]<-data[pol1, Fields][F,]

	F<-p_max_dx==data[pol2,'max_dx.4']
	F[is.na(F)]<-FALSE
	p_max_dx_pol[F,]<-data[pol2, Fields][F,]

	F<-p_max_dx==data[circ,'max_dx.4']
	F[is.na(F)]<-FALSE
	p_max_dx_pol[F,]<-data[circ, Fields][F,]

	signalX<-cos(p_max_dx_pol[, 'max_dx.1'])*cos(p_max_dx_pol[, 'max_dx.2'])
	signalY<-sin(p_max_dx_pol[, 'max_dx.1'])*cos(p_max_dx_pol[, 'max_dx.2'])
	signalZ<-sin(p_max_dx_pol[, 'max_dx.2'])
	
	signalCos<-BandAxis[1]*signalX+BandAxis[2]*signalY+BandAxis[3]*signalZ
	
	Factor <- no60hz[plus] & !ksVetoP
	Factor[is.na(Factor)]<-FALSE

	BluePoints<-data.frame(Band=data[plus,'band.3'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown.2'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor, 'max_dx.3'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("BluePoints_", BandName, ".txt", sep=""))
	write.table(BluePoints, file)
	close(file)
	
	
	
	Factor<-no60hz[plus] & ksVetoP & highResMax[plus]
	Factor[is.na(Factor)]<-FALSE

	RedDiamonds<-data.frame(Band=data[plus,'band.3'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown.2'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor, 'max_dx.3'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("RedDiamonds_", BandName, ".txt", sep=""))
	write.table(RedDiamonds, file)
	close(file)
	
	Factor<-no60hz[plus] & ksVetoP & !highResMax[plus] & highDxAny
	Factor[is.na(Factor)]<-FALSE
	
	RedPoints<-data.frame(Band=data[plus,'band.3'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown.2'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor,'max_dx.3'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("RedPoints_", BandName, ".txt", sep=""))
	write.table(RedPoints, file)
	close(file)

	# Write out info for all bands
	AllPoints<-data.frame(Band=data[plus, 'band.3'],
		pulsarUL=pulsarUL[plus],
		circularUL=circularUL,
		resMax=highRes[plus],
		fdot=data[plus, 'spindown.2'],
		maxDx=p_max_dx,
		RA=p_max_dx_pol[, 'max_dx.1'],
		DEC=p_max_dx_pol[, 'max_dx.2'],
		f0=p_max_dx_pol[, 'max_dx.7'],
		ul=p_max_dx_pol[, 'max_dx.5'],
		pol=p_max_dx_pol[,'max_dx.3'],
		thetaCos=signalCos
		)

	file<-ofn(paste("AllPoints_", BandName, ".txt", sep=""))
	write.table(AllPoints, file)
	close(file)

	}
