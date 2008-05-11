source("params.R")

require("lattice")
require("RMySQL")

p<-function(...) {
	return(paste(..., sep=""))
	}

#SkyBandNamesField<-"grid_points"
SkyBandNamesField<-"max_dx"

FieldsUsed<-c("band", "spindown", "hist_residuals.3", "hist_residuals.4", "median",
		"max_dx.1", "max_dx.2", "max_dx_pol_0", "max_dx.4",
		"max_dx.5", "max_dx.7", "masked", "cputime")

for(skyband in 0:(NBands-1)) {
	FieldsUsed<-c(FieldsUsed, 
		p("ks_hist_", skyband, "_3"), 
		p("max_band_", skyband, "_5"),
		p("max_band_", skyband, "_6"),
		p("max_high_ul_band_", skyband, "_2"),
		p("max_ratio_", skyband, "_3")
		)
	}

NAStrings<-c("NA", "NaN", "NAN")

vn<-function(field, i, position) {
	return(data[,p(var_prefix, field, "_", i, "_", position)])
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}

if(0) {
	fn<-p(prefix, "stat", suffix)
	cat("Loading data from ", fn, "\n", sep="")
	# Cheat - get a small sample first and use it to find out types and such.
	data<-read.table(pipe(p("head --lines=1001 ", fn)), header=TRUE, na.strings=NAStrings, sep=" ")
	cat("Found header with", dim(data)[2], "columns\n")
	
	Types<-lapply(data, class)
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

cat("Connecting to the database\n")
con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

cat("Loading table", DataSet, "\n")
header<-dbGetQuery(con, p("SELECT * FROM ", DataSet, " LIMIT 10"))

BandLine<-unique(dbGetQuery(con, p("SELECT ", p(grep(p("^", SkyBandNamesField, "_sky_band_name"), names(header), value=TRUE), collapse=", "), " FROM ", DataSet, " LIMIT 10")))

NBands<- dim(BandLine)[2]

BandData<-data.frame(i=1:(NBands+length(CompositeBands)), idx=NA, Name=NA)
k<-1
for(i in 0:(NBands-1)) {
	BandData[k, 'idx']<-i
	BandData[k, 'Name']<- as.character(BandLine[1, p(SkyBandNamesField, "_sky_band_name_", i)])
	k<-k+1
	}

for(i in NBands:(NBands+length(CompositeBands)-1)) {
	BandData[k, 'idx']<- i
	BandData[k, 'Name']<- names(CompositeBands)[i-NBands+1]
	k<-k+1
	}

rownames(BandData)<-BandData[,'Name']


Fields<-c("band", "spindown", "hist_residuals_max", "masked", "cputime", "median")
for(i in 0: (NBands-1)) {
	MAXFields<-c(p("max_dx_snr_", i),
		p("max_dx_freq_", i),
		p("ks_hist_", i, "_3"),
		p("max_high_ul_band_", i, "_2"),
		p("max_band_", i, "_5"),
		p("max_band_", i, "_6"),
		p("max_ratio_", i, "_3")
		)

	SUMFields<-c(p("grid_points_", i, "_3")
			)

# 		p("max_dx_pol_", i),
# 		p("max_dx_ra_", i),
# 		p("max_dx_dec_", i),
	Fields<- c(Fields, 
			p("MAX(", MAXFields, ") as ", MAXFields),
			p("SUM(", SUMFields, ") as ", SUMFields),
			p("MAX(IF(max_band_", i, "_4='circular', max_band_", i, "_6, NULL)) as max_circ_ul_band_", i, "_2")
		)
	}

#data<-dbGetQuery(con, p("SELECT ", p(Fields, collapse=", ")," FROM ", DataSet, " WHERE Band>50 GROUP BY band, spindown ORDER BY band, spindown"))
data<-dbGetQuery(con, p("SELECT ", p(Fields, collapse=", ")," FROM ", DataSet, " GROUP BY band, spindown"))

data[,'total_grid_points']<- 0
for(i in 1:NBands) {
	data[,'total_grid_points']<-data[,'total_grid_points']+data[, p("grid_points_", BandData[i, 'idx'], "_3")]
	}

BandVars<-data.frame(prefix=c("max_dx_snr_", "max_dx_freq_", "ks_hist_", "max_high_ul_band_", "max_circ_ul_band_", "max_band_", "max_band_", "max_ratio_"),
		     suffix=c("", "", "_3", "_2", "_2", "_5", "_6", "_3"))

FormVar<-function(k, i) return(p(BandVars[k, 'prefix'], i, BandVars[k, 'suffix']))

for(band in names(CompositeBands)) {
	for(k in 1:dim(BandVars)[1]) {
		data[,FormVar(k, BandData[band, 'idx'])]<-do.call(pmax, data[,FormVar(k, BandData[CompositeBands[[band]],'idx'])])
		}
	data[, p("grid_points_", BandData[band, 'idx'], "_3")]<-do.call(`+`, data[, p("grid_points_", BandData[CompositeBands[[band]], 'idx'], "_3")])
	}


highRes<-data[,'hist_residuals_max']
highResMax<-highRes>ResLarge

nomasked<-data[,'masked']==0
# Assume 1 Hz within 60 Hz is 
no60hz <- (((data[,'band']+1)%% 60)>2) & (((data[,'band']+1.25)%% 60)>2)

Freq<-data[, 'band']
no60hzReduced<-no60hz
highResMaxReduced<-highResMax

	
FancyPlot <- function(UL, f0=NULL, f1=NULL, ylab="y", title=ylab, median_shift=2, median_scale=1) {
	cat("Plotting \"", title, "\"\n", sep="")

	if(is.null(f0)) {
		f0<- min(Freq, na.rm=TRUE)
		}
	if(is.null(f1)) {
		f1<- max(Freq, na.rm=TRUE)
		}

	plot(Freq, UL, type="n", main=title, xlim=c(f0, f1), xlab="Frequency", ylab=ylab)
	grid()

	#points(Freq, median_scale*(log10(2.0*data[plus,'median'])*0.5-log10(1800*16384))+median_shift, col="pink", type="s")

	points(Freq[!no60hzReduced], UL[!no60hzReduced], col="light blue", pch=19)
	Factor<-no60hzReduced

	points(Freq[Factor & !ksVetoP], UL[Factor & !ksVetoP], col="blue", pch=20)
	Factor<-Factor & ksVetoP

	points(Freq[Factor & highResMaxReduced], UL[Factor & highResMaxReduced], col="red", pch=23)
	Factor<-Factor & !highResMaxReduced

	points(Freq[Factor & !highDxAny], UL[Factor & !highDxAny], col="dark green", type="s")

	points(Freq[Factor & highDxAny], UL[Factor & highDxAny], col="red", pch=20)
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

dir.create(output_dir)	

#
# Output raw data file
#
write.csv(data, paste(output_dir, "/raw_data.csv", sep=""), col.names=TRUE, row.names=FALSE)

#
# Describe loaded dataset
# 

Log<-file(paste(output_dir, "/report.log", sep=""), "w")
cat("Dataset loaded from table", DataSet, "\n", file=Log)


BandCounts<-aggregate(data.frame(Count=rep(1, dim(data)[1])), data[,'band', drop=FALSE], sum)
BandCounts[,'Incomplete']<- BandCounts[,'Count']<max(BandCounts[,'Count'], na.rm=TRUE)
BandCounts[,'band']<-as.numeric(as.character(BandCounts[,'band']))
BandCounts<-BandCounts[order(BandCounts[,'band']),,drop=FALSE]
N<-dim(BandCounts)[1]
Step<-median(BandCounts[2:N,'band']-BandCounts[1:(N-1), 'band'])
AllBands<-seq(from=BandCounts[1,'band'], to=BandCounts[N,'band'], by=Step)

cat("First band:", BandCounts[1,'band'], "\n", file=Log)
cat("Last band:", BandCounts[N,'band'], "\n", file=Log)
cat("The following bands are incomplete:", paste(BandCounts[BandCounts[,'Incomplete'], 'band'], collapse=", "), "\n", file=Log)
X<-paste(setdiff(AllBands, BandCounts[,'band']), collapse=", ")
if(X=="")X<-"NONE"
cat("The following bands are missing:", X, "\n", file=Log)

LogFiles<-dbGetQuery(con, p("SELECT DISTINCT LogFile FROM ", DataSet))
Instances<-as.integer(gsub(".*/output/([^/]*)/.*", "\\1", LogFiles[,1]))

cat("First instance:", min(Instances), "\n", file=Log)
cat("Last instance:", max(Instances), "\n", file=Log)
Missing<-sort(setdiff(seq(from=min(Instances), to=max(Instances), by=1), Instances))
X<- paste(Missing, collapse=", ")
if(X=="")X<-"NONE"
cat("The following instances are missing:", X, "\n", file=Log)

IncompleteLogFiles<-dbGetQuery(con, p("SELECT DISTINCT LogFile FROM ", DataSet, " WHERE cputime IS NULL"))
Incomplete<-as.integer(gsub(".*/output/([^/]*)/.*", "\\1", IncompleteLogFiles[,1]))
X<- paste(Incomplete, collapse=", ")
if(X=="")X<-"NONE"
cat("The following instances are incomplete:", X, "\n", file=Log)

Todo<-sort(unique(c(Missing, Incomplete)))

if(length(Todo)>0) {
	cat("--------- todo.dag -------------\n", file=Log)
	cat(paste("JOB A", Todo, " condor\nVARS A", Todo, " PID=\"", Todo, "\"\n", sep="", collapse="\n"), file=Log)
	cat("-----------------------------------\n", file=Log)
	}


close(Log)

#
# Compute time estimates
#
#plot(band[,4],cputime[,3], xlab="Frequency", ylab="seconds", main="Compute time")

skyband<-0
BandName<- ""

start_plot("HighResHist")
plot(hist(highRes, 50))
dev.off()

start_plot("HighRes")
plot(data[,'band'], highRes)
points(data[,'band'], rep(ResLarge, dim(data)[1]), col="red", type="l")
dev.off()

start_plot("timing")
plot(Freq, data[,'cputime']/3600.0, xlab="Frequency", ylab="hours", main="Compute time")
dev.off()

#start_plot("instance_timing")
#plot(data[,"Instance"], data[,'cputime']/3600.0, xlab="Instance", ylab="hours", main="Compute time")
#dev.off()

#CPUDays<-sum(cputime[,3]/(3600.0*5.0*24.0))
#
#plot(band[,4],masked[,3], xlab="Frequency", ylab="Count", main="Number of masked points")

		

for(i in 1:dim(BandData)[1]) {
	skyband <- BandData[i, 'idx']
	BandName <- BandData[i, 'Name']
	cat("Processing skyband", BandName, "\n")
	#
	# Load large data sets first
	#

	#
	# Read ks_hist here as it takes too much RAM otherwise
	#
	ksVetoP <- data[,p("ks_hist_", skyband, "_3")]<ksLarge
	#data[,'max_dx_pol_0']
		
	highulFactor <- 2.0/cos(pi/8)
	circulFactor <- 1.0
	
	# All sky limits
	#pulsarUL<-max_high_ul[,4]*highulFactor
	#circularUL<-max_circ_ul[,4]*circulFactor

	pulsarUL<-vn("max_high_ul_band", skyband, 2)
	circularUL<-vn("max_circ_ul_band", skyband, 2)
	DxMaxBand <- vn("max_band", skyband, 5)
	
	pulsarUL <- highulFactor*pulsarUL
	circularUL <- circulFactor*circularUL

	DxMaxBand <- data[, p("max_dx_snr_", skyband)]
	DxMaxBandCapped <- pmin(DxMaxBand, DxCap)
	highDxAny <- DxMaxBand > DxLarge
	
	if(sum(!is.na(pulsarUL))<1) {
		cat("no data for skyband", skyband, "\n");
		next
		}

	write.csv(data.frame(
		Frequency=Freq,
		Spindown=data[,'spindown'],
		SkyAreaPct=100*data[, p("grid_points_", skyband, "_3")]/data[,"total_grid_points"],
		UpperLimit=pulsarUL,
		CircUpperLimit=circularUL,
		Line60=as.integer(!no60hz),
		LineWandering=as.integer(highResMaxReduced),
		KSVeto= as.integer(!ksVetoP)),
			paste(output_dir, "/data_", BandName, ".csv", sep=""), col.names=TRUE, row.names=FALSE)
	
	#
	# Histogram of detected strength
	#
	Factor<-no60hz & ksVetoP & !highResMax
	#start_plot("detection_strength_hist")
	#plot(hist(pmin(DxMaxBandCapped, 15)[Factor], 200), main="Histogram of capped detection strength", xlab="Detection strength capped at 15", ylab="Count")
	#dev.off()
	
	
	start_plot("grid_points")
	plot(data[,'band'], vn("grid_points", skyband, 3), ylab="Grid points", xlab="Frequency")
	dev.off()

	start_plot("grid_points_pct")
	plot(data[,'band'], 100.0*vn("grid_points", skyband, 3)/data[,'total_grid_points'], ylab="Sky area (percent)", xlab="Frequency")
	dev.off()

	# Plot of general limits on h0
	start_plot("h0UL")
	FancyPlot(log10(pulsarUL), title=h0ULtitle, ylab="Log10 Strain")
	dev.off()
	#
 	start_plot("circUL")
 	FancyPlot(log10(circularUL), title=circULtitle, ylab="Log10 Strain")
 	dev.off()
	#
	# Detection strength plot, over all polarizations sampled
	start_plot("max_dx")
	FancyPlot(DxMaxBandCapped, title=DXtitle, ylab="Capped detection strength", median_shift=370, median_scale=10)
	dev.off()
	
	
	#
	# 200 Hz chunks
	#
	f0<-floor(min(data[,'band'], na.rm=TRUE))
	f1<-max(data[,'band'], na.rm=TRUE)
	while(f0<f1){
		start_plot(paste("max_dx_",f0,sep=""))
		FancyPlot(DxMaxBandCapped, title=DXtitle, ylab="Capped detection strength", median_shift=370, median_scale=10, f0=f0,f1=f0+200.0)
		dev.off()
	
		start_plot(paste("max_dx_map_",f0,sep=""))
		FancyMap(DxMaxBandCapped, title=DXtitle, levels=c(min(DxMaxBandCapped, na.rm=TRUE),7,10,15,20,30,50.0), f0=f0,f1=f0+200.0)
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
	
	# start_plot("max_dx_map")
	# FancyMap(DxMaxBandCapped, title="Detection strength", levels=c(min(DxMaxBandCapped),7,10,15,20,30,50))
	# dev.off()
	
	start_pdf_plot("max_dx_map", title="Detection strength")
	FancyMap(DxMaxBandCapped, title="Detection strength", levels=c(min(DxMaxBandCapped, na.rm=TRUE),7,10,15,20,30,50))
	dev.off()
	
	start_pdf_plot("max_ratio_map")
	FancyMap(vn("max_ratio", skyband,3), title="Maximum veto ratio (by weight)")
	dev.off()


	next
	#
	# Write out special points data
	#
	
	p_max_dx<-pmax(data[cross, 'max_dx.4'], 
			data[plus,'max_dx.4'], 
			data[pol1,'max_dx.4'], 
			data[pol2,'max_dx.4'], 
			data[circ,'max_dx.4'])
			
	Fields<-c('max_dx.1','max_dx.2','max_dx_pol_0', 'max_dx.5', 'max_dx.7')
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

	BluePoints<-data.frame(Band=data[plus,'band'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor, 'max_dx_pol_0'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("BluePoints_", skyband, "_txt", sep=""))
	write.table(BluePoints, file)
	close(file)
	
	
	
	Factor<-no60hz[plus] & ksVetoP & highResMax[plus]
	Factor[is.na(Factor)]<-FALSE

	RedDiamonds<-data.frame(Band=data[plus,'band'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor, 'max_dx_pol_0'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("RedDiamonds_", skyband, "_txt", sep=""))
	write.table(RedDiamonds, file)
	close(file)
	
	Factor<-no60hz[plus] & ksVetoP & !highResMax[plus] & highDxAny
	Factor[is.na(Factor)]<-FALSE
	
	RedPoints<-data.frame(Band=data[plus,'band'][Factor],
		pulsarUL=pulsarUL[plus][Factor],
		circularUL=circularUL[Factor],
		resMax=highRes[plus][Factor],
		fdot=data[plus,'spindown'][Factor],
		maxDx=p_max_dx[Factor],
		RA=p_max_dx_pol[Factor, 'max_dx.1'],
		DEC=p_max_dx_pol[Factor, 'max_dx.2'],
		f0=p_max_dx_pol[Factor, 'max_dx.7'],
		ul=p_max_dx_pol[Factor, 'max_dx.5'],
		pol=p_max_dx_pol[Factor,'max_dx_pol_0'],
		thetaCos=signalCos[Factor]
		)
	
	file<-ofn(paste("RedPoints_", skyband, "_txt", sep=""))
	write.table(RedPoints, file)
	close(file)

	# Write out info for all bands
	AllPoints<-data.frame(Band=data[plus, 'band'],
		pulsarUL=pulsarUL[plus],
		circularUL=circularUL,
		resMax=highRes[plus],
		fdot=data[plus, 'spindown'],
		maxDx=p_max_dx,
		RA=p_max_dx_pol[, 'max_dx.1'],
		DEC=p_max_dx_pol[, 'max_dx.2'],
		f0=p_max_dx_pol[, 'max_dx.7'],
		ul=p_max_dx_pol[, 'max_dx.5'],
		pol=p_max_dx_pol[,'max_dx_pol_0'],
		thetaCos=signalCos
		)

	file<-ofn(paste("AllPoints_", skyband, "_txt", sep=""))
	write.table(AllPoints, file)
	close(file)

	}
