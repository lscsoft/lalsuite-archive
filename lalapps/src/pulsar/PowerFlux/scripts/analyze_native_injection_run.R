library("RMySQL")
p<-function(...)paste(..., sep="")

source("params.R")

Data<-read.table(DataFile, header=TRUE)
Params<-read.table(ParameterFile, header=TRUE)

if(!is.null(CandidatesTable)) {
	con<-dbConnect("MySQL", user=User, dbname=Database, host=Host)
	OptCand<-dbGetQuery(con, p("SELECT * FROM ", CandidatesTable))
	} else {
	# Create dataframe with 0 rows
	OptCand<-NULL
	}

LOG<-file("run.log", open="w")
Labels<-unique(Data[,"label.1"])
Missing<- sort(setdiff(1:max(Labels), Labels))
LabelCount<-aggregate(Data[,"cputime.2"], Data[,"label.1", drop=FALSE], function(x)sum(!is.na(x)))
F<-LabelCount[,"x"]!=5
Incomplete<-sort(as.integer(as.character(LabelCount[F, "label.1"])))

cat(file=LOG, length(Missing), "missing instances:", p(Missing, collapse=" "), "\n")
cat(file=LOG, length(Incomplete), "incomplete instances:", p(Incomplete, collapse=" "), "\n")
redo<-sort(c(Missing, Incomplete))
if(length(redo)>0) {
	cat(file=LOG, "----------------------- dag.redo ----------------------\n")
	cat(file=LOG, p("JOB A", redo, " condor\nVARS A", redo, " PID=\"", redo, "\"", collapse="\n"))
	}
close(LOG)


LargeStrain<-Data[,'fake_strain.2']>1e-24
smallSpindown<-Data[,'fake_spindown.2']< SpindownCutoff
highStrain<- Data[,'fake_strain.2']>StrainCutoff

highulFactor <- 1.0/cos(pi/8)
circulFactor <- 1.0

pulsarUL<-Data[,'max_high_ul.3']*highulFactor

injectedStrain<- Data[,'fake_strain.2']

excessStrain<-2.0*pulsarUL-injectedStrain

goodBand<-array(TRUE, dim(Data)[1])
equatorialBand<- abs(Data[,'fake_dec.2'])<0.305

if(!is.null(OptCand)) {
	OptCand[,'i']<-as.integer(OptCand[,'label'])
	
	Merged<-merge(Params, OptCand, by="i", suffixes=c("", "_opt"))
	Merged[,'FDist']<- 1800.0*abs(Merged[,'frequency']-Merged[,'f0'])
	Merged[,'SDist']<- abs(Merged[,'spindown']-Merged[,'spindown_opt'])
	Merged[,'Dist']<- acos(sin(Merged[,'dec'])*sin(Merged[,'dec_opt'])+cos(Merged[,'dec'])*cos(Merged[,'dec_opt'])*cos(Merged[,'ra']-Merged[,'ra_opt']))
	
	Merged[,'Detected']<- (Merged[,'Dist']<0.05)*1.1 + (Merged[,'FDist']<4)*1.2 + (Merged[,'SDist']<1e-10)*1.4
	Merged[,'Detected2']<- (Merged[,'Dist']<0.07)*1.1 + (Merged[,'FDist']<5)*1.2 + (Merged[,'SDist']<2e-10)*1.4
	Merged[,'Detected3']<- (Merged[,'Dist']<0.07)*1.01 + (Merged[,'FDist']<5)*1.02 + (Merged[,'SDist']<2e-10)*1.04+(Merged[,'power_cor']>0.06)*1.08
	
	
	MSmallSp<-Merged[,'spindown']< SpindownCutoff
	
	M_Detected<-merge(Params, aggregate(Merged[,'Detected', drop=FALSE], Merged[,'i', drop=FALSE], max, na.rm=TRUE), by="i")
	M_Detected2<-merge(Params, aggregate(Merged[,'Detected2', drop=FALSE], Merged[,'i', drop=FALSE], max, na.rm=TRUE), by="i")
	M_Detected3<-merge(Params, aggregate(Merged[,'Detected3', drop=FALSE], Merged[,'i', drop=FALSE], max, na.rm=TRUE), by="i")
	
	#
	# 0.05 - good cutoff
	#
	M_MinDist<-merge(Params, aggregate(Merged[,'Dist', drop=FALSE], Merged[,'i', drop=FALSE], min), by="i")
	#
	# 4 - good cutoff
	#
	M_MinFDist<-merge(Params, aggregate(Merged[,'FDist', drop=FALSE], Merged[,'i', drop=FALSE], min), by="i")
	#
	# 5e-11 - 1e-10
	#
	M_MinSDist<-merge(Params, aggregate(Merged[,'SDist', drop=FALSE], Merged[,'i', drop=FALSE], min), by="i")
	#
	# 
	#
	M_MaxPowerCor<-merge(Params, aggregate(Merged[,'power_cor', drop=FALSE], Merged[,'i', drop=FALSE], max), by="i")
	}

start_plot<-function(name) {
	png(filename=name, width=800, height=800, bg="white")
	}

start_plot("cos_iota_dist.png")
plot(hist(cos(Params[,'iota']),50), xlab=expression(cos(iota)), ylab="Count", main="Distribution of cos(iota)")
dev.off()


start_plot("spindown_sensitivity.png")
plot(Data[goodBand, 'fake_spindown.2'], excessStrain[goodBand], xlab="Spindown", ylab="Detected strain-injected strain", main="Spindown sensitivity, skyband 0", log="x")
dev.off()

start_plot("high_strain_spindown_sensitivity.png")
plot(Data[goodBand & highStrain, 'fake_spindown.2'], excessStrain[goodBand & highStrain], xlab="Spindown", ylab="Detected strain-injected strain", main=p("Spindown sensitivity, skyband 0, strain > ", StrainCutoff), log="x")
dev.off()

start_plot("strain_sensitivity.png")
plot(injectedStrain[goodBand & smallSpindown], excessStrain[goodBand & smallSpindown], xlab="Injected strain", ylab="Detected strain-injected strain", main=p("Strain sensitivity, skyband 0, spindown<", SpindownCutoff))
dev.off()

start_plot("strain_sensitivity_equator.png")
plot(injectedStrain[goodBand & smallSpindown & equatorialBand], excessStrain[goodBand & smallSpindown & equatorialBand], xlab="Injected strain", ylab="Detected strain-injected strain", main=p("Strain sensitivity, equatorial sky band, spindown<", SpindownCutoff))
dev.off()

start_plot("iota_sensitivity.png")
plot(Data[goodBand & smallSpindown, 'fake_iota.3'], excessStrain[goodBand & smallSpindown], xlab="Iota", ylab="Detected strain-injected strain", main=p("Iota sensitivity, skyband 0, spindown<", SpindownCutoff))
dev.off()

start_plot("psi_sensitivity.png")
plot(Data[goodBand & smallSpindown, 'fake_psi.3'], excessStrain[goodBand & smallSpindown], xlab="Psi", ylab="Detected strain-injected strain", main=p("Psi sensitivity, skyband 0, spindown<", SpindownCutoff))
dev.off()

start_plot("phi_sensitivity.png")
plot(Data[goodBand & smallSpindown, 'fake_phi.3'], excessStrain[goodBand & smallSpindown], xlab="Phi", ylab="Detected strain-injected strain", main=p("Phi sensitivity, skyband 0, spindown<", SpindownCutoff))
dev.off()

start_plot("max_dx.png")
plot(injectedStrain[goodBand & smallSpindown], pmin(Data[goodBand & smallSpindown, 'max_dx.4'], SNRCap), xlab="Injected strain", ylab="Detection strength", main=p("skyband 0, spindown<", SpindownCutoff), log="x")
dev.off()

start_plot("max_dx_equator.png")
plot(injectedStrain[goodBand & smallSpindown & equatorialBand], pmin(Data[goodBand & smallSpindown & equatorialBand, 'max_dx.4'], 15), xlab="Injected strain", ylab="Detection strength", main=p("Equatorial sky band, spindown<", SpindownCutoff), log="x")
dev.off()

start_plot("f0_sensitivity.png")
plot(Data[goodBand & smallSpindown, 'fake_frequency.2'], excessStrain[goodBand & smallSpindown], xlab="f0", ylab="Detected strain-injected strain", main=p("f0 sensitivity, skyband 0, spindown<", SpindownCutoff))
dev.off()

F<- 2.0*pulsarUL<injectedStrain

X<-pi+(Data[,'fake_ra.3']-pi)*cos(Data[,'fake_dec.2'])
Y<-Data[,'fake_dec.2']
start_plot("sky_underestimate_distribution.png")
plot( X, Y, type='n')
points(X[goodBand & smallSpindown], Y[goodBand & smallSpindown], col="blue", pch="+")
points(X[goodBand & smallSpindown & F], Y[goodBand & smallSpindown & F], col="red", pch="+")
points(X[!goodBand & smallSpindown], Y[!goodBand & smallSpindown], col="black", pch="+")
dev.off()

X<- (Data[,'fake_frequency.2']*1800.0) %% 1.0
Y<-Data[,'fake_spindown.2']
start_plot("f_sp_underestimate_distribution.png")
plot( X, Y, type='n', log="y")
points(X[goodBand], Y[goodBand], col="blue", pch="+")
points(X[goodBand & F], Y[goodBand & F], col="red", pch="+")
points(X[!goodBand], Y[!goodBand], col="black", pch="+")
dev.off()

start_plot("optimization_time.png")
plot(Data[,'optimized_candidates_count.1'], Data[,'second_pass_time.4'])
dev.off()

if(!is.null(OptCand)) {
	start_plot("final_cand_dist.png")
	A<-M_MinDist
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'Dist'], main="Minimum distance to injection point", ylab="Radians", xlab="strain")
	dev.off()
	
	start_plot("final_cand_fdist.png")
	A<-M_MinFDist
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'FDist'], main="Minimum frequency difference", ylab="Frequency bins", xlab="strain", ylim=c(0, 20))
	dev.off()
	
	start_plot("final_cand_sdist.png")
	A<-M_MinSDist
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'SDist'], main="Minimum spindown difference", ylab="Spindown difference", xlab="strain")
	dev.off()
	
	start_plot("final_cand_power_cor.png")
	A<-M_MaxPowerCor
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'power_cor'], main="Maximum power correlation", ylab="Spindown difference", xlab="strain")
	dev.off()
	
	start_plot("final_cand_detection.png")
	A<-M_Detected
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'Detected'], main="(Dist<0.05)*1.1 +(FDist<4)*1.2+(SDist<1e-10)*1.4", ylab="Detection ?", xlab="strain")
	dev.off()
	
	start_plot("final_cand_detection2.png")
	A<-M_Detected2
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'Detected2'], main="(Dist<0.07)*1.1 +(FDist<5)*1.2+(SDist<2e-10)*1.4", ylab="Detection ?", xlab="strain")
	dev.off()
	
	start_plot("final_cand_detection3.png")
	A<-M_Detected3
	plot(A[abs(A[,'spindown'])<2e-10,'h0'], A[abs(A[,'spindown'])<2e-10,'Detected3'], main="(Merged[,'Dist']<0.07)*1.01 + (Merged[,'FDist']<5)*1.02 + (Merged[,'SDist']<2e-10)*1.04+(Merged[,'power_cor']>0.06)*1.08", ylab="Detection ?", xlab="strain")
	dev.off()
	}

# plot(Merged[MSmallSp,'h0'], Merged[MSmallSp,'snr'], log="xy")
# plot(Merged[MSmallSp,'h0'], Merged[MSmallSp,'Dist'], ylim=c(0, 0.1), log="x")
# plot(Merged[MSmallSp,'h0'], Merged[MSmallSp,'FDist'], ylim=c(0, 10), log="x")
# plot(Merged[MSmallSp,'h0'], Merged[MSmallSp,'SDist']/1800.0, ylim=c(0, 1e-9), log="x")
# plot(Merged[MSmallSp,'snr'], Merged[MSmallSp,'SDist']/1800.0, ylim=c(0, 1e-9), log="x")
# plot(Merged[MSmallSp,'snr'], Merged[MSmallSp,'FDist'], log="x", ylim=c(0,10))
# plot(Merged[MSmallSp,'snr'], Merged[MSmallSp,'Dist'], log="x", ylim=c(0,0.1))

#
# ROC curve
#
#ROCData<-data.frame(InjStrain=injectedStrain[goodBand & smallSpindown], SNR=Data[goodBand & smallSpindown, 'max_dx.4'], i=Data[goodBand & smallSpindown, 'i'])
ROCData<-data.frame(InjStrain=injectedStrain[goodBand & smallSpindown], SNR=Data[goodBand & smallSpindown, 'max_dx.4'], i=as.integer(Data[goodBand & smallSpindown, 'LogFile']))
ROCData<-ROCData[!is.na(ROCData[,'SNR']), ]
ROCData[, 'Detected']<- ROCData[,'SNR']>LargeSNR

ROCData<-ROCData[order(ROCData[,'InjStrain'], ROCData[,'i']), ]
ROCData[,'CumDet']<-cumsum(ROCData[,'Detected'])

First<-(1:(dim(ROCData)[1]))[!duplicated(ROCData[,c('InjStrain', 'i')])]
Last<-c(First[2:length(First)], dim(ROCData)[1])

ROCData2<-ROCData[First, c('InjStrain', 'i')]
ROCData2[,'Det']<-ROCData[Last, 'CumDet']-ROCData[First, 'CumDet']+ROCData[First, 'Detected']
#ROCData2<-aggregate(data.frame(Det=ROCData[,'Detected']), ROCData[,c("InjStrain", "i")], sum)
ROCData2[, 'Detected']<- ROCData2[,'Det']>0
ROCData2[, 'InjStrain']<- as.numeric(as.character(ROCData2[,'InjStrain']))
ROCData2<-ROCData2[order(ROCData2[,'InjStrain'], decreasing=TRUE),]
ROCData2[,'A']<-cumsum(ROCData2[,'Detected'])

N<-dim(ROCData2)[1]

MaxStrain<-max(ROCData2[,'InjStrain'])
MinStrain<-min(ROCData2[,'InjStrain'])
NGroups<-15
Group<- floor(NGroups*log10(ROCData2[,'InjStrain']/MinStrain)/(log10(MaxStrain/MinStrain)+1e-5))
First<-(1:N)[!duplicated(Group)]
Last<- ((N:1)[!duplicated(Group[N:1])])
Last<-Last[length(Last):1]

LevelLo<- exp(log(10)*(Group[First])*(log10(MaxStrain/MinStrain)+1e-5)/NGroups)*MinStrain
LevelHi<- exp(log(10)*(Group[First]+1)*(log10(MaxStrain/MinStrain)+1e-5)/NGroups)*MinStrain


write.table(ROCData2[,c('i', 'InjStrain', 'Detected')], "Injections.dat", sep="\t", row.names=FALSE, col.names=TRUE)

ROCData3<-data.frame(InjStrainHi=ROCData2[First, 'InjStrain'], 
	InjStrainLo=ROCData2[Last, 'InjStrain'], 
	Ratio=(ROCData2[Last, 'A']-ROCData2[First, 'A']+ROCData2[First, 'Detected'])/(Last-First+1),
	Count=Last-First+1)
ROCData3[,'Sd']<-sqrt( (ROCData3[,'Ratio']-ROCData3[,'Ratio']^2)*ROCData3[,'Count']/(ROCData3[,'Count']-1))
A<-sqrt(1/(ROCData3[,'Count']-1))
A<-1- 0.01^(1/ROCData3[,'Count'])
F<-ROCData3[,'Sd']<A
ROCData3[F,'Sd']<-A[F]

write.table(ROCData3, "ROCData.dat", sep="\t", row.names=FALSE, col.names=TRUE)

start_plot("ROC.png")
plot(LevelHi, ROCData3[,'Ratio'], xlab="Injected strain", ylab="Detection efficiency", main=p("Ratio of injections with spindown < ", SpindownCutoff, " with SNR>", LargeSNR), type="n")
grid()
#rect(ROCData3[,'InjStrainLo'], ROCData3[,'Ratio'], ROCData3[,'InjStrainHi'], ROCData3[,'Ratio'])
rect(LevelLo, ROCData3[,'Ratio']-ROCData3[,'Sd'], LevelHi, ROCData3[,'Ratio']+ROCData3[,'Sd'], density=-10, col="gray90", angle=45)
points(0.5*(LevelLo+LevelHi), ROCData3[,'Ratio'], type="l", col="blue")
dev.off()
