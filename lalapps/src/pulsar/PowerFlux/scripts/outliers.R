library("lattice")

A<-read.csv("outliers6.csv")

FPulsar0<- abs(A[,'H.frequency']-265.58)<0.1
FPulsar1<- abs(A[,'H.frequency']-849.08)<1
FPulsar2<- abs(A[,'H.frequency']-575.16)<0.1
FPulsar3<- abs(A[,'H.frequency']-108.86)<0.1
FPulsar6<- abs(A[,'H.frequency']-148.72)<2
FPulsar8<- abs(A[,'H.frequency']-193.72)<0.1
FPulsar9<- abs(A[,'H.frequency']-763.85)<0.1

F128 <- abs(A[,'H.frequency']-128.00)<0.1
F60 <-  ((A[,'H.frequency']+1.25) %% 60.0)<2.50

F64 <-  abs(A[,'H.frequency']-64.00)<1

LinesList<- sort(unique(c(64*1:5, 60*1:30, 56*1:5)))
FLines<- A[,'H.frequency']==0
for(line in LinesList)FLines<- FLines | (abs(A[,'H.frequency']-line)<0.1)

FPulsars<-FPulsar0 | FPulsar1 | FPulsar2 | FPulsar3 | FPulsar6 | FPulsar8 | FPulsar9

SNRMatch<- abs(log10(A[,'H.snr']/A[,'L.snr']))<0.3

FBin<- round((A[,'H.frequency']+A[,'L.frequency'])*5)/10.0
SBin<- round((A[,'H.spindown']+A[,'L.spindown'])*5e9)/10.0e9

Unknown<- sort(unique(FBin[!FPulsars & SNRMatch & !FLines]))

FUnknown<- A[,'H.frequency']==0
for(line in Unknown)FUnknown<- FUnknown | (abs(A[,'H.frequency']-line)<0.1)

A[,'FBin']<- FBin
A[,'FBinC']<-as.character(FBin)

A[,'SBin']<- SBin
A[,'SBinC']<-as.character(SBin)

A[,'CombinedSNR']<-sqrt(A[,'H.snr']^2+A[,'L.snr']^2)

A[,'MaxSNR']<- pmax(A[,'H.snr'], A[,'L.snr'])
A[,'MinSNR']<- pmin(A[,'H.snr'], A[,'L.snr'])

#xyplot(H.snr~L.snr|FBin, A[FUnknown,])

B<- aggregate(A[FUnknown, c('CombinedSNR', 'MaxSNR', 'MinSNR'), drop=FALSE], A[FUnknown, c("FBin", "SBin"),drop=FALSE], max, na.rm=TRUE)
B[,'FBin']<-as.numeric(as.character(B[,'FBin']))
B[,'SBin']<-as.numeric(as.character(B[,'SBin']))
B<-B[order(-B[,'MinSNR'], B[,'CombinedSNR'], B[,'FBin'], B[,'SBin']),,drop=FALSE]
B[,'IndexSNR']<-1:(dim(B)[1])
B[,'SortIndex']<-B[,'IndexSNR']+10^(floor(log10(dim(B)[1]))+1)*duplicated(B[,'FBin'])
B<-B[order(B[,'SortIndex']),,drop=FALSE]

dir.create("followup")
write.csv(B, "followup/followup.csv", row.names=FALSE)

for(i in 1:(dim(B)[1])) {
	config<-file(paste("followup/config.", i, sep=""), "w")
	cat(file=config, "first-bin", (B[i,"FBin"]-0.125)*1800.0, "\n")
	cat(file=config, "spindown-start", B[i,"SBin"], "\n")
	close(config)
	}


