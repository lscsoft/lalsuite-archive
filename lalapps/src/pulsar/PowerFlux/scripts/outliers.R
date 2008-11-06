library("lattice")

p<-function(...)paste(sep="", ...)

Prefix<-""
#Prefix<-"narrow_"

cat("Reading data\n")

A<-read.csv(p(Prefix, "outliers5.csv"))
F<-A[,'H.ifo_freq']<1
F[is.na(F)]<-FALSE
A[F, 'H.ifo_freq']<-NA
F<-A[,'L.ifo_freq']<1
F[is.na(F)]<-FALSE
A[F, 'L.ifo_freq']<-NA

OF1<-read.csv("orig_followup1.csv", header=TRUE)
OF1[,"Followup1"]<-p("#", 1:dim(OF1)[1])

cat("Processing\n")

# FPulsar0<- abs(A[,'H.frequency']-265.58)<0.1
# FPulsar1<- abs(A[,'H.frequency']-849.06)<1
# FPulsar2<- abs(A[,'H.frequency']-575.16)<0.5
# FPulsar3<- abs(A[,'H.frequency']-108.86)<0.1
# FPulsar4<- abs(A[,'H.frequency']-1401.55)<0.5
# FPulsar5<- abs(A[,'H.frequency']-52.80)<0.1
# FPulsar6<- abs(A[,'H.frequency']-148.29)<2
# FPulsar7<- abs(A[,'H.frequency']-1220.90)<0.2
# FPulsar8<- abs(A[,'H.frequency']-193.75)<0.1
# FPulsar9<- abs(A[,'H.frequency']-763.85)<0.1
# 
# F128 <- abs(A[,'H.frequency']-128.00)<0.1
# F60 <-  ((A[,'H.frequency']+1.25) %% 60.0)<2.50
# 
# F64 <-  abs(A[,'H.frequency']-64.00)<1
# 
# LinesList<- sort(unique(c(64*1:5, 60*1:30, 56*1:5)))
# FLines<- A[,'H.frequency']==0
# for(line in LinesList)FLines<- FLines | (abs(A[,'H.frequency']-line)<0.1)
# 
# FPulsars<-FPulsar0 | FPulsar1 | FPulsar2 | FPulsar3 | FPulsar6 | FPulsar8 | FPulsar9

# SNRMatch<- abs(log10(A[,'H.snr']/A[,'L.snr']))<0.3

FBin<- round((A[,'H.frequency']+A[,'L.frequency'])*5)/10.0
SBin<- round((A[,'H.spindown']+A[,'L.spindown'])*5e9)/10.0e9

#Unknown<- sort(unique(FBin[SNRMatch & !FLines]))
#Unknown<- sort(unique(FBin[!FLines]))

#FUnknown<- A[,'H.frequency']==0
#for(line in Unknown)FUnknown<- FUnknown | (abs(A[,'H.frequency']-line)<0.1)

#A[,"Comment"]<-""

#for(i in 0:9) {
#	A[get(paste("FPulsar", i, sep="")),"Comment"]<-paste("Inj pulsar", i)
#	}

A[,'FBin']<- FBin
#A[,'FBinC']<-as.character(FBin)

A[,'SBin']<- SBin
#A[,'SBinC']<-as.character(SBin)

A[,'CombinedSNR']<-sqrt(A[,'H.snr']^2+A[,'L.snr']^2)

A[,'MaxSNR']<- pmax(A[,'H.snr'], A[,'L.snr'])
A[,'MinSNR']<- pmin(A[,'H.snr'], A[,'L.snr'])

#xyplot(H.snr~L.snr|FBin, A[FUnknown,])

B<- aggregate(A[, c('CombinedSNR', 'MaxSNR', 'MinSNR'), drop=FALSE], A[, c("FBin", "SBin"),drop=FALSE], max, na.rm=TRUE)

C<-merge(A, B[,c("FBin", "SBin", "MinSNR", "MaxSNR"), drop=FALSE], by=c("FBin", "SBin"), suffixes=c("_orig", ""))

F<- (C[,"MinSNR_orig"]>=C[,"MinSNR"])
C<-C[F,,drop=FALSE]
rm(F)

cat("C has", dim(C)[1], "rows\n")
# C[,"FBin"]<-as.factor(C[,"FBin"])
# C[,"SBin"]<-as.factor(C[,"SBin"])
#rm(A)

# SNRLevel<-5.25
# F<-C[,"MinSNR"]>=SNRLevel
# C2<-C
# C<-C[F,,drop=FALSE]
# cat("C has", dim(C)[1], "rows after restricting to MinSNR>", SNRLevel, "\n")


H1only<- unique(C[,c("H.frequency", "H.spindown", "H.ra", "H.dec", "H.iota", "H.psi", "H.ifo_freq", "H.snr", "FBin", "SBin"), drop=FALSE])
H1Count<-aggregate(list(H.count=H1only[,"FBin"]), H1only[,c("FBin", "SBin"), drop=FALSE], length)

L1only<- unique(C[,c("L.frequency", "L.spindown", "L.ra", "L.dec", "L.iota", "L.psi", "L.ifo_freq", "L.snr", "FBin", "SBin"), drop=FALSE])
L1Count<-aggregate(list( L.count= L1only[,"FBin"]),  L1only[,c("FBin", "SBin"), drop=FALSE], length)

Index<-split(1:dim(C)[1], C[,c("FBin", "SBin"), drop=FALSE], drop=TRUE)

L<-list()

for(col in c("FBin", "SBin")) {
	cat("COPY:", col, "\n")
	X<-C[,col]
	L[[col]]<-unlist(lapply(Index, function(idx)X[[idx[[1]]]]))
	}

cat("Log10 H.snr/L.snr\n")
X<- abs(log10(C[,'H.snr']/C[,'L.snr']))
L[["LogSNRRatio"]]<-unlist(lapply(Index, function(idx)min(X[idx], na.rm=TRUE)))

for(col in c("H.frequency", "H.spindown", "H.ra", "H.dec", "H.iota", "H.psi", "H.ifo_freq", "H.snr",
	"L.frequency", "L.spindown", "L.ra", "L.dec", "L.iota", "L.psi", "L.ifo_freq", "L.snr", "MinSNR", "MaxSNR")) {
	cat("MEAN:", col, "\n")
	X<-C[,col]
	L[[col]]<-unlist(lapply(Index, function(idx)mean(X[idx], na.rm=TRUE)))
	}

for(col in c("H.frequency", "H.spindown", "H.ra", "H.dec", "H.iota", "H.psi", "H.ifo_freq", "H.snr",
	"L.frequency", "L.spindown", "L.ra", "L.dec", "L.iota", "L.psi", "L.ifo_freq", "L.snr")) {
	cat("SD:  ", col, "\n")
	X<-C[,col]
	L[[p(col, "_sd")]]<-unlist(lapply(Index, function(idx)sd(X[idx], na.rm=TRUE)))
	}



# re
# Means<- aggregate(unclass(C[,c("H.frequency", "H.spindown", "H.ra", "H.dec", "H.iota", "H.psi", "H.ifo_freq", "H.snr",
# 			"L.frequency", "L.spindown", "L.ra", "L.dec", "L.iota", "L.psi", "L.ifo_freq", "L.snr", "MinSNR", "MaxSNR"), drop=FALSE]), C[,c("FBin", "SBin"), drop=FALSE], mean, na.rm=TRUE)
# 
# Sds<- aggregate(C[,c("H.frequency", "H.spindown", "H.ra", "H.dec", "H.iota", "H.psi", "H.ifo_freq", "H.snr",
# 			"L.frequency", "L.spindown", "L.ra", "L.dec", "L.iota", "L.psi", "L.ifo_freq", "L.snr"), drop=FALSE], C[,c("FBin", "SBin"), drop=FALSE], sd, na.rm=TRUE)

Data0<-merge(H1Count, L1Count, by=c("FBin", "SBin"))
#Data1<- merge(Means, Sds, by=c("FBin", "SBin"), suffixes=c("", "_sd"))
Data<- merge(Data0, as.data.frame(L), by=c("FBin", "SBin"))


Data[,'FBin']<-as.numeric(as.character(Data[,'FBin']))
Data[,'SBin']<-as.numeric(as.character(Data[,'SBin']))
Data<-Data[order(-Data[,'MinSNR'], -Data[,'MaxSNR'], Data[,'FBin'], Data[,'SBin']),,drop=FALSE]
Data[,'IndexSNR']<-1:(dim(Data)[1])
Data[,'SortIndex']<-Data[,'IndexSNR']+10^(floor(log10(dim(Data)[1]))+1)*duplicated(Data[,'FBin'])
Data<-Data[order(Data[,'SortIndex']),,drop=FALSE]

Data[,"Index"]<- -1
Data[,"Comment"]<-""

Data<-merge(Data, Data[!duplicated(Data[,"FBin"]),c("FBin", "SortIndex")], by="FBin", all.x=TRUE, suffixes=c("", "_primary"))


Data[,'Frequency']<-0.5*(Data[,"H.frequency"]+Data[,"L.frequency"])
Data[,'Spindown']<-0.5*(Data[,"H.spindown"]+Data[,"L.spindown"])


line.list<-read.table(textConnection("Frequency Width Comment
 265.58 0.1 InjPulsar0
 849.06 1   InjPulsar1
 575.16 0.5 InjPulsar2
 108.86 0.1 InjPulsar3
1401.55 0.5 InjPulsar4
  52.80 0.1 InjPulsar5
 148.29 2   InjPulsar6
1220.90 0.2 InjPulsar7
 193.75 0.1 InjPulsar8
 763.85 0.1 InjPulsar9"), header=TRUE)

Data[,"Comment"]<-""
for(i in 1:dim(line.list)[1]) {
	F<- abs(Data[,"Frequency"]-line.list[i, "Frequency"])<line.list[i,"Width"]
	Data[F, "Comment"]<- as.character(line.list[i, "Comment"])
	}

F60 <-  ((Data[,'Frequency']+1.25) %% 60.0)<2.50
Data[F60, "Comment"]<- "60 Hz multiple"

for(f in (1:5)*64) {
	F<- abs(Data[,"Frequency"]-f)<0.1
	Data[F, "Comment"]<- "64 Hz multiple"
	}

F<-Data[,'LogSNRRatio']>0.3
Data[F,"Comment"]<-p(Data[F, "Comment"], " SNR ratio>=", round(10*10^Data[F, 'LogSNRRatio'])/10.0)

FUnknown<- Data[,"Comment"]==""

Data<-merge(Data, OF1[,c("FBin", "SBin", "Followup1"),drop=FALSE], by=c("FBin", "SBin"), all.x=TRUE)
Data<-merge(Data, OF1[!duplicated(OF1[,"FBin"]),c("FBin", "Followup1"),drop=FALSE], by=c("FBin"), all.x=TRUE, suffixes=c("", "_primary"))

for(col in c("Followup1", "Followup1_primary")) {
	Data[is.na(Data[,col]), col]<-""
	}

Data<-Data[order(Data[,"SortIndex"]),,drop=FALSE]
Data[,"Index"]<-1:dim(Data)[1]
#
# Rearrange columns
#
#CSd<-regexpr("_sd$", names(Data))>=0
Cols<-names(Data)
CParams<- (regexpr("^[H|L][.]", Cols)>=0)
NewCols<- c(Cols[!CParams], Cols[CParams])
Data<-Data[,NewCols, drop=FALSE]


cat(dim(Data)[1], "total bins,", sum(FUnknown), "unknown\n")

cat("Writing files\n")

ODIR<-p(Prefix, "followup")

dir.create(ODIR)
dir.create(p(ODIR, "/bins"))
dir.create(p(ODIR, "/H"))
dir.create(p(ODIR, "/L"))
write.csv(Data, p(ODIR, "/followup.csv"), row.names=FALSE)

for(i in 1:(dim(Data)[1])) {
	config<-file(p(ODIR, "/bins/config.", i), "w")
	cat(file=config, "first-bin", round((Data[i,"Frequency"]-0.125)*1800.0), "\n")
	cat(file=config, "spindown-start", Data[i,"Spindown"], "\n")
	close(config)

	for(prefix in c("H", "L")) {
		config<-file(p(ODIR, "/", prefix, "/config.", i), "w")
		cat(file=config, "first-bin", round((Data[i, p(prefix, ".frequency")]-0.125)*1800.0), "\n")
		cat(file=config, "spindown-start", Data[i, p(prefix, ".spindown")], "\n")
		cat(file=config, "focus-ra", Data[i, p(prefix, ".ra")], "\n")
		cat(file=config, "focus-dec", Data[i, p(prefix, ".dec")], "\n")
		close(config)
		}
	}


