library("lattice")
p<-function(...)paste(sep="", ...)

source("params.R")

epsilon<-23.439281*pi/180
# x, y, z
ecliptic_pole<-c(0, -sin(epsilon), cos(epsilon))
dot2<-function(ra, dec, v) {
	return(v[1]*cos(ra)*cos(dec)+v[2]*sin(ra)*cos(dec)+v[3]*sin(dec))
	}

ecliptic_dist<-function(ra1, dec1, ra2, dec2) {
	x<-cos(ra1)*cos(dec1)-cos(ra2)*cos(dec2)
	y<-sin(ra1)*cos(dec1)-sin(ra2)*cos(dec2)
	z<-sin(dec1)-sin(dec2)
	
	a<- x*ecliptic_pole[1]+y*ecliptic_pole[2]+z*ecliptic_pole[3]
	return(sqrt((x-a*ecliptic_pole[1])^2+(y-a*ecliptic_pole[2])^2+(z-a*ecliptic_pole[3])^2))
	}

cluster_contiguous_values<-function(v, tolerance) {
	N<-length(v)
	P<-order(v)
	x<-v[P]
	FNewGroup<-c(TRUE, abs(x[1:(N-1)]-x[2:N])>tolerance)
	Group<-cumsum(FNewGroup)
	# return to given ordering
	Group[P]<-Group
	return(Group)
	}

A<-read.table("input_matches.csv", header=TRUE, stringsAsFactors=FALSE)
names(A)<-p(names(A), "_orig")
names(A)<-gsub("_orig_orig", "_inj", names(A))
names(A)<-gsub("h0_orig", "h0_inj", names(A))
cat("Loaded original input with", dim(A)[1], "rows\n")

col_names<-read.table(pipe("grep tag: 0/*/powerflux.log | head --lines=1"), header=FALSE, sep=" ")
B<-read.table("band_info.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(B)<-gsub(":", "", unlist(col_names[1,]))

data<-merge(A, B[B$kind=="snr" & B[,"set"]==p(Segment, "_all"),,drop=FALSE], by.x="line_id_orig", by.y="label", all.x=TRUE, suffixes=c("_orig", ""))
if(MinSNRfactor>0) {
	data<-merge(data, B[B$kind=="snr" & B[,"set"]==p(Segment, "_LHO"),,drop=FALSE], by.x="line_id_orig", by.y="label", all.x=TRUE, suffixes=c("", "_H1"))
	data<-merge(data, B[B$kind=="snr" & B[,"set"]==p(Segment, "_LLO"),,drop=FALSE], by.x="line_id_orig", by.y="label", all.x=TRUE, suffixes=c("", "_L1"))
	}
cat("Merged data has", dim(data)[1], "rows\n")

data[,"dist"]<-ecliptic_dist(data[,"ra_orig"], data[,"dec_orig"], data[,"ra"], data[,"dec"])
data[,"Comment"]<-data[,"Comment_orig"]

#
# Find masked entries
#
F<-data[,"frequency"]==0
F[is.na(F)]<-TRUE
cat("Found", sum(F), "masked entries\n")
data[F, "Comment"]<-p(data[F, "Comment"], " masked")
data[F, "frequency"]<-0
data[F, "spindown"]<- -100
data[F, "ra"]<- -10
data[F, "dec"]<- -10
data[F, "snr"]<- -10
data[F, "ul"]<- -1
data[F, "iota"]<- -10
data[F, "psi"]<- -10

data.same.segment<-data[data$set==data$set_orig,,drop=FALSE]

C<-data[order(-data[, "snr"]),,drop=FALSE]
data.max.snr<-C[!duplicated(C[,"line_id_orig"]),,drop=FALSE]

data.combined<- unique(rbind(data.same.segment, data.max.snr))
data.combined<-data.combined[order(data.combined[,"line_id_orig"], -data.combined[,"snr"]),,drop=FALSE]

cat("Combined data has", dim(data.combined)[1], "rows\n")

F<-(abs(data.combined[,"frequency"]-data.combined[,"frequency_orig"])<FrequencyTolerance) & (data.combined[,"frequency"]>0)
F[is.na(F)]<-FALSE

cat(sum(F), "have passed frequency cut\n")
data.combined[,"F0_cut"]<-F

F<-data.combined[,"snr"]>data.combined[,"snr_orig"]*LooseSNRfactor
F[is.na(F)]<-FALSE
data.combined[,"SNR_increase"]<-F
cat(sum(data.combined[,"SNR_increase"]), "outliers have passed SNR_increase cut\n")

if(MinSNRfactor<0) {
	data.combined[,"SNR_min"]<-TRUE
	cat("SNR_min cut disabled\n")
	} else {
	F<-(MinSNRfactor*data.combined[,"snr"]<data.combined[,"snr_H1"]) & (MinSNRfactor*data.combined[,"snr"]<data.combined[,"snr_L1"])
	F[is.na(F)]<-FALSE
	data.combined[,"SNR_min"]<-F
	cat(sum(data.combined[,"SNR_min"]), "outliers have passed SNR_min cut\n")
	}

keep.columns<-c(names(data.combined)[regexpr("(^dist|_orig|tag.)$", names(data.combined))<0], "line.f0_orig", "line.comment_orig", "min_gps_orig", "max_gps_orig", "nchunks_orig", "snr_orig", "i_orig", "line_id_orig", "tag")

# Drop columns we might not have
F<- ! keep.columns %in% names(data.combined)
if(sum(F)>0) {
	cat("*** missing columns:", p(keep.columns[F], collapse=" "), "\n")
	keep.columns<-keep.columns[!F]
	}

new_outliers<-unique(data.combined[,keep.columns,drop=FALSE])

#
# More messing with column names to make output suitable for input to the same script
#
names(new_outliers)<-gsub("_orig", "", names(new_outliers))
names(new_outliers)<-gsub("_inj", "_orig", names(new_outliers))
names(new_outliers)<-gsub("h0_orig", "h0", names(new_outliers))
new_outliers[,"line_id_orig"]<-new_outliers[,"line_id"]
new_outliers[,"line_id"]<-p(gsub(" ", "0", formatC(new_outliers[,"frequency"], digits=4, width=9, format="f")), "_", new_outliers[,"set"], "_", 1:(dim(new_outliers)[1]))

#
# Cluster outliers by 0.1 Hz groups which would all be using the same background plot
#
new_outliers[,"background_group"]<-cluster_contiguous_values(new_outliers[,"frequency"], 0.05)

# Add columns to fit Joe's script
L<-list(f0="frequency", f1="spindown")
for(col in names(L))
	new_outliers[,col]<-new_outliers[,L[[col]]]

write.table(new_outliers, "followup_matches.raw.csv", sep="\t", col.names=TRUE, row.names=FALSE)
reduced_outliers<-new_outliers[new_outliers[,"SNR_increase"] & new_outliers[,"SNR_min"] & new_outliers[,"F0_cut"],,drop=FALSE]
cat(dim(reduced_outliers)[1], "reduced outliers\n")

#
# Identify largest SNR in each background group
#
reduced_outliers<-reduced_outliers[order(reduced_outliers[,"background_group"],  reduced_outliers[,"i"], reduced_outliers[,"snr"]),,drop=FALSE]
reduced_outliers[,"primary"]<- !duplicated(reduced_outliers[,c("background_group", "i")])

reduced_outliers<-reduced_outliers[order(-reduced_outliers[,"primary"], reduced_outliers[,"frequency"]),,drop=FALSE]

cat(sum(reduced_outliers[,"primary"]), "primary outliers\n")

write.table(reduced_outliers, "followup_matches.csv", sep="\t", col.names=TRUE, row.names=FALSE)
