require("lattice")
require("RMySQL")
require("parallel")

p<-function(...) {
	return(paste(..., sep=""))
	}


sft_length<-as.numeric(read.table(pipe("grep 'SFT coherence time' 0/*/powerflux.log"), header=FALSE, sep=":")[1,2])
if(is.na(sft_length))sft_length<-1800.0
cat("Using sft coherence length of", sft_length, "\n")

source("params.R")

cosdist<-function(ra1, dec1, ra2, dec2) (sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))
dist<-function(ra1, dec1, ra2, dec2) {
	a<-cosdist(ra1, dec1, ra2, dec2)
	return(ifelse(a>=1.0,  0.0 , ifelse(a<= -1, -pi, acos(a))))
	}

epsilon<-23.439281*pi/180
# x, y, z
ecliptic_pole<-c(0, -sin(epsilon), cos(epsilon))
ecliptic_x_axis<-c(1, 0, 0)
ecliptic_y_axis<-c(0, cos(epsilon), sin(epsilon))
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

con<-dbConnect(dbDriver("MySQL"), host=MYSQL_HOST, user=MYSQL_USER, password="", dbname=MYSQL_DB)

output.dir<-p("report/coinc_", Segment)
	
dir.create(output.dir)

LOG<-file(p(output.dir, "/inj_coinc.log"), open="w")
cat(file=LOG, "Using sft coherence length of", sft_length, "\n")

start_plot<-function(name, width=800, height=800) {
	png(p(output.dir, "/", name, ".png"), width=width, height=height, res=125)
	}

params<-read.table("params.txt", header=TRUE)

firstbin<-read.table("firstbin.txt", header=FALSE)
names(firstbin)<-c("tag", "v2", "firstbin")
firstbin[,"i"]<-as.integer(gsub(".*output/([0-9]*)/.*$", "\\1", firstbin[,"tag"]))
params<-merge(params, firstbin[,c("i", "firstbin"),drop=FALSE], by="i", all.x=TRUE)

side_cut<-read.table("side_cut.txt", header=FALSE)
names(side_cut)<-c("tag", "v2", "side_cut")
side_cut[,"i"]<-as.integer(gsub(".*output/([0-9]*)/.*$", "\\1", side_cut[,"tag"]))
params<-merge(params, side_cut[,c("i", "side_cut"),drop=FALSE], by="i", all.x=TRUE)

params[,"f0_bin"]<-params[,"f0"]*sft_length-params[,"firstbin"]-params[,"side_cut"]

# override standard dbQuery as it has issues
dbGetQuery<-function(con, query) {
	res<-dbSendQuery(con, query)
	L<-list()
	while(!dbHasCompleted(res)) {
		a<-fetch(res, FetchChunk)
                data<-a

                if(length(data)[1]<1 || dim(data)[1]<1)next

                L[[length(L)+1]]<-data
		}
	dbClearResult(res)
	if(length(L)<1) return(invisible(NULL))
		else
	if(length(L)>1) {
		return(do.call(rbind, L))
		} else return(L[[1]])
	}

# override standard dbQuery as it has issues
dbGetQueryCoinc<-function(con, query) {
        res<-dbSendQuery(con, query)
        L<-list()
        while(!dbHasCompleted(res)) {
                a<-fetch(res, FetchChunk)

                coincidences<-a

                if(length(coincidences)[1]<1 || dim(coincidences)[1]<1)next

		for(i in 1:(dim(coincidences)[2])) {
			if(class(coincidences[[i]])=="factor") {
				coincidences[[i]]<-as.character(coincidences[[i]])
				}
			}

                names(coincidences)<-p(names(coincidences), rep(c("", "_H1", "_L1"), each=dim(coincidences)[2]/3))
                #FDist<- ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_H1"], coincidences[,"dec_H1"])<(LocationTolerance(coincidences[,"frequency"])) & ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_L1"], coincidences[,"dec_L1"])<(LocationTolerance(coincidences[,"frequency"])) 
                #FDist[is.na(FDist)]<-FALSE
                #coincidences<-coincidences[FDist,,drop=FALSE]

                #cat("Found", dim(coincidences)[1],"coincidences out of", dim(a)[1], "rows fetched\n")

                #if(dim(coincidences)[1]<1)next

#               if(dim(coincidences)[1]>1000) {
#                       cat("truncating\n")
#                       coincidences<-coincidences[order(-coincidences$snr, -pmin(coincidences$snr_H1, coincidences$snr_L1), -pmax(coincidences$snr_H1, coincidences$snr_L1)),,drop=FALSE][1,,drop=FALSE]
#                       }

                L[[length(L)+1]]<-coincidences
                }
        dbClearResult(res)
        if(length(L)<1) return(invisible(NULL))
                else
        if(length(L)>1) {
                return(do.call(rbind, L))
                } else return(L[[1]])
        }

PrimaryMatches<-read.table(p(output.dir, "/primary_matches.csv"), header=TRUE, sep="\t")

PrimaryMatches<-merge(params, PrimaryMatches, by.x="i", by.y="Instance", all.x=TRUE, suffixes=c("_orig", ""))

Found<- abs(PrimaryMatches[,"f0"]-PrimaryMatches[,"f0_orig"])<FoundFrequencyTolerance & ecliptic_dist(PrimaryMatches[,"ra"], PrimaryMatches[,"dec"], PrimaryMatches[,"ra_orig"], PrimaryMatches[,"dec_orig"])<FoundLocationTolerance(PrimaryMatches[,"f0"]) & abs(PrimaryMatches[,"spindown"]-PrimaryMatches[,"spindown_orig"])<FoundSpindownTolerance

#Found<-  abs(PrimaryMatches[,"spindown"]-PrimaryMatches[,"spindown_orig"])<100*SpindownTolerance


Found[is.na(Found)]<-FALSE

names(PrimaryMatches)<-gsub("^phi$", "phi_orig", names(PrimaryMatches))
PrimaryMatches[,"line_id"]<-p(PrimaryMatches[,"i"], "_", 1:dim(PrimaryMatches)[1])
write.table(PrimaryMatches, p(output.dir, "/primary_matches_inj.csv"), sep="\t", col.names=TRUE, row.names=FALSE)

# Compute ranks
for(i in 1:(dim(params)[1])) {
	F<-PrimaryMatches[,"i"]==params[i,"i"]
	PrimaryMatches[F,"rank"]<-rank(-PrimaryMatches[F, "SNR.H1L1"])
	}

PrimaryMatches[is.na(PrimaryMatches[,"rank"]), "rank"]<-1

# Trim rank
F<-PrimaryMatches[,"rank"]<=MaxInjectionOutlierRank
Found<-Found[F]
PrimaryMatches<-PrimaryMatches[F,,drop=FALSE]
write.table(PrimaryMatches, p(output.dir, "/primary_matches_inj_trimmed.csv"), sep="\t", col.names=TRUE, row.names=FALSE)

EfficiencyTable<-data.frame(PrimaryMatches[,c("i", "h0", "spindown_orig", "f0_bin", "f0_orig", "ra_orig", "dec_orig", "psi_orig", "iota_orig", "SNR.H1L1", "SNR.L1", "SNR.H1", "f0", "ra", "dec", "spindown", "iota", "psi", "rank"),drop=FALSE], Found=as.integer(Found))
EfficiencyTable<-EfficiencyTable[order(-EfficiencyTable[,"Found"], -EfficiencyTable[,"SNR.H1L1"]),,drop=FALSE]
EfficiencyTable<-EfficiencyTable[!duplicated(EfficiencyTable[,"i"]),,drop=FALSE]

for(i in 1:dim(EfficiencyTable)[1]) {
	cat(i, "\n")
	EfficiencyTable[i,"H1L1.outliers"]<-dbGetQuery(con, p("SELECT COUNT(*) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_all' AND `label`=", EfficiencyTable[i, "i"]))[1,1]
	EfficiencyTable[i,"H1.outliers"]<-dbGetQuery(con, p("SELECT COUNT(*) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_LHO' AND `label`=",  EfficiencyTable[i, "i"]))[1,1]
	EfficiencyTable[i,"L1.outliers"]<-dbGetQuery(con, p("SELECT COUNT(*) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_LLO' AND `label`=",  EfficiencyTable[i, "i"]))[1,1]

	EfficiencyTable[i,"H1L1.max_snr"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_all' AND `label`=", EfficiencyTable[i, "i"]))[1,1]
	EfficiencyTable[i,"H1.max_snr"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_LHO' AND `label`=", EfficiencyTable[i, "i"]))[1,1]
	EfficiencyTable[i,"L1.max_snr"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_LLO' AND `label`=", EfficiencyTable[i, "i"]))[1,1]
	}

EfficiencyTable[,"data_present"]<-  !is.na(EfficiencyTable[,"H1L1.max_snr"]) & !is.na(EfficiencyTable[,"H1.max_snr"]) & !is.na(EfficiencyTable[,"L1.max_snr"])

EfficiencyTable[,"h0_adj"]<-EfficiencyTable[,"h0"]*sqrt(cos(EfficiencyTable[, "iota_orig"])^2+0.25)

EfficiencyTable[,"fdist"]<-abs(EfficiencyTable[,"f0"]-EfficiencyTable[,"f0_orig"])
EfficiencyTable[,"sdist"]<-abs(EfficiencyTable[,"spindown"]-EfficiencyTable[,"spindown_orig"])
EfficiencyTable[,"ecldist"]<- ecliptic_dist(EfficiencyTable[,"ra"], EfficiencyTable[,"dec"], EfficiencyTable[,"ra_orig"], EfficiencyTable[,"dec_orig"])


#print(xyplot(dec_orig~ra_orig|abs(dot2(ra_orig, dec_orig, ecliptic_pole))<0.25, EfficiencyTable[EfficiencyTable$Found<1 & EfficiencyTable$h0>1e-24,,drop=FALSE]))

start_plot("coincidence_rank_histrogram")
print(histogram(~rank, EfficiencyTable[EfficiencyTable[,"Found"]>0,,drop=FALSE], nint=200, xlim=c(0, 20)))
dev.off()


X<-EfficiencyTable[order(EfficiencyTable[,'h0']),,drop=FALSE]
#X<-X[X[,"f0_bin"]>=25 & X[,"f0_bin"]<=475,,drop=FALSE]
Summary_h0<-aggregate(X[,c("h0", "Found", "data_present"),drop=FALSE], list( floor((1:dim(X)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_h0")
print(xyplot(Found+data_present~log10(h0), Summary_h0))
dev.off()

EfficiencyTable<-EfficiencyTable[order(EfficiencyTable[,'h0_adj']),,drop=FALSE]
Summary_h0_adj<-aggregate(EfficiencyTable[,c("h0_adj", "Found", "data_present"),drop=FALSE], list( floor((1:dim(EfficiencyTable)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_h0_adj")
print(xyplot(Found+data_present~log10(h0_adj), Summary_h0_adj))
dev.off()

EfficiencyTable<-EfficiencyTable[order(EfficiencyTable[,'spindown_orig']),,drop=FALSE]
Summary_fdot<-aggregate(EfficiencyTable[,c("spindown_orig", "Found", "data_present"),drop=FALSE], list( floor((1:dim(EfficiencyTable)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_spindown")
print(xyplot(Found~spindown_orig, Summary_fdot))
dev.off()
#print(xyplot(Found~log10(-spindown_orig), Summary_fdot))

EfficiencyTable<-EfficiencyTable[order(EfficiencyTable[,'f0_orig']),,drop=FALSE]
Summary_f0<-aggregate(EfficiencyTable[,c("f0_orig", "Found", "data_present"),drop=FALSE], list( floor((1:dim(EfficiencyTable)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_f0")
print(xyplot(Found+data_present~f0_orig, Summary_f0))
dev.off()

EfficiencyTable<-EfficiencyTable[order(EfficiencyTable[,'H1L1.max_snr']),,drop=FALSE]
Summary_snr<-aggregate(EfficiencyTable[,c("H1L1.max_snr", "Found"),drop=FALSE], list( floor((1:dim(EfficiencyTable)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_H1L1_snr")
print(xyplot(Found~pmin(H1L1.max_snr,30), Summary_snr))
dev.off()

EfficiencyTable<-EfficiencyTable[order(EfficiencyTable[,'f0_bin']),,drop=FALSE]
Summary_f0_bin<-aggregate(EfficiencyTable[,c("f0_bin", "Found", "data_present"),drop=FALSE], list( floor((1:dim(EfficiencyTable)[1]-1) / 100)), mean)
start_plot("coincidence_efficiency_vs_f0_bin")
print(xyplot(Found+data_present~f0_bin, Summary_f0_bin))
dev.off()

#
# Plots relative to established upper limit value
#
UL<-read.table("upper_limits.csv", header=TRUE)

if(0) {
	RelEfficiencyTable<-EfficiencyTable

	for(i in 1:(dim(RelEfficiencyTable)[1])) {
		f<-RelEfficiencyTable[i,"f0_orig"]
		X<-UL[((UL$band+10.0/sft_length)<f) & ((UL$band+491.0/sft_length)>f),,drop=FALSE]
		RelEfficiencyTable[i, "UL"]<-max(X[,"ul"], na.rm=TRUE)
		RelEfficiencyTable[i, "non_gaussian"]<-max(X[,"non_gaussian"], na.rm=TRUE)
		}
	} else {

	RelEfficiencyTable<-merge(EfficiencyTable, UL, by.x="i", by.y="instance", suffix=c("", "_UL"))
	RelEfficiencyTable[,"UL"]<-RelEfficiencyTable[,"ul"]
	RelEfficiencyTable[,"non_gaussian"]<-(RelEfficiencyTable[,"max_m1_neg"]>0.44 | RelEfficiencyTable[,"min_m4"]<1.6)
	}

RelEfficiencyTable[,"rel_h0"]<-RelEfficiencyTable[,"h0"]/RelEfficiencyTable[,"UL"]

start_plot("coincidence_localization_f0")
print(xyplot(f0-f0_orig~rel_h0|ifelse(Found, "Found", "Not found"), RelEfficiencyTable, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
dev.off()

start_plot("coincidence_localization_f0_zoomed")
print(xyplot(f0-f0_orig~rel_h0|ifelse(Found, "Found", "Not found"), RelEfficiencyTable, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}, ylim=c(-3*FoundFrequencyTolerance, 3*FoundFrequencyTolerance)))
dev.off()

start_plot("coincidence_localization_spindown")
print(xyplot(spindown-spindown_orig~rel_h0|ifelse(Found, "Found", "Not found"), RelEfficiencyTable, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
dev.off()

start_plot("coincidence_localization_ecliptic_distance")
print(xyplot(I(ecldist*f0/400.0)~rel_h0|ifelse(Found, "Found", "Not found"), RelEfficiencyTable, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
dev.off()

X<-RelEfficiencyTable[order(RelEfficiencyTable[,'rel_h0']),,drop=FALSE]
X<-X[!is.na(X$rel_h0),,drop=FALSE]
#X<-X[X[,"f0_bin"]>=25 & X[,"f0_bin"]<=475,,drop=FALSE]
group<-floor((1:dim(X)[1]-1)/100)
if(sum(group==max(group))<100)group[group==max(group)]<-max(group)-1
Summary_rel_h0<-aggregate(X[,c("rel_h0", "Found", "data_present"),drop=FALSE], list(group), mean, na.rm=TRUE)
start_plot("coincidence_efficiency_vs_rel_h0")
print(xyplot(Found+data_present~rel_h0, Summary_rel_h0, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
dev.off()

Y<-X[X$non_gaussian>0,,drop=FALSE]
if(dim(Y)[1]>0) {
	group<-floor((1:dim(Y)[1]-1)/100)
	if(sum(group==max(group))<100)group[group==max(group)]<-max(group)-1
	Summary_rel_h0a<-aggregate(Y[,c("rel_h0", "Found", "data_present"),drop=FALSE], list(group), mean, na.rm=TRUE)
	png("coincidence_efficiency_vs_rel_h0_non_gaussian.png", width=600, height=600, res=125)
	print(xyplot(Found+data_present~rel_h0, Summary_rel_h0a, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
	dev.off()
	}

Y<-X[X$non_gaussian==0,,drop=FALSE]
if(dim(Y)[1]>0) {
	group<-floor((1:dim(Y)[1]-1) / 100)
	if(sum(group==max(group))<100)group[group==max(group)]<-max(group)-1
	Summary_rel_h0b<-aggregate(Y[,c("rel_h0", "Found", "data_present"),drop=FALSE], list(group), mean, na.rm=TRUE)
	png("coincidence_efficiency_vs_rel_h0_gaussian.png", width=600, height=600, res=125)
	print(xyplot(Found+data_present~rel_h0, Summary_rel_h0b, scales=list(x=list(log=TRUE, at=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10))), panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red"); panel.abline(v=0, col="red")}))
	dev.off()
	}
	
X<-RelEfficiencyTable[order(RelEfficiencyTable[,'f0_orig']),,drop=FALSE]
X<-X[!is.na(X$rel_h0),,drop=FALSE]
Y<-X[X$rel_h0>2,,drop=FALSE]
group<-floor((1:dim(Y)[1]-1)/100)
if(sum(group==max(group))<100)group[group==max(group)]<-max(group)-1
Summary_rel_h0a<-aggregate(Y[,c("f0_orig", "Found", "data_present"),drop=FALSE], list(group), mean, na.rm=TRUE)
start_plot("coincidence_efficiency_strong_vs_f0")
print(xyplot(Found+data_present~f0_orig, Summary_rel_h0a, panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red")}))
dev.off()

X<-RelEfficiencyTable[order(RelEfficiencyTable[,'spindown_orig']),,drop=FALSE]
X<-X[!is.na(X$rel_h0),,drop=FALSE]
Y<-X[X$rel_h0>2,,drop=FALSE]
group<-floor((1:dim(Y)[1]-1)/100)
if(sum(group==max(group))<100)group[group==max(group)]<-max(group)-1
Summary_rel_h0a<-aggregate(Y[,c("spindown_orig", "Found", "data_present"),drop=FALSE], list(group), mean, na.rm=TRUE)
start_plot("coincidence_efficiency_strong_vs_spindown")
print(xyplot(Found+data_present~spindown_orig, Summary_rel_h0a, panel=function(x,y,...){panel.xyplot(x,y,...); panel.abline(h=0.95, col="red")}))
dev.off()


X<-RelEfficiencyTable[order(RelEfficiencyTable[,'spindown_orig']),,drop=FALSE]
X<-X[!is.na(X$rel_h0) & X$rel_h0>2 & X$Found==0,,drop=FALSE]
start_plot("coincidence_missed_injections_f0_vs_spindown")
print(xyplot(spindown_orig~f0_orig, X))
dev.off()

start_plot("coincidence_missed_injections_ra_vs_dec")
print(xyplot(dec_orig~ra_orig, X))
dev.off()
