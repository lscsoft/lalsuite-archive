require("lattice")
require("RMySQL")

p<-function(...) {
	return(paste(..., sep=""))
	}

universal_statistics<-FALSE

sft_length<-as.numeric(read.table(pipe("grep 'SFT coherence time' 0/*/powerflux.log"), header=FALSE, sep=":")[1,2])
if(is.na(sft_length))sft_length<-1800.0
cat("Using sft coherence length of", sft_length, "\n")

source("params.R")

cosdist<-function(ra1, dec1, ra2, dec2) (sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))

con<-dbConnect(dbDriver("MySQL"), host=MYSQL_HOST, user=MYSQL_USER, password="", dbname=MYSQL_DB)

Cputime<-read.table(p("runtime.txt"), header=FALSE)
RunningInstances<-list()
try({
	Running<-read.table(p("running.txt"), header=FALSE)
	RunningInstances<- as.integer(as.character(gsub(".*err.", "", Running[,3])))
	})
VetoedInstances<-(1:2)[0]
try({
	Vetoed<-read.table(p("vetoed.txt"), header=FALSE)
	VetoedInstances<- as.integer(as.character(gsub(".*err.(.*):.*", "\\1", Vetoed[,1])))
	})

output.dir<-p("report/ul_", Segment)

dir.create(output.dir, recursive=TRUE)

LOG<-file(p(output.dir, "/band_info.log"), open="w")

names(Cputime)<-c("tag", "x", "seconds")
Cputime[,"Band"]<- as.numeric(gsub("^.*-(.*)/powerflux.*$", "\\1", Cputime[,"tag"]))
Cputime[,"Hours"]<-Cputime[,'seconds']/3600

print(xyplot(Hours~Band, Cputime, auto.key=list(columns=2, lines=TRUE, points=FALSE)))

#runtime<-function(band)(1.5*(band/100)^2)

runtime<-function(band)(12.35*(band/100)^2-4.1*band/100+25)

start_plot<-function(name, width=1024, height=1024, res=125, ...) {
	png(p(output.dir, "/", name, ".png"), width=width, height=height, res=res, ...)
	}

start_plot("timing")
print(xyplot(Hours~Band, Cputime))
dev.off()

grid_points<-read.table("grid_points.txt", header=FALSE)
names(grid_points)<-c("tag", "skyband", "skyband_name", "count")
grid_points[,"i"]<-as.numeric(gsub("/.*", "", gsub("output/", "", grid_points[,"tag"])))

total_grid_points<-aggregate(grid_points[,"count", drop=FALSE], grid_points[,"i", drop=FALSE], sum)

spindown_counts<-read.table("spindown_count.txt", header=FALSE)
names(spindown_counts)<-c("tag", "x", "spindown_count")
spindown_counts[,"i"]<-as.numeric(gsub("/.*", "", gsub("output/", "", spindown_counts[,"tag"])))

band_info_counts<-dbGetQuery(con, p("SELECT SUM(valid_count) as valid_count, SUM(masked_count) as masked_count, SUM(template_count) as template_count, label, `set`, kind, COUNT(*) as entry_count FROM `", BandInfoTableName, "` WHERE kind='ul' AND `set`='", Segment, "_all' GROUP BY label"))

nshifts<-read.table(pipe("grep nfshift: 0/*/powerflux.log"), header=FALSE)[1,2]

counts_summary<-merge(total_grid_points, spindown_counts, by="i", all.x=TRUE)
counts_summary<-merge(counts_summary, band_info_counts, by.x="i", by.y="label", all=TRUE)

counts_summary[,"scheduled_units"]<-counts_summary[,"spindown_count"]*counts_summary[,"count"]*nshifts
counts_summary[,"completed_units"]<-counts_summary[,"valid_count"]+counts_summary[,"masked_count"]
cat(file=LOG, "Total work units scheduled", sprintf("%.0f", sum(as.numeric(counts_summary[,"scheduled_units"]))), "\n")
cat(file=LOG, "Total work units completed", sprintf("%.0f", sum(as.numeric(counts_summary[,"completed_units"]))), "\n")

RedoInstances<-counts_summary[,"scheduled_units"]!=counts_summary[,"completed_units"]
RedoInstances[is.na(RedoInstances)]<-TRUE
RedoInstances<-counts_summary[RedoInstances, "i"]

CompletedInstances<- sort(as.integer(as.character(gsub(".*output/([^/]*)(/.*)?/powerflux.log.*", "\\1", Cputime[,"tag"])))) 
cat(file=LOG, "Last instance completed:", max(CompletedInstances), "\n")
cat(file=LOG, "Found", length(RunningInstances), "running instances\n")
cat(file=LOG, "Found", length(VetoedInstances), "vetoed instances\n")
IncompleteInstances<-unique(c(RedoInstances, sort(setdiff(1:max(CompletedInstances), CompletedInstances))))
cat(file=LOG, "Found", length(IncompleteInstances), "incomplete instances:", p(IncompleteInstances, collapse=" "), "\n")
MissingInstances<-sort(unique(c(RedoInstances, sort(setdiff(1:max(CompletedInstances), c(CompletedInstances, RunningInstances, VetoedInstances))))))
cat(file=LOG, "Found", length(MissingInstances), "missing instances:", p(MissingInstances, collapse=" "), "\n")

cat(file=LOG, "-------------------- missing.dag ----------------------\n")
cat(file=LOG, p("JOB A", MissingInstances, " condor\nVARS A", MissingInstances, " PID=\"", MissingInstances, "\"\n", collapse=""))
cat(file=LOG, "-------------------------------------------------------\n")


#runtime<-function(band)(0.5+1.5*(band/100)^2)
# 675144 hours for single spindown (neglecting startup costs)
# sames as 70 days on 400 nodes

#SkyBands<-dbGetQuery(con, p("SELECT DISTINCT skyband_name FROM `", BandInfoTableName, "`"))[,1]

UpperLimits_all<-dbGetQuery(con, p("SELECT first_bin, skyband_name, kind, MAX(ul) as ul, MIN(min_m4) as min_m4, MAX(max_m4) as max_m4, MAX(max_m1_neg) as max_m1_neg FROM `", BandInfoTableName, "` WHERE  `set`='", Segment,  "_all' AND kind IN ('ul', 'circ') GROUP BY first_bin, skyband_name, kind"))

UpperLimits_LHO<-dbGetQuery(con, p("SELECT first_bin, skyband_name, kind, MAX(ul) as ul, MIN(min_m4) as min_m4, MAX(max_m4) as max_m4, MAX(max_m1_neg) as max_m1_neg FROM `", BandInfoTableName, "` WHERE  `set`='", Segment,  "_LHO' AND kind IN ('ul', 'circ') GROUP BY first_bin, skyband_name, kind"))

UpperLimits_LLO<-dbGetQuery(con, p("SELECT first_bin, skyband_name, kind, MAX(ul) as ul, MIN(min_m4) as min_m4, MAX(max_m4) as max_m4, MAX(max_m1_neg) as max_m1_neg FROM `", BandInfoTableName, "` WHERE  `set`='", Segment,  "_LLO' AND kind IN ('ul', 'circ') GROUP BY first_bin, skyband_name, kind"))

UpperLimits<-merge(UpperLimits_all, UpperLimits_LHO, by=c("first_bin", "skyband_name", "kind"), suffixes=c("", "_LHO"), all=TRUE)
UpperLimits<-merge(UpperLimits, UpperLimits_LLO, by=c("first_bin", "skyband_name", "kind"), suffixes=c("_all", "_LLO"), all=TRUE)

# In SQL MAX() of no entries reports 0. Patch up !
UpperLimits[,"ul_all"]<-ifelse(UpperLimits[,"ul_all"]<=0.0, NA, UpperLimits[,"ul_all"])
UpperLimits[,"ul_LHO"]<-ifelse(UpperLimits[,"ul_LHO"]<=0.0, NA, UpperLimits[,"ul_LHO"])
UpperLimits[,"ul_LLO"]<-ifelse(UpperLimits[,"ul_LLO"]<=0.0, NA, UpperLimits[,"ul_LLO"])
if(!universal_statistics) {
	UpperLimits[,"NonGaussian_all"]<-(UpperLimits[,"max_m1_neg_all"]>0.44 | UpperLimits[,"min_m4_all"]<1.95)
	UpperLimits[,"NonGaussian_LHO"]<-(UpperLimits[,"max_m1_neg_LHO"]>0.44 | UpperLimits[,"min_m4_LHO"]<1.95)
	UpperLimits[,"NonGaussian_LLO"]<-(UpperLimits[,"max_m1_neg_LLO"]>0.44 | UpperLimits[,"min_m4_LLO"]<1.95)
	} else {
	UpperLimits[,"NonGaussian_all"]<-FALSE
	UpperLimits[,"NonGaussian_LHO"]<-FALSE
	UpperLimits[,"NonGaussian_LLO"]<-FALSE
	}

Contaminated<-UpperLimits[,"NonGaussian_all"] & UpperLimits[,"NonGaussian_LLO"] & UpperLimits[,"NonGaussian_LHO"]

UpperLimits[,"ul"]<-pmin(ifelse(UpperLimits[,"NonGaussian_all"], NA, UpperLimits[,"ul_all"]), ifelse(UpperLimits[,"NonGaussian_LHO"], NA, UpperLimits[,"ul_LHO"]), ifelse(UpperLimits[,"NonGaussian_LLO"], NA, UpperLimits[,"ul_LLO"]), na.rm=TRUE)

sft_length<-as.numeric(read.table(pipe("grep 'SFT coherence time' 0/*/powerflux.log"), header=FALSE, sep=":")[1,2])
if(is.na(sft_length))sft_length<-1800.0
cat("Using sft coherence length of", sft_length, "\n")

UpperLimits[,"band"]<-UpperLimits[,"first_bin"]/sft_length
UpperLimits[,"F60Hz"]<- ((UpperLimits[,"band"]+1.25) %% 60.0)<2.5

#
# Remove very high upper limits caused by instrumental lines - we are likely not correct there anyway
#
F<-UpperLimits[,"ul"]>MaxUL
F[is.na(F)]<-FALSE
if(sum(F)>1) {
	cat(file=LOG, "Found", sum(F), "very high upper limits in bands", p(sort(unique(UpperLimits[F, "first_bin"])), collapse=" "), "\n")
	UpperLimits[F,"ul"]<-NA
	}

F<-is.na(UpperLimits[,"ul"]) & UpperLimits[,"skyband_name"]!="Lines" & !UpperLimits[,"F60Hz"]
if(sum(F)>1) {
	cat(file=LOG, "Found", sum(F), "bands without an upper limit (bins):", p(sort(unique(UpperLimits[F, "first_bin"])), collapse=" "), "\n")
	cat(file=LOG, "Found", sum(F), "bands without an upper limit (freq):", p(sort(unique(UpperLimits[F, "band"])), collapse=" "), "\n")
	}

write.table(UpperLimits, p(output.dir, "/upper_limits.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

ul_plot<-function(data_orig, title="UL", correction=ULfactor, ...) {
data<-data_orig
data[,"ul"]<-data[,"ul"]*correction

data<-data[order(-data[,"ul"]),,drop=FALSE]
F<- !duplicated(data[,"band"])

data<-data[F,,drop=FALSE]

data<-data[order(data[,"band"]),,drop=FALSE]

plot((data$band), log10(data$ul), main=title, ylab="log10(h0)", xlab="Band", type="n", ...)

FUL60Hz<- abs((data$band+1.25) %% 60)<2.50
points((data$band[FUL60Hz]), log10(data$ul[FUL60Hz]), col="cyan")

F<- !FUL60Hz
points((data$band[F]), log10(data$ul[F]), col="green4", type="l")
}

start_plot("ul_allsky")
ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"]!="Lines",], title="All-sky upper limits")
dev.off()

start_plot("ul_allsky_125_200")
ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"]!="Lines",], title="All-sky upper limits", xlim=c(125, 200), ylim=c(-25, -23))
dev.off()

#start_plot("ul_equator")
#ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"] %in% c("equator"),], title="Equatorial upper limits")
#dev.off()

#start_plot("ul_mid")
#ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"] %in% c("midnorth", "midsouth"),], title="Mid-sky upper limits")
#dev.off()

#start_plot("ul_polar")
#ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"] %in% c("north", "south"),], title="Polar upper limits")
#dev.off()

#start_plot("ul_circ")
#ul_plot(UpperLimits[UpperLimits[,"kind"]=="circ" & UpperLimits[,"skyband_name"]!="Lines",], title="Best case all-sky upper limits", correction=1.0)
#dev.off()

start_plot("ul_lines")
try({
ul_plot(UpperLimits[UpperLimits[,"kind"]=="ul" & UpperLimits[,"skyband_name"] %in% c("Lines"),], title="Artifact affected upper limits")})
dev.off()

aggregate_plot<-function(table=BandInfoTableName, quantity="SUM(template_count)", xlab="Band", where.clause=p("kind='ul' "), ...) {
	Data<-dbGetQuery(con, p("SELECT first_bin, skyband_name, ", quantity, " as agg_quantity FROM `", table, "` WHERE  `set`='", Segment, "_all' AND (", where.clause, ") GROUP BY first_bin, skyband_name"))
	X<-reshape(Data, idvar="first_bin", timevar="skyband_name", direction="wide")
	names(X)<-gsub("agg_quantity.", "",names(X))
	X[,"band"]<-X[,"first_bin"]/sft_length
	X<-X[order(X$band),,drop=FALSE]
	
	print(xyplot(as.formula(p(p(sort(unique(Data[,"skyband_name"])), collapse="+"), "~band")), X, auto.key=list(columns=6), type="l", xlab=xlab, ...))
	}

start_plot("template_count")
aggregate_plot(ylab="Template count")
dev.off()

start_plot("valid_count")
aggregate_plot(quantity="SUM(valid_count)", ylab="Valid count")
dev.off()

start_plot("masked_count")
aggregate_plot(quantity="LOG(SUM(masked_count)+1)", ylab="log10(Masked count+1)")
dev.off()

start_plot("masked_ratio")
aggregate_plot(quantity="SUM(masked_count)/(SUM(masked_count)+SUM(valid_count))", ylab="Masked ratio")
dev.off()

start_plot("snr")
aggregate_plot(quantity="IF(MAX(snr)<100, MAX(snr), 100)", ylab="SNR (clamped)", where.clause=p("kind='snr'"))
dev.off()

Data<-dbGetQuery(con, p("SELECT first_bin, skyband_name, spindown, snr, masked_count FROM `", BandInfoTableName, "` WHERE  `set`='", Segment, "_all' AND kind='snr' AND skyband_name!='Lines'"))

BinData<-function(x, group=c(0, 2, 2.5, 3, 5, 6, 7,10,15,20,50)) {
	y<-x
	F0<-FALSE
	for(g in sort(group,decreasing=TRUE)) {
		F<- (y>=g) & ! F0
		y[F]<-g
		F0<-F0 | F
		}
	y[!F0]<-g
	return(y)
	}
Data[,"group"]<-BinData(Data[,"snr"])
Data[,"frequency"]<-Data[,"first_bin"]/sft_length

start_plot("snr_map")
print(xyplot(spindown~frequency, Data, groups=Data[,"group"], pch="o", auto.key=list(columns=5), par.settings=list(superpose.symbol=list(col=rev(rainbow(11))))))
dev.off()

Data_ul<-dbGetQuery(con, p("SELECT first_bin, skyband_name, spindown, ul, masked_count FROM `", BandInfoTableName, "` WHERE  `set`='", Segment, "_all' AND kind='ul' AND skyband_name!='Lines'"))

Data_ul[,"group"]<-BinData(Data_ul[,"ul"]/cos(pi/16), c(0, 9e-25, 1e-24, 2e-24, 3e-24, 4e-24, 5e-24, 1e-23, 2e-23, 3e-23) )
Data_ul[,"frequency"]<-Data_ul[,"first_bin"]/sft_length

start_plot("equator_ul_map")
print(xyplot(spindown~frequency, Data_ul, groups=Data_ul[,"group"], pch="o", auto.key=list(columns=5), par.settings=list(superpose.symbol=list(col=rev(rainbow(length(unique(Data_ul[,"group"]))))))))
dev.off()

re


start_plot("gauss_test_max_m1_neg")
aggregate_plot(quantity="MAX(IF(template_count=0, NULL, max_m1_neg))", ylab="max_m1_neg")
dev.off()

start_plot("gauss_test_min_m1_neg")
aggregate_plot(quantity="MIN(IF(template_count=0, NULL, min_m1_neg))", ylab="min_m1_neg")
dev.off()

start_plot("gauss_test_max_m3_neg")
aggregate_plot(quantity="MAX(IF(template_count=0, NULL, max_m3_neg))", ylab="max_m3_neg")
dev.off()

start_plot("gauss_test_min_m3_neg")
aggregate_plot(quantity="MIN(IF(template_count=0, NULL, min_m3_neg))", ylab="min_m3_neg")
dev.off()

start_plot("gauss_test_max_m4")
aggregate_plot(quantity="MAX(IF(template_count=0, NULL, log10(max_m4)))", ylab="log10(max_m4)")
dev.off()

start_plot("gauss_test_min_m4")
aggregate_plot(quantity="MIN(IF(template_count=0, NULL, min_m4))", ylab="min_m4")
dev.off()

ul_plot<-function(table=BandInfoTableName, where.clause="skyband_name!='Lines' AND kind='ul'", title="UL", correction=ULfactor) {
data<- dbGetQuery(con, p("SELECT first_bin, MAX(ul) as ul, MIN(min_m4) as min_m4, MAX(max_m4) as max_m4, MAX(max_m1_neg) as max_m1_neg FROM `", table, "` WHERE  `set`='", Dataset, "' AND (", where.clause, ") GROUP BY first_bin"))
data[,"ul"]<-data[,"ul"]*correction
plot((data$first_bin/sft_length), log10(data$ul), main=title, ylab="log10(h0)", xlab="Band", type="n")


FUL60Hz<- abs((data$first_bin/sft_length+1.25) %% 60)<2.50
points((data$first_bin[FUL60Hz]/sft_length), log10(data$ul[FUL60Hz]), col="cyan")
FUL<- !FUL60Hz

FNonGaussian<-(data[,"max_m1_neg"]>0.44 | data[,"min_m4"]<1.95)
F<-FUL & FNonGaussian
# F<- FUL & (band_info_full_ul$first_bin %in% SSmallKSValueBins)
points((data$first_bin[F]/sft_length), log10(data$ul[F]), col="blue", pch=23)
FUL<-FUL & !FNonGaussian

FDisturbedSpectrum<-data[,"max_m4"]>10
F<-FUL & FDisturbedSpectrum
# F<- FUL & (band_info_full_ul$first_bin %in% SSmallKSValueBins)
points((data$first_bin[F]/sft_length), log10(data$ul[F]), col="red", pch=23)
FUL<-FUL & !FDisturbedSpectrum

# 
# F<- FUL & (band_info_full_ul$first_bin %in% SLargeKSValueBins)
# points((band_info_full_ul$first_bin[F]/1800), log10(band_info_full_ul$ul[F]), col="blue", pch=23)
# FUL<-FUL & !(band_info_full_ul$first_bin %in% SLargeKSValueBins)

F<-FUL
points((data$first_bin[F]/sft_length), log10(data$ul[F]), col="green4", type="l")
}

snr_plot<-function(table=BandInfoTableName, where.clause="skyband_name!='Lines' AND kind='snr'", title="SNR", clamp=20, snr.high=7) {
data<- dbGetQuery(con, p("SELECT first_bin, MAX(snr) as snr, MIN(min_m4) as min_m4, MAX(max_m1_neg) as max_m1_neg FROM `", table, "` WHERE  `set`='", Dataset, "' AND (", where.clause, ") GROUP BY first_bin"))
data[,"snr"]<- pmin(data[,"snr"], clamp)
plot((data$first_bin/sft_length), data$snr, main=title, ylab="SNR", xlab="Band", type="n")


FUL60Hz<- abs((data$first_bin/sft_length+1.25) %% 60)<2.50
points((data$first_bin[FUL60Hz]/sft_length), data$snr[FUL60Hz], col="cyan")
FUL<- !FUL60Hz

FNonGaussian<-(data[,"max_m1_neg"]>0.44 | data[,"min_m4"]<1.95) 
F<-FUL & FNonGaussian
# F<- FUL & (band_info_full_ul$first_bin %in% SSmallKSValueBins)
points((data$first_bin[F]/sft_length), data$snr[F], col="blue", pch=23)
FUL<-FUL & !FNonGaussian
# 
# F<- FUL & (band_info_full_ul$first_bin %in% SLargeKSValueBins)
# points((band_info_full_ul$first_bin[F]/1800), log10(band_info_full_ul$ul[F]), col="blue", pch=23)
# FUL<-FUL & !(band_info_full_ul$first_bin %in% SLargeKSValueBins)

FHigh<-data[,"snr"]>snr.high
F<-FUL & FHigh
# F<- FUL & (band_info_full_ul$first_bin %in% SSmallKSValueBins)
points((data$first_bin[F]/sft_length), data$snr[F], col="red", pch="+")
FUL<-FUL & !FHigh

F<-FUL
points((data$first_bin[F]/sft_length), data$snr[F], col="green4", type="l")
}

value_map<-function(value, table=BandInfoTableName, kind="ul", where.clause=p("skyband_name!='Lines' AND kind='", kind, "'"), title=value, quantile.breaks=NULL, breaks=NULL, z.max.clamp=NULL, z.min.clamp=NULL, sql.statement=p("SELECT first_bin, spindown, ", value, " as Z FROM `", table, "` WHERE `set`='", Dataset, "' AND (", where.clause, ") ")) {
data<- dbGetQuery(con, sql.statement)

#browser()

if(!is.null(z.max.clamp))data[,"Z"]<-pmin(data[,"Z"], z.max.clamp)
if(!is.null(z.min.clamp))data[,"Z"]<-pmax(data[,"Z"], z.min.clamp)

try({
plot( (data$first_bin/sft_length), data$spindown, type="n", main=title, xlab="Band", ylab="Spindown")
})

if(is.null(breaks)) {
	if(!is.null(quantile.breaks))breaks<-quantile(data$Z, probs=quantile.breaks)
		else breaks<- (0:20)*(max(data$Z, na.rm=TRUE)-min(data$Z, na.rm=TRUE))/20.0+min(data$Z, na.rm=TRUE)
	}

breaks<-breaks[!is.na(breaks)]

Group<-findInterval(data$Z, breaks, all.inside=TRUE)
Colors<-topo.colors(length(breaks)-1)

for(i in 1:(length(breaks)-1)) {
	F<-Group==i
	F[is.na(F)]<-FALSE
	points( (data$first_bin/sft_length)[F], data$spindown[F], pch=19, col=Colors[[i]])
	}

}

ul_plot()

start_plot("ul")
ul_plot()
dev.off()

start_plot("circ", correction=1.0, title="Circular polarization UL")
ul_plot(where.clause="skyband_name!='Lines' AND kind='circ'")
dev.off()

start_plot("circ_polar", correction=1.0, title="Circular polarization UL")
ul_plot(where.clause="skyband_name IN ('north', 'south') AND kind='circ'")
dev.off()

start_plot("ul_equator")
ul_plot(where.clause="skyband_name='equator' AND kind='ul'", title="Equatorial band UL")
dev.off()

start_plot("ul_lines")
ul_plot(where.clause="skyband_name='Lines' AND kind='ul'", title="Lines sky band UL")
dev.off()

start_plot("ul_midband")
ul_plot(where.clause="skyband_name IN ('midsouth', 'midnorth') AND kind='ul'", title="Mid sky bands UL")
dev.off()

start_plot("ul_polar")
ul_plot(where.clause="skyband_name IN ('south', 'north') AND kind='ul'", title="Polar UL")
dev.off()

snr_plot()

start_plot("snr")
snr_plot()
dev.off()

start_plot("snr_equator")
snr_plot(where.clause=p("skyband_name='equator' AND kind='snr'"), title="Equatorial band SNR")
dev.off()

start_plot("snr_lines")
snr_plot(where.clause="skyband_name='Lines' AND kind='snr'", title="Lines sky band SNR")
dev.off()

start_plot("snr_midband")
snr_plot(where.clause="skyband_name IN ('midsouth', 'midnorth') AND kind='snr'", title="Mid sky bands SNR")
dev.off()

start_plot("snr_polar")
snr_plot(where.clause="skyband_name IN ('south', 'north') AND kind='snr'", title="Polar SNR")
dev.off()

value_map("snr", z.max.clamp=20, kind="snr")

start_plot("snr_map")
value_map("snr", z.max.clamp=20, kind="snr")
dev.off()

start_plot("circ_ul_map")
value_map("log10(ul)", kind="circ")
dev.off()

start_plot("ul_map")
value_map("log10(ul)", kind="ul")
dev.off()

start_plot("max_m1_neg_map")
value_map("max_m1_neg", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("min_m1_neg_map")
value_map("-min_m1_neg", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("max_m3_neg_map")
value_map("max_m3_neg", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("min_m3_neg_map")
value_map("-min_m3_neg", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("max_m4_map")
value_map("max_m4", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("min_m4_map")
value_map("-min_m4", kind="ul", quantile.breaks=c(0, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()

start_plot("max_weight_loss_fraction_map")
value_map("max_weight_loss_fraction", kind="ul", breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1))
dev.off()

start_plot("masked_count_map")
value_map("masked_count", kind="ul", breaks=c(0, 1, 10, 20, 50, 100, 1e10))
dev.off()

start_plot("masked_ratio_map")
value_map("masked_count/(masked_count+valid_count)", breaks=c(0, 0.1, 0.5, 0.8, 0.9, 1), kind="ul")
dev.off()

high_bands<- sort(dbGetQuery(con, p("SELECT DISTINCT first_bin FROM `", BandInfoTableName, "` WHERE kind='snr' AND `set`='", Dataset, "' AND snr>7.0 AND frequency_bin>10 AND frequency_bin<491"))[,1])


L<-list()
for(first_bin in high_bands) {
	coincidences<-dbGetQuery(con, p("SELECT * FROM `", BandInfoTableName, "` a, `", BandInfoTableName, "` b, `", BandInfoTableName, "` c WHERE a.kind='snr' AND b.kind='snr' AND c.kind='snr' AND a.`set`='", Segment, "_all'  AND b.`set`='", segment, "_LHO'  AND c.`set`='", segment, "_LLO' AND a.snr>7.0 AND b.snr>5.0 AND c.snr>5.0 AND a.snr>b.snr AND a.snr>c.snr AND abs(a.frequency-b.frequency)<", FrequencyTolerance, " AND abs(a.frequency-c.frequency)<", FrequencyTolerance, " AND abs(a.`dec`-b.`dec`)<", LocationTolerance, " AND abs(a.`dec`-c.`dec`)<", LocationTolerance, " AND abs(a.spindown-b.spindown)<", SpindownTolerance, " AND abs(a.spindown-c.spindown)<", SpindownTolerance, " AND a.first_bin=", first_bin, " AND b.first_bin=", first_bin, " AND c.first_bin=", first_bin))

	if(dim(coincidences)[1]<1)next
	names(coincidences)<-p(names(coincidences), rep(c("", "_H1", "_L1"), each=dim(coincidences)[2]/3))
	FDist<- cosdist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_H1"], coincidences[,"dec_H1"])>cos(LocationTolerance) &
		cosdist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_L1"], coincidences[,"dec_L1"])>cos(LocationTolerance)
	FDist[is.na(FDist)]<-FALSE
	coincidences<-coincidences[FDist,,drop=FALSE]
	if(dim(coincidences)[1]<1)next

	cat("Found", dim(coincidences)[1], "matches in band", first_bin/sft_length.0, "\n")
	L[[length(L)+1]]<-coincidences
	}

coincidences<-do.call(rbind, L)

cat(con=LOG, "Found", dim(coincidences)[1], "coincidences\n")

#
# Note: full S5 uses spindown reference time of 846885755
#  Early S5 used 818880650
#


line.list<-read.table(textConnection("Frequency fdot Width Comment
 46.70 0 0.1 H1CalibrationLine46.7
 54.70 0 0.1 L1CalibrationLine54.7
 265.58 -4.15e-12 0.1 InjPulsar0
 849.06 -3e-10 1   InjPulsar1
 575.16 -1.37e-13 0.5 InjPulsar2
 108.86 -1.46e-17 0.1 InjPulsar3
1401.55 -2.54e-8 0.5 InjPulsar4
  52.80 -4.03e-18 0.1 InjPulsar5
 148.29 -6.73e-9 2   InjPulsar6
1220.90 -1.12e-9 0.2 InjPulsar7
 193.75  -8.65e-9 1.25 InjPulsar8
 763.85 -1.4E-17 0.1 InjPulsar9
 128.00 0 0.1 128Hz
 256.10 0 0.25 256Hz
 350 0 20 ViolinMode
 393.10 0 0.25 H1CalibrationLine393.1
 396.70 0 0.25 L1CalibrationLine396.7
 1144.3000 0 0.25 H1CalibrationLine1144.3
 1151.9118 0 0.25 L1CalibrationLine1151.9"), header=TRUE)

# Update frequencies taking changed reference time into account:
line.list[,"Frequency"]<-line.list[,"Frequency"]+line.list[,"fdot"]*(846885755-818880650)


ExtraStats<-data.frame(Idx=0, Comment="", SNR.H1L1=coincidences$snr,
		SNR.H1=coincidences$snr_H1,
		SNR.L1=coincidences$snr_L1,
		f0=coincidences$frequency,
		fdot=coincidences$spindown,
		ra=coincidences$ra,
		dec=coincidences$dec,
		fdist=pmax(abs(coincidences$frequency-coincidences$frequency_H1), abs(coincidences$frequency-coincidences$frequency_H1), na.rm=TRUE),
		sdist=pmax(abs(coincidences$spindown-coincidences$spindown_H1), abs(coincidences$spindown-coincidences$spindown_H1), na.rm=TRUE),
		dist=acos(pmax(cosdist(coincidences$ra, coincidences$dec, coincidences$ra_H1, coincidences$dec_H1), cosdist(coincidences$ra, coincidences$dec, coincidences$ra_L1, coincidences$dec_L1), na.rm=TRUE))
		)

ExtraStats[,"Idx"]<-1:dim(ExtraStats)[1]

ExtraStats[,"Comment"]<-""

for(i in 1:dim(line.list)[1]) {
        F<- abs(ExtraStats[,"f0"]-line.list[i, "Frequency"])<line.list[i,"Width"]
        ExtraStats[F, "Comment"]<- p(as.character(line.list[i, "Comment"]), " ")
        }

F<-coincidences[,"skyband_name"]=="Lines"
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "Lines ")

F<- abs( (coincidences$frequency+1.25) %% 60.0 ) < 2.5
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "60 Hz ")

F<- abs( (coincidences$frequency-64)) < 0.01
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "64 Hz ")

F<- abs( (coincidences$frequency-128)) < 0.01
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "64 Hz ")

F<- coincidences$max_m1_neg>0.44
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m1 ")

F<- coincidences$min_m4<1.95
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m4 ")

F<- coincidences$weight_loss_fraction_H1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_H1 ")

F<- coincidences$weight_loss_fraction_L1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_L1 ")

F<- coincidences$frequency_bin<20
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "leftmost bin ")

F<- coincidences$frequency_bin>480
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "rightmost bin ")

write.table(data.frame(ExtraStats, coincidences), p(output.dir, "/band_info_snr_coincidences.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

start_plot("band_info_snr_coincidences_f0fdot")
print(xyplot(spindown~frequency, coincidences))
dev.off()

start_plot("band_info_snr_coincidences_sky")
print(xyplot(I(-dec)~I(-ra), coincidences))
dev.off()

close(LOG)
