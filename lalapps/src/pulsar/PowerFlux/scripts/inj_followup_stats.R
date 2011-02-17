library("lattice")

p<-function(...)paste(..., sep="")

source("params.R")

Params<-read.table("../params.txt", header=TRUE)
Input<-read.table("input_matches.csv", header=TRUE)
Input<-merge(Params, Input, by="i", all=TRUE, suffixes=c("_inj", ""))
Input[,"line_id"]<-as.character(Input[,"line_id"])


Output<-read.table("followup_matches.csv", header=TRUE)
Output<-merge(Params, Output, by="i", all=TRUE, suffixes=c("_inj", ""))
Output[,"line_id"]<-as.character(Output[,"line_id"])

cosdist<-function(ra1, dec1, ra2, dec2) (sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))
dist<-function(ra1, dec1, ra2, dec2) {
	a<-cosdist(ra1, dec1, ra2, dec2)
	return(ifelse(a>=1.0,  0.0 , ifelse(a<= -1, -pi, acos(a))))
	}

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

Input[,"dist"]<-dist(Input[,"ra"], Input[,"dec"], Input[,"ra_inj"], Input[,"dec_inj"])
Input[is.na(Input[,"dist"]), "dist"]<-0.3
Input[,"h0_rel"]<-Input[,"h0_inj"]/FoundUpperLimit

Output[,"dist"]<-dist(Output[,"ra_orig"], Output[,"dec_orig"], Output[,"ra_inj"], Output[,"dec_inj"])
Output[is.na(Output[,"dist"]), "dist"]<-0.3
Output[,"h0_rel"]<-Output[,"h0_inj"]/FoundUpperLimit

ROC_table<-function(table, col="h0_inj", group.func=function(x)return(x), groups) {
	table[,"Found"]<- as.integer(!is.na(table[,"snr"]))
	X<-table[order(table[,"Found"], decreasing=TRUE),,drop=FALSE]
	X<-X[!duplicated(X[,"i"]),,drop=FALSE]
	Gy<-group.func(X[,col])
	Group<-findInterval(Gy, groups, rightmost.closed=TRUE)

	Y<-aggregate(X[,"Found",drop=FALSE], list(Group=Group), mean)
	return(Y)
	}

ROC_plot<-function(col="h0_inj", group.func=function(x)return(x), group.inv.func=function(x)return(x), groups=10, ...) {
	Gy<-group.func(c(Input[,col], Output[,col]))
	Groups<-seq(min(Gy, na.rm=TRUE), max(Gy, na.rm=TRUE), length.out=groups)

	ROCInput<-ROC_table(Input, col=col, group.func=group.func, groups=Groups)
	ROCOutput<-ROC_table(Output, col=col, group.func=group.func, groups=Groups)

	C<-merge(ROCInput, ROCOutput, by="Group", suffixes=c("_input", "_output"), all=TRUE)
	C[,col]<-group.inv.func(Groups[C[,"Group"]])

	C[,"Found_input"]<-C[,"Found_input"]*100
	C[,"Found_output"]<-C[,"Found_output"]*100
	print(xyplot(as.formula(p("Found_input+Found_output~", col)), C, pch=c(3, 1), cex=1, ...))
	}

ComparisonPlot<-function(formula, decreasing=TRUE, best.snr=FALSE, auto.key=list(text=c("Input", "Output"), columns=2), pch=c(3, 1), ...) {
	C<-merge(Input, Output, by.x=c("i", "line_id"), by.y=c("i", "line_id_orig"), suffixes=c("_input", "_output"))
	C<-C[!is.na(C[,"line_id"]),]
	if(best.snr) {
		C<-C[order(C[,"snr_output"], decreasing=TRUE),,drop=FALSE]
		C<-C[!duplicated(C[,"i"]),,drop=FALSE]
		}
	print(xyplot(formula, C, auto.key=auto.key, pch=pch, ...))
	}

make_plot<-function(name, width=600, height=600, dpi=100, pointsize=18, ...) {
	png(p(name, ".png"), width=width, height=height, res=dpi, pointsize=18, ...)
	}

# make_plot<-function(name, width=600, height=600, dpi=100, pointsize=18, ...) {
# 	pdf(p(name, ".pdf"), width=width*5/600, height=height*5/600, bg="white", ...)
# 	}


make_plot("injection_recovery")
ROC_plot(group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0", ylab="% found")
dev.off()

make_plot("injection_recovery_rel")
ROC_plot(col="h0_rel", group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0 relative to upper limit", ylab="% found")
dev.off()

make_plot("injection_recovery_by_f0")
ROC_plot(col="f0_inj", auto.key=list(columns=2), xlab="Injection frequency", ylab="% found")
dev.off()

make_plot("injection_recovery_by_spindown")
ROC_plot(col="spindown_inj", auto.key=list(columns=2), xlab="Injection spindown", ylab="% found")
dev.off()

make_plot("injection_recovery_by_distance")
ROC_plot(col="dist", auto.key=list(columns=2), xlab="Distance from injection", ylab="% found")
dev.off()


make_plot("snr_improvement")
ComparisonPlot(I(snr_output/snr_input)~h0_input, auto.key=FALSE, xlab="h0", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_rel")
ComparisonPlot(I(snr_output/snr_input)~h0_rel_input, auto.key=FALSE, xlab="h0 relative to upper limit", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_zoomed")
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~h0_input, auto.key=FALSE, xlab="h0", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_rel_zoomed")
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~h0_rel_input, auto.key=FALSE, xlab="h0 relative to upper limit", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_vs_snr_zoomed")
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~log10(snr_input), auto.key=FALSE, xlab="log10(snr)", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_vs_iota_zoomed")
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~iota_inj_input, auto.key=FALSE, xlab="iota", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("snr_improvement_vs_dist_zoomed")
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~dist(ra_output, dec_output, ra_inj_input, dec_inj_input), auto.key=FALSE, xlab="dist_output", ylab="SNR Output/SNR input", pch=3)
dev.off()

make_plot("f0_improvement")
ComparisonPlot(I(abs(f0_input-f0_inj_input))+I(abs(f0_output-f0_inj_input))~h0_input, xlab="h0", ylab="Frequency mismatch, Hz", best.snr=TRUE)
dev.off()

make_plot("f0_improvement_rel")
ComparisonPlot(I(abs(f0_input-f0_inj_input))+I(abs(f0_output-f0_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Frequency mismatch, Hz", best.snr=TRUE)
dev.off()

make_plot("f0_improvement_zoomed")
ComparisonPlot(I(abs(f0_input-f0_inj_input))+I(abs(f0_output-f0_inj_input))~h0_input, xlab="h0", ylab="Frequency mismatch, Hz", best.snr=TRUE, ylim=c(0, 2*FrequencyTolerance))
dev.off()

make_plot("f0_improvement_rel_zoomed")
ComparisonPlot(I(abs(f0_input-f0_inj_input))+I(abs(f0_output-f0_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Frequency mismatch, Hz", best.snr=TRUE, ylim=c(0, 2*FrequencyTolerance))
dev.off()

make_plot("spindown_improvement")
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_input, xlab="h0", ylab="Spindown mismatch, Hz/s", best.snr=TRUE)
dev.off()

make_plot("spindown_improvement_rel")
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Spindown mismatch, Hz/s", best.snr=TRUE)
dev.off()

make_plot("spindown_improvement_zoomed")
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_input, xlab="h0", ylab="Spindown mismatch, Hz/s", best.snr=TRUE, ylim=c(0, 2*SpindownTolerance))
dev.off()

make_plot("spindown_improvement_rel_zoomed")
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Spindown mismatch, Hz/s", best.snr=TRUE, ylim=c(0, 2*SpindownTolerance))
dev.off()

make_plot("distance_improvement")
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_input, xlab="h0", ylab="Distance from injection, rad", best.snr=TRUE)
dev.off()

make_plot("distance_improvement_rel")
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Distance from injection, rad", best.snr=TRUE)
dev.off()

make_plot("distance_improvement_zoomed")
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_input, xlab="h0", ylab="Distance from injection, rad", best.snr=TRUE, ylim=c(0, 2*LocationTolerance))
dev.off()

make_plot("distance_improvement_rel_zoomed")
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Distance from injection, rad", best.snr=TRUE, ylim=c(0, 2*LocationTolerance))
dev.off()

#
# Try to extract phases if they exist in this run
#
try({
Output[,"phase"]<-as.integer(as.character(gsub("^[^_]*_([^_/]*)/.*", "\\1", Output[,"tag"])))
make_plot("phase_primary")
print(xyplot(phase~h0_inj, Output[Output[,"primary"],,drop=FALSE], scales=list(x=list(log=TRUE))))
dev.off()
})
