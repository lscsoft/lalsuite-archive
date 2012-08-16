library("lattice")

p<-function(...)paste(..., sep="")
ExcludeMissingUL<-TRUE

source("params.R")

cosdist<-function(ra1, dec1, ra2, dec2) (sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))
dist<-function(ra1, dec1, ra2, dec2) {
	a<-cosdist(ra1, dec1, ra2, dec2)
	return(ifelse(a>=1.0,  0.0 , ifelse(a<= -1.0, -pi, acos(a))))
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

Params<-read.table("../params.txt", header=TRUE)
UL<-read.table("../upper_limits.csv", header=TRUE)
UL[,"non_gaussian"]<- UL[,"max_m1_neg"]>0.44 | UL[,"min_m4"]<1.6

Params<-merge(Params, UL, by.x="i", by.y="instance", suffixes=c("", "_UL"), all.x=TRUE)
Params[,"h0_rel"]<-Params[,"h0"]/Params[,"ul"]
Params[,"ecl_cos"]<-dot2(Params[,"ra"], Params[,"dec"], ecliptic_pole)

Params[,"h0_rel_gaussian"]<-ifelse(Params[,"non_gaussian"], NA, Params[,"h0_rel"])
Params[,"h0_rel_gaussian_B"]<-ifelse(Params[,"non_gaussian"] | abs(Params[,"f0"]-615)<10, NA, Params[,"h0_rel"])
Params[,"h0_rel_non_gaussian"]<-ifelse(!Params[,"non_gaussian"], NA, Params[,"h0_rel"])

Params[,"ecl_cos_large"]<-ifelse(Params[,"h0_rel"]>=1, Params[,"ecl_cos"], NA)
Params[,"ecl_cos_large_B"]<-ifelse(Params[,"h0_rel"]>=1 & abs(Params[,"f0"]-615)>10, Params[,"ecl_cos"], NA)
Params[,"spindown_large"]<-ifelse(Params[,"h0_rel"]>=1, Params[,"spindown"], NA)
Params[,"f0_large"]<-ifelse(Params[,"h0_rel"]>=1, Params[,"f0"], NA)

Params[,"ecl_cos_large2"]<-ifelse(Params[,"h0_rel"]>=2, Params[,"ecl_cos"], NA)
Params[,"spindown_large2"]<-ifelse(Params[,"h0_rel"]>=2, Params[,"spindown"], NA)
Params[,"f0_large2"]<-ifelse(Params[,"h0_rel"]>=2, Params[,"f0"], NA)

Params[,"ecl_cos_large2_gaussian"]<-ifelse(Params[,"h0_rel"]>=2 & !Params[,"non_gaussian"], Params[,"ecl_cos"], NA)
Params[,"spindown_large2_gaussian"]<-ifelse(Params[,"h0_rel"]>=2 & !Params[,"non_gaussian"], Params[,"spindown"], NA)
Params[,"f0_large2_gaussian"]<-ifelse(Params[,"h0_rel"]>=2 & !Params[,"non_gaussian"], Params[,"f0"], NA)

Input<-read.table("input_matches.csv", header=TRUE)
Input<-merge(Params, Input, by="i", all=TRUE, suffixes=c("_inj", ""))
Input[,"line_id"]<-as.character(Input[,"line_id"])


Output<-read.table("followup_matches.csv", header=TRUE)
Output<-merge(Params, Output, by="i", all=TRUE, suffixes=c("_inj", ""))
Output[,"line_id"]<-as.character(Output[,"line_id"])

if(ExcludeMissingUL) {
	Params<-Params[!is.na(Params[,"h0_rel"]),,drop=FALSE]
	Input<-Input[!is.na(Input[,"h0_rel"]),,drop=FALSE]
	Output<-Output[!is.na(Output[,"h0_rel"]),,drop=FALSE]
	}

Input[,"dist"]<-dist(Input[,"ra"], Input[,"dec"], Input[,"ra_inj"], Input[,"dec_inj"])
Input[is.na(Input[,"dist"]), "dist"]<-0.3

Output[,"dist"]<-dist(Output[,"ra_orig"], Output[,"dec_orig"], Output[,"ra_inj"], Output[,"dec_inj"])
Output[is.na(Output[,"dist"]), "dist"]<-0.3

C<-merge(Input[,setdiff(names(Input), "line_id_orig"),drop=FALSE], Output, by.x=c("i", "line_id"), by.y=c("i", "line_id_orig"), suffixes=c("_input", "_output"), all=TRUE)
missing_injections<-C[is.na(C[,"snr_output"]) & !is.na(C[,"snr_input"]),,drop=FALSE]

write.table(missing_injections, "missing_injections.csv", col.names=TRUE, row.names=FALSE)


ROC_table<-function(table, col="h0_inj", group.func=function(x)return(x), groups) {
	table[,"Found"]<- as.integer(!is.na(table[,"snr"]))
	X<-table[order(table[,"Found"], decreasing=TRUE),,drop=FALSE]
	X<-X[!duplicated(X[,"i"]),,drop=FALSE]
	Gy<-group.func(X[,col])
	Group<-findInterval(Gy, groups, rightmost.closed=TRUE, all.inside=TRUE)

	# Aggregate orders Group column right now - but just to be sure in case later R versions do it differently
	Y<-aggregate(X[,"Found",drop=FALSE], list(Group=Group), mean)
	Y<-Y[order(Y[,1]),,drop=FALSE]
	Y2<-aggregate(X[,"Found",drop=FALSE], list(Group=Group), sd)
	Y2<-Y2[order(Y2[,1]),,drop=FALSE]
	Y3<-aggregate(X[,"Found",drop=FALSE], list(Group=Group), length)
	Y3<-Y3[order(Y3[,1]),,drop=FALSE]

	return(cbind(Y, data.frame(Found_sd=Y2[,2]/sqrt(Y3[,2]))))
	}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	panel.arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
	}

ROC_plot<-function(col="h0_inj", group.func=function(x)return(x), group.inv.func=function(x)return(x), groups=10, error.bars=FALSE, ...) {
	Gy<-group.func(c(Input[,col], Output[,col]))
	#Groups<-seq(min(Gy, na.rm=TRUE), max(Gy, na.rm=TRUE), length.out=groups)
	Groups<-quantile(Gy, (0:groups)/groups, na.rm=TRUE)

	ROCInput<-ROC_table(Input, col=col, group.func=group.func, groups=Groups)
	ROCOutput<-ROC_table(Output, col=col, group.func=group.func, groups=Groups)

	C<-merge(ROCInput, ROCOutput, by="Group", suffixes=c("_input", "_output"), all=TRUE)
	C[,col]<-group.inv.func(Groups[C[,"Group"]])

	C[,"Found_input"]<-C[,"Found_input"]*100
	C[,"Found_output"]<-C[,"Found_output"]*100
	C[,"Found_sd_input"]<-C[,"Found_sd_input"]*100
	C[,"Found_sd_output"]<-C[,"Found_sd_output"]*100
	print(xyplot(as.formula(p("Found_input+Found_output~", col)), C, pch=c(3, 1), par.settings=list(superpose.symbol=list(pch=c(3,1,4))), cex=1, col=c("blue", "magenta"), ...))
	trellis.focus("panel", row=1, column=1)
	if(error.bars) {
		print(error.bar(C[,col], C[,"Found_input"], C[,"Found_sd_input"], col="blue"))
		print(error.bar(C[,col], C[,"Found_output"], C[,"Found_sd_output"], col="magenta"))
		}
	trellis.unfocus()
	}

ComparisonPlot<-function(formula, decreasing=TRUE, best.snr=FALSE, omit.found=FALSE, auto.key=list(text=c("Input", "Output"), columns=2), pch=c(3, 1), ...) {
	C<-merge(Input[,setdiff(names(Input), "line_id_orig"),drop=FALSE], Output, by.x=c("i", "line_id"), by.y=c("i", "line_id_orig"), suffixes=c("_input", "_output"), all=omit.found)
	if(!omit.found)C<-C[!is.na(C[,"line_id"]),]
	if(best.snr) {
		C<-C[order(C[,"snr_output"], decreasing=TRUE),,drop=FALSE]
		C<-C[!duplicated(C[,"i"]),,drop=FALSE]
		}
	if(omit.found) {
		C<-C[is.na(C[,"snr_output"]) & !is.na(C[,"snr_input"]),,drop=FALSE]
		}
	print(xyplot(formula, C, auto.key=auto.key, pch=pch, par.settings=list(superpose.symbol=list(pch=c(3,1,4))), ...))
	}

 make_plot<-function(name, width=600, height=600, dpi=100, pointsize=18, ...) {
 	png(p(name, ".png"), width=width, height=height, res=dpi, pointsize=18, ...)
 	}

#make_plot<-function(name, width=600, height=600, dpi=100, pointsize=18, ...) {
#	pdf(p(name, ".pdf"), width=width*5/600, height=height*5/600, bg="white", ...)
#	}


make_plot("injection_recovery")
ROC_plot(group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0", ylab="% found", groups=60)
dev.off()

make_plot("injection_recovery_rel")
ROC_plot(col="h0_rel", group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0 relative to upper limit", ylab="% found", groups=40, panel=function(x,y,...) { panel.abline(v=1, col="green4"); panel.abline(h=95, col="green4"); panel.xyplot(x,y,...) })
dev.off()

make_plot("injection_recovery_rel_gaussian")
ROC_plot(col="h0_rel_gaussian", group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0 relative to upper limit", ylab="% found", groups=40, panel=function(x,y,...) { panel.abline(v=1, col="green4"); panel.abline(h=95, col="green4"); panel.xyplot(x,y,...) })
dev.off()

make_plot("injection_recovery_rel_gaussian_B")
ROC_plot(col="h0_rel_gaussian_B", group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0 relative to upper limit", ylab="% found", groups=40, panel=function(x,y,...) { panel.abline(v=1, col="green4"); panel.abline(h=95, col="green4"); panel.xyplot(x,y,...) })
dev.off()

make_plot("injection_recovery_rel_non_gaussian")
ROC_plot(col="h0_rel_non_gaussian", group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0 relative to upper limit", ylab="% found", groups=20, panel=function(x,y,...) { panel.abline(v=1, col="green4"); panel.abline(h=75, col="green4"); panel.xyplot(x,y,...) })
dev.off()

make_plot("injection_recovery_by_f0")
ROC_plot(col="f0_inj", auto.key=list(columns=2), xlab="Injection frequency", ylab="% found")
dev.off()

make_plot("injection_recovery_by_f0_large")
ROC_plot(col="f0_large", auto.key=list(columns=2), xlab="Injection frequency", ylab="% found")
dev.off()

make_plot("injection_recovery_by_f0_large2")
ROC_plot(col="f0_large2", auto.key=list(columns=2), xlab="Injection frequency", ylab="% found")
dev.off()

make_plot("injection_recovery_by_spindown")
ROC_plot(col="spindown_inj", auto.key=list(columns=2), xlab="Injection spindown", ylab="% found")
dev.off()

make_plot("injection_recovery_by_spindown_large")
ROC_plot(col="spindown_large", auto.key=list(columns=2), xlab="Injection spindown", ylab="% found")
dev.off()

make_plot("injection_recovery_by_spindown_large2")
ROC_plot(col="spindown_large2", auto.key=list(columns=2), xlab="Injection spindown", ylab="% found")
dev.off()

make_plot("injection_recovery_by_distance")
ROC_plot(col="dist", auto.key=list(columns=2), xlab="Distance from injection", ylab="% found")
dev.off()

make_plot("injection_recovery_by_ecliptic_pole_projection")
ROC_plot(col="ecl_cos", auto.key=list(columns=2), xlab="Projection on ecliptic axis", ylab="% found")
dev.off()

make_plot("injection_recovery_by_ecliptic_pole_projection_large")
ROC_plot(col="ecl_cos_large", auto.key=list(columns=2), xlab="Projection on ecliptic axis", ylab="% found")
dev.off()

make_plot("injection_recovery_by_ecliptic_pole_projection_large_B")
ROC_plot(col="ecl_cos_large_B", auto.key=list(columns=2), xlab="Projection on ecliptic axis", ylab="% found")
dev.off()

make_plot("injection_recovery_by_ecliptic_pole_projection_large2")
ROC_plot(col="ecl_cos_large2", auto.key=list(columns=2), xlab="Projection on ecliptic axis", ylab="% found")
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
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_input, xlab="h0", ylab="Distance from injection, rad", best.snr=TRUE, ylim=c(0, 2*LocationToleranceMax))
dev.off()

make_plot("distance_improvement_rel_zoomed")
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_rel_input, xlab="h0 relative to upper limit", ylab="Distance from injection, rad", best.snr=TRUE, ylim=c(0, 2*LocationToleranceMax))
dev.off()

#
# Characterize missing injections, if any
#
make_plot("missing_injections_sky")
ComparisonPlot(dec_inj_input~ra_inj_input, xlab="RA", ylab="DEC", best.snr=TRUE, omit.found=TRUE, pch=3)
dev.off()

make_plot("missing_injections_f0_fdot")
ComparisonPlot(spindown_inj_input~f0_inj_input, xlab="f0", ylab="spindown", best.snr=TRUE, omit.found=TRUE, pch=3)
dev.off()

make_plot("missing_injections_alignment")
ComparisonPlot(iota_inj_input~psi_inj_input, xlab="psi", ylab="iota", best.snr=TRUE, omit.found=TRUE, pch=3)
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
