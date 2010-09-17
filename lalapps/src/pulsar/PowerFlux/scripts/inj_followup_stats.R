library("lattice")

p<-function(...)paste(..., sep="")

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

	C<-merge(ROCInput, ROCOutput, by="Group", suffixes=c("_input", "_output"))
	C[,col]<-group.inv.func(Groups[C[,"Group"]])

	C[,"Found_input"]<-C[,"Found_input"]*100
	C[,"Found_output"]<-C[,"Found_output"]*100
	print(xyplot(as.formula(p("Found_input+Found_output~", col)), C, ...))
	}

ComparisonPlot<-function(formula, decreasing=TRUE, best.snr=FALSE, ...) {
	C<-merge(Input, Output, by.x=c("i", "line_id"), by.y=c("i", "line_id_orig"), suffixes=c("_input", "_output"))
	C<-C[!is.na(C[,"line_id"]),]
	if(best.snr) {
		C<-C[order(C[,"snr_output"], decreasing=TRUE),,drop=FALSE]
		C<-C[!duplicated(C[,"i"]),,drop=FALSE]
		}
	print(xyplot(formula, C, ...))
	}

png("injection_recovery.png", width=600, height=600)
ROC_plot(group.func=log10, group.inv.func=function(x)return(10^x), auto.key=list(columns=2), xlab="h0", ylab="% found")
dev.off()

png("snr_improvement.png", width=600, height=600)
ComparisonPlot(I(snr_output/snr_input)~h0_input, xlab="h0", ylab="SNR Output/SNR input")
dev.off()

png("snr_improvement_zoomed.png", width=600, height=600)
ComparisonPlot(I(pmin(snr_output/snr_input, 3))~h0_input, xlab="h0", ylab="SNR Output/SNR input")
dev.off()

png("f0_improvement.png", width=600, height=600)
ComparisonPlot(I(abs(f0_input-f0_inj_input))+I(abs(f0_output-f0_inj_input))~h0_input, auto.key=list(columns=2), xlab="h0", ylab="Frequency mismatch, Hz", best.snr=TRUE)
dev.off()

png("spindown_improvement.png", width=600, height=600)
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_input, auto.key=list(columns=2), xlab="h0", ylab="Spindown mismatch, Hz/s", best.snr=TRUE)
dev.off()

png("spindown_improvement.png_zoomed", width=600, height=600)
ComparisonPlot(I(abs(spindown_input-spindown_inj_input))+I(abs(spindown_output-spindown_inj_input))~h0_input, auto.key=list(columns=2), xlab="h0", ylab="Spindown mismatch, Hz/s", best.snr=TRUE, ylim=c(0, 2e-11))
dev.off()

png("distance_improvement.png", width=600, height=600)
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_input, auto.key=list(columns=2), xlab="h0", ylab="Distance from injection, rad", best.snr=TRUE)
dev.off()

png("distance_improvement_zoomed.png", width=600, height=600)
ComparisonPlot(I(dist(ra_input, dec_input, ra_inj_input, dec_inj_input))+I(dist(ra_output, dec_output, ra_inj_input, dec_inj_input))~h0_input, auto.key=list(columns=2), xlab="h0", ylab="Distance from injection, rad", best.snr=TRUE, ylim=c(0, 0.1))
dev.off()
