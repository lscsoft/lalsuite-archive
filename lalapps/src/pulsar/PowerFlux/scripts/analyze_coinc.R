source("params.R")

require("lattice")
require("RMySQL")

p<-function(...) {
	return(paste(..., sep=""))
	}

cosdist<-function(ra1, dec1, ra2, dec2) (sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))
dist<-function(ra1, dec1, ra2, dec2) {
	a<-cosdist(ra1, dec1, ra2, dec2)
	return(ifelse(a>=1.0,  0.0 , ifelse(a<= -1.0, pi, acos(a))))
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

# Adapted from Joe's matlab script
ifo_A2<-function(ra, dec, iota, psi, gamma, lambda) {
	# equations from JKS
	j3<- (1.0/128)*(28-44*cos(lambda)^2+5*sin(2*gamma)^2*cos(lambda)^4)
	j2<- (1.0/1024)*(68-20*cos(lambda)^2-13*sin(2*gamma)^2*cos(lambda)^4)
	j1<- (1.0/256)*(4-20*cos(lambda)^2+35*sin(2*gamma)^2*cos(lambda)^4)

	e2<- 4*j2-j3*cos(2*dec)+j1*cos(2*dec)^2
	e1<- 4*j1*cos(dec)^4

	G2<- (1.0/4)*(1+6*cos(iota)^2+cos(iota)^4)
	F2<- (1.0/4)*(1-cos(iota)^2)^2

	A2<- F2*e1*cos(4.0*psi)+G2*e2

	return(A2)
	}

estimated_fstat_timebase<-function(h0, ra, dec, iota, psi, M.H1, M.L1,  twoFcutoff=60, H1.uptime=0.75, L1.uptime=0.75, sft.coherence=1800) {
	# LHO position
	gamma_LHO<-171.8*pi/180
	lambda_LHO<-46.45*pi/180

	A2_LHO<- ifo_A2(ra, dec, iota, psi, gamma_LHO, lambda_LHO)
	
	# LHO position
	gamma_LLO<-243.0*pi/180
	lambda_LLO<- 30.56*pi/180

	A2_LLO<- ifo_A2(ra, dec, iota, psi, gamma_LLO, lambda_LLO)

	Sh.H1<- 0.5*M.H1^2*sft.coherence
	Sh.L1<- 0.5*M.L1^2*sft.coherence

	timebase<- twoFcutoff/(A2_LHO*h0^2*H1.uptime/Sh.H1+A2_LLO*h0^2*L1.uptime/Sh.L1)
	return(timebase)
	}

con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

# override standard dbQuery as it has issues
dbGetQuery<-function(con, query) {
	res<-dbSendQuery(con, query)
	L<-list()
	while(!dbHasCompleted(res)) {
		a<-fetch(res, 10000)

		L[[length(L)+1]]<-a
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
		a<-fetch(res, 20000)

		coincidences<-a

		if(length(coincidences)[1]<1 || dim(coincidences)[1]<1)next

		names(coincidences)<-p(names(coincidences), rep(c("", "_H1", "_L1"), each=dim(coincidences)[2]/3))
		FDist<- ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_H1"], coincidences[,"dec_H1"])<(LocationTolerance) &
			ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_L1"], coincidences[,"dec_L1"])<(LocationTolerance) 
		FDist[is.na(FDist)]<-FALSE
		coincidences<-coincidences[FDist,,drop=FALSE]
		if(dim(coincidences)[1]<1)next

#		if(dim(coincidences)[1]>1000) {
#			cat("truncating\n")
#			coincidences<-coincidences[order(-coincidences$snr, -pmin(coincidences$snr_H1, coincidences$snr_L1), -pmax(coincidences$snr_H1, coincidences$snr_L1)),,drop=FALSE][1,,drop=FALSE]
#			}

		L[[length(L)+1]]<-coincidences
		}
	dbClearResult(res)
	if(length(L)<1) return(invisible(NULL))
		else
	if(length(L)>1) {
		return(do.call(rbind, L))
		} else return(L[[1]])
	}

create_index<-function() {
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ADD COLUMN F0INDEX INTEGER NOT NULL"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " MODIFY COLUMN `set` VARCHAR(255)"))
	dbGetQuery(con, p("UPDATE ", CoincTableName, " SET F0INDEX=ROUND(FREQUENCY/", FrequencyTolerance, ")"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ADD INDEX F0IDX(F0INDEX, `set`)"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ORDER BY F0INDEX"))
	}


output.dir<-p("report/coinc_", Segment)

dir.create(output.dir, recursive=TRUE)

LOG<-file(p(output.dir, "/coinc.log"), open="w")

Lines_skyband<-dbGetQuery(con, p("SELECT skyband FROM `", BandInfoTableName, "` WHERE skyband_name='Lines' LIMIT 1"))[1,1]

strong_lines<-dbGetQuery(con, p("SELECT dataset, strength, line_bin, flag FROM `", LineInfoTableName, "` WHERE strength>", StrongLine))

start_png<-function(name, width=600, height=600) {
	png(p(output.dir, "/", name), width=width, height=height)
	}

start_png("strong_lines.png", width=1200, height=600)
print(xyplot(I(strength*LineStrengthToSNR)~I(line_bin/1800.0)|dataset, strong_lines, xlab="Frequency", ylab="SNR"))
dev.off()

OC<-dbGetQuery(con, p("SELECT `label`,  `set`, SUM(`count`) as `count`, SUM(IF(`index`>3, `count`, 0)) as `count4`, SUM(IF(`index`>5, `count`, 0)) as `count6` FROM ", OutlierHistogramTableName, " WHERE `set` LIKE '", Segment, "%' AND skyband != ", Lines_skyband, " GROUP BY `label`, `set`"))

start_png("OC4_ratio.png", width=1200, height=400)
xyplot(I(count4/count)~label|set, OC)
dev.off()

start_png("OC6_ratio.png", width=1200, height=400)
xyplot(I(count6/count)~label|set, OC)
dev.off()

#InstanceList<- as.integer(as.character(dbGetQuery(con, p("SELECT DISTINCT label FROM `", CoincTableName, "`"))[,1]))

#SkyBands<-dbGetQuery(con, p("SELECT DISTINCT skyband_name FROM `", TableName, "`"))[,1]


high_bands<- dbGetQuery(con, p("SELECT F0INDEX, MAX(snr) as max_snr, COUNT(*) as `count` FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_all' AND snr>", AllSNRCutoff, " GROUP BY F0INDEX"))

cat("Found", dim(high_bands)[1], "indices with SNR >", AllSNRCutoff, "\n")
print(high_bands[order(-high_bands[,"count"]),,drop=FALSE][1:10, ])

#ExcludeBands<-c(108.9, 193.4, 393.4)
ExcludeBands<-c()
for(f0 in ExcludeBands) {
	F<- abs(high_bands[,"F0INDEX"]*FrequencyTolerance-f0)<0.25
	high_bands<-high_bands[!F,,drop=FALSE]
	}
cat(dim(high_bands)[1], "remain after excluding", p(ExcludeBands, collapse=" "), "bands\n")
high_bands<-high_bands[order(high_bands[,"F0INDEX"]),,drop=FALSE]
# F<- abs( (high_bands*FrequencyTolerance+1.25) %% 60.0 ) < 2.5
# high_bands<-high_bands[!F]
# cat(length(high_bands), "remain after excluding 60 Hz bands\n")

start_png("coinc_max_snr.png", width=600, height=600)
print(xyplot(pmin(max_snr, 50)~I(F0INDEX*FrequencyTolerance), high_bands))
dev.off()

start_png("coinc_count.start_png", width=600, height=600)
print(xyplot(count~I(F0INDEX*FrequencyTolerance), high_bands))
dev.off()


high_bands[,"snr_cutoff_all"]<-pmax(high_bands[,"max_snr"]/2.0, AllSNRCutoff)
high_bands[,"snr_cutoff_ifo"]<-pmax(high_bands[,"max_snr"]/4.0, IFOSNRCutoff)

cat(sum(high_bands[,"snr_cutoff_all"]>AllSNRCutoff), "indices with elevated snr level for _all\n")
cat(sum(high_bands[,"snr_cutoff_ifo"]>IFOSNRCutoff), "indices with elevated snr level for ifos\n")


L<-list()
for(i in 1:dim(high_bands)[1]) {
	fidx<-high_bands[i, "F0INDEX"]

	a_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_all' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_all"], " AND F0INDEX=", fidx, " AND frequency_bin>=20 AND frequency_bin<=480")

	b_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LHO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND F0INDEX IN (", fidx, ", ", fidx+1, ", ", fidx-1, ")",  " AND frequency_bin>=20 AND frequency_bin<=480")

	c_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LLO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND F0INDEX IN (", fidx, ", ", fidx+1, ", ", fidx-1, ")",  " AND frequency_bin>=20 AND frequency_bin<=480")

	coincidences<-dbGetQueryCoinc(con, p("SELECT * FROM (", a_table, ") a, (", b_table, ") b, (", c_table, ") c WHERE a.snr>b.snr AND a.snr>c.snr AND abs(a.frequency-b.frequency)<", FrequencyTolerance, " AND abs(a.frequency-c.frequency)<", FrequencyTolerance, " AND abs(a.`dec`-b.`dec`)<", DecTolerance, " AND abs(a.`dec`-c.`dec`)<", DecTolerance, " AND abs(a.spindown-b.spindown)<", SpindownTolerance, " AND abs(a.spindown-c.spindown)<", SpindownTolerance))
	if(length(coincidences)<1 || dim(coincidences)[1]<1)next

# 	names(coincidences)<-p(names(coincidences), rep(c("", "_H1", "_L1"), each=dim(coincidences)[2]/3))
# 	FDist<- ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_H1"], coincidences[,"dec_H1"])<(LocationTolerance) &
# 		ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_L1"], coincidences[,"dec_L1"])<(LocationTolerance) 
# 	FDist[is.na(FDist)]<-FALSE
# 	coincidences<-coincidences[FDist,,drop=FALSE]

#	if(dim(coincidences)[1]<1)next

	cat("Found", dim(coincidences)[1], "matches in band", fidx*FrequencyTolerance, "\n")

#	cat("truncating\n")
#	coincidences<-coincidences[order(-coincidences$snr, -pmin(coincidences$snr_H1, coincidences$snr_L1), -pmax(coincidences$snr_H1, coincidences$snr_L1)),,drop=FALSE][1,,drop=FALSE]

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
 54.496 0 0.01 VME_CPU_clock
108.992 0 0.01 VME_CPU_clock
 128.00 0 0.1 128Hz
 256.10 0 0.25 256Hz
 350 0 20 ViolinMode
 393.10 0 0.25 H1CalibrationLine393.1
 396.70 0 0.25 L1CalibrationLine396.7
 1144.3000 0 0.25 H1CalibrationLine1144.3
 1151.9118 0 0.25 L1CalibrationLine1151.9"), header=TRUE)

# Update frequencies taking changed reference time into account:

ExtraStats<-data.frame(Idx=0, Comment="", SNR.H1L1=coincidences$snr,
		SNR.H1=coincidences$snr_H1,
		SNR.L1=coincidences$snr_L1,
		Max.SNR.H1L1=0,
		Max.SNR.H1=0,
		Max.SNR.L1=0,
		f0=coincidences$frequency,
		fdot=coincidences$spindown,
		ra=coincidences$ra,
		dec=coincidences$dec,
		fdist=pmax(abs(coincidences$frequency-coincidences$frequency_H1), abs(coincidences$frequency-coincidences$frequency_H1), na.rm=TRUE),
		sdist=pmax(abs(coincidences$spindown-coincidences$spindown_H1), abs(coincidences$spindown-coincidences$spindown_H1), na.rm=TRUE),
		dist=pmax(ecliptic_dist(coincidences$ra, coincidences$dec, coincidences$ra_H1, coincidences$dec_H1), ecliptic_dist(coincidences$ra, coincidences$dec, coincidences$ra_L1, coincidences$dec_L1), na.rm=TRUE),
		line.f0=0,
		line.comment="",
		OC.H1L1=0,
		OC.H1=0,
		OC.L1=0,
		FStat60.timebase=estimated_fstat_timebase(h0=coincidences[, "snr"]*coincidences[, "S"], ra=coincidences[, "ra"], dec=coincidences[, "dec"], iota=coincidences[, "iota"], psi=coincidences[, "psi"], M.H1=coincidences[, "M_H1"], M.L1=coincidences[, "M_L1"], twoFcutoff=60)/(24*3600),
		Group=0,
		Primary=0
		)

ExtraStats[,"Idx"]<-1:dim(ExtraStats)[1]

ExtraStats[,"Comment"]<-""
ExtraStats[,"line.comment"]<-""

for(i in 1:dim(line.list)[1]) {
        F<- abs(ExtraStats[,"f0"]-line.list[i, "Frequency"])<line.list[i,"Width"]
        ExtraStats[F, "Comment"]<- p(as.character(line.list[i, "Comment"]), " ")
        }

#F<-coincidences[,"skyband"] %in% Lines_bands
#ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "Lines ")

F<- abs( (coincidences$frequency+1.25) %% 60.0 ) < 2.5
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "60 Hz ")

F<- abs( (coincidences$frequency+0.05) %% 16.0 ) < 0.1
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "16 Hz ")

F<- coincidences$max_m1_neg>0.44
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m1 ")

F<- coincidences$min_m4<1.95
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m4 ")

F<- coincidences$weight_loss_fraction_H1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_H1 ")

F<- coincidences$weight_loss_fraction_L1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_L1 ")

F<- coincidences$skyband==Lines_skyband | coincidences$skyband_H1==Lines_skyband | coincidences$skyband_L1==Lines_skyband
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "lines ")

P<-order(coincidences[, "F0INDEX"])
Idx<-coincidences[P, "F0INDEX"]
Group<-Idx
for(i in 2:length(Group)) {
	if(Idx[i-1]>=Idx[i]-1)Group[i]<-Group[i-1]
	}

ExtraStats[P, "Group"]<-Group

P<-order(-ExtraStats$SNR.H1L1, -pmin(ExtraStats$SNR.H1, ExtraStats$SNR.L1), ExtraStats$sdist, ExtraStats$dist, ExtraStats$fdist)


ExtraStats[P,"Primary"]<- !duplicated(ExtraStats[P, "Group"])

for(i in 1:dim(ExtraStats)[1]) {

	if(ExtraStats[i, "Comment"]!="" && ExtraStats[i, "Primary"]<1)next
	cat(".")
	FOINDEX_list<- p((coincidences[i, "F0INDEX"]-MaxF0INDEXTolerance):(coincidences[i, "F0INDEX"]+MaxF0INDEXTolerance), collapse=", ")
	ra<-coincidences[i, "ra"]
	dec<-coincidences[i, "dec"]

	ExtraStats[i, "Max.SNR.H1L1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", coincidences[i, 'set'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND F0INDEX IN (", FOINDEX_list, ") AND abs(frequency-", coincidences[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", coincidences[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]

	ExtraStats[i, "Max.SNR.H1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", coincidences[i, 'set_H1'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND F0INDEX IN (", FOINDEX_list, ") AND abs(frequency-", coincidences[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", coincidences[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]

	ExtraStats[i, "Max.SNR.L1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", coincidences[i, 'set_L1'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND F0INDEX IN (", coincidences[i, "F0INDEX"], ", ", coincidences[i, "F0INDEX"]+1, ", ", coincidences[i, "F0INDEX"]-1, ") AND abs(frequency-", coincidences[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", coincidences[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]

	ExtraStats[i, "OC.H1L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set'], "'AND `skyband`!=", Lines_skyband))[1,1]

	ExtraStats[i, "OC.H1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set_H1'], "' AND `skyband`!=", Lines_skyband))[1,1]

	ExtraStats[i, "OC.L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set_L1'], "'AND `skyband`!=", Lines_skyband))[1,1]

	bin_width<-1800.0*(coincidences[i, "frequency"]*1e-4+abs(coincidences[i, "spindown"]*0.5*run_timebase))

	joint_bin_width<-max(bin_width, VeryStrongLineBinWidth)

	lines<-dbGetQuery(con, p("SELECT dataset, strength, line_bin, flag FROM `", LineInfoTableName, "` WHERE ((strength>", ExtraStats[i, "SNR.H1L1"]*OutlierSNRToLineStrength, " AND line_bin>", coincidences[i, "frequency"]*1800-bin_width, 
		" AND line_bin<", coincidences[i, "frequency"]*1800+bin_width, ") OR (strength>", VeryStrongLine, " AND line_bin>",  coincidences[i, "frequency"]*1800-VeryStrongLineBinWidth, 
		" AND line_bin<", coincidences[i, "frequency"]*1800+VeryStrongLineBinWidth, ")) AND line_bin>", coincidences[i, "frequency"]*1800-joint_bin_width, 
		" AND ( line_bin<", coincidences[i, "frequency"]*1800+joint_bin_width, ")"))

	if(is.null(lines) || (dim(lines)[1]<1)) {
		ExtraStats[i, "line.f0"]<- 0
		ExtraStats[i, "line.comment"]<-""
		} else {
		ExtraStats[i, "line.f0"]<- lines[which.min(abs(lines[,'line_bin']-coincidences[i, "frequency"]*1800)), "line_bin"]/1800.0
		X<-aggregate(lines[,'strength'], lines[,'dataset', drop=FALSE], max, na.rm=TRUE)
		F<-X[,'x']>StrongLine
		F[is.na(F)]<-FALSE
		X<-X[F,,drop=FALSE]
		
		ExtraStats[i, "line.comment"]<-p(X[,'dataset'], "=", X[,'x']*LineStrengthToSNR, collapse=" ")
		}
	}
cat("\n")

# 	ExtraStats[i, "OC.H1L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set'], "'AND `skyband`=", coincidences[i, 'skyband']))[1,1]
# 
# 	ExtraStats[i, "OC.H1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set_H1'], "' AND `skyband`=", coincidences[i, 'skyband_H1']))[1,1]
# 
# 	ExtraStats[i, "OC.L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", coincidences[i, 'label'], "' AND `set`='", coincidences[i, 'set_L1'], "'AND `skyband`=", coincidences[i, 'skyband_L1']))[1,1]


FullData<-data.frame(ExtraStats, coincidences)

write.table(FullData, p(output.dir, "/coinc_snr_coincidences.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

write.table(FullData[ExtraStats$Comment=="",,drop=FALSE], p(output.dir, "/coinc_snr_coincidences_unmarked.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

write.table(FullData[ExtraStats$Comment=="" & pmin(coincidences$snr_L1, coincidences$snr_H1)>CoincMinSNRCutoff,,drop=FALSE], p(output.dir, "/coinc_snr_coincidences_largesnr_unmarked.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

ReducedData<-FullData[FullData$Primary>0, ,drop=FALSE]

F<-pmax(ReducedData[,"Max.SNR.H1"], ReducedData[,"Max.SNR.L1"], na.rm=TRUE) > ReducedData[,"SNR.H1L1"]
F[is.na(F)]<-FALSE
ReducedData[F,"Comment"]<-p(ReducedData[F,"Comment"], " single_ifo_max")

ReducedData<-ReducedData[order(ReducedData[,"SNR.H1L1"], pmin(ReducedData[,"SNR.H1"], ReducedData[,"SNR.L1"]), decreasing=TRUE),,drop=FALSE]

write.table(ReducedData, p(output.dir, "/coinc_snr_coincidences_reduced.csv"), row.names=FALSE, col.names=TRUE, sep="\t")

dag<-file("dag.followup", open="w")
for(i in 1:dim(ReducedData)[1]) {
	firstbin<-round(ReducedData[i, "f0"]*1800-250)
	cat(file=dag, "JOB D", i, " condor\n", sep="")
	cat(file=dag, "VARS D", i, " PID=\"", i, "\" RA=\"", ReducedData[i, "ra"], "\" DEC=\"", ReducedData[i, "dec"], "\" SPINDOWN=\"", ReducedData[i, "fdot"], "\" SPINDOWN_SWEEP_START=\"", ReducedData[i, "fdot"]-DagSpindownOffset, "\" F0=\"", ReducedData[i, "f0"], "\"  FIRSTBIN=\"", firstbin, "\" DATASET=\"", DagDataset(firstbin), "\"\n", sep="")
	}
close(dag)

start_png("coinc_snr_coincidences_f0fdot.start_png", width=600, height=600)
print(xyplot(spindown~frequency, coincidences))
dev.off()

start_png("coinc_snr_coincidences_sky.start_png", width=600, height=600)
print(xyplot(I(-dec)~I(-ra), coincidences))
dev.off()


start_png("coinc_max_m1_neg_hist.start_png", width=600, height=600)
print(histogram(~max_m1_neg, coincidences, nint=50))
dev.off()

start_png("coinc_min_m1_neg_hist.start_png", width=600, height=600)
print(histogram(~min_m1_neg, coincidences, nint=50))
dev.off()

start_png("coinc_max_m3_neg_hist.start_png", width=600, height=600)
print(histogram(~max_m3_neg, coincidences, nint=50))
dev.off()

start_png("coinc_min_m3_neg_hist.start_png", width=600, height=600)
print(histogram(~min_m3_neg, coincidences, nint=50))
dev.off()

start_png("coinc_max_m4_hist.start_png", width=600, height=600)
print(histogram(~max_m4, coincidences, nint=50))
dev.off()

start_png("coinc_min_m4_hist.start_png", width=600, height=600)
print(histogram(~min_m4, coincidences, nint=50))
dev.off()

start_png("coinc_snr_coincidences_sky_marked_by_mark.start_png", width=1024, height=1024)
print(xyplot(I(-dec)~I(-ra)|ExtraStats$Comment[ExtraStats$Comment!=""], coincidences[ExtraStats$Comment!="",,drop=FALSE], pch="+"))
dev.off()

start_png("coinc_snr_coincidences_f0snr_unmarked.start_png", width=600, height=600)
print(xyplot(snr~frequency, coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_f0minsnr_unmarked.start_png", width=600, height=600)
print(xyplot(pmin(snr_H1, snr_L1)~frequency, coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_f0fdot_unmarked.start_png", width=600, height=600)
print(xyplot(spindown~frequency, coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_sky_unmarked.start_png", width=600, height=600)
print(xyplot(I(-dec)~I(-ra), coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_sky_unmarked_by_f0.start_png", width=1024, height=1024)
print(xyplot(I(-dec)~I(-ra)|as.character(round(frequency*4)/4), coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_snr_unmarked_by_f0.start_png", width=1024, height=1024)
print(xyplot(snr_H1~snr_L1|as.character(round(frequency*4)/4), coincidences[ExtraStats$Comment=="",,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_snr_hist_unmarked.start_png", width=600, height=600)
print(histogram(~pmin(snr_H1, snr_L1), coincidences[ExtraStats$Comment=="",,drop=FALSE], nint=20))
dev.off()

start_png("coinc_snr_coincidences_sky_unmarked_by_f0_overlayed.start_png", width=1024, height=1024)
X<-FullData[ExtraStats$Comment=="",,drop=FALSE]
X[,"fgroup"]<-round(X$frequency*4)/4
plot(-X$ra, -X$dec, type="n")
Groups<-unique(X$fgroup)
cols<-rainbow(length(Groups), v=0.75)
for(i in 1:length(Groups)) {
	Y<-X[X$fgroup==Groups[[i]],,drop=FALSE]
	points(-Y$ra, -Y$dec, col=cols[[i]], pch=19)
	}
dev.off()

start_png("coinc_snr_coincidences_sky_unmarked_largesnr_by_f0.start_png", width=1024, height=1024)
print(xyplot(I(-dec)~I(-ra)|as.character(round(frequency*4)/4), coincidences[ExtraStats$Comment=="" & pmin(coincidences$snr_L1, coincidences$snr_H1)>CoincMinSNRCutoff,,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_snr_unmarked_largesnr_by_f0.start_png", width=1024, height=1024)
print(xyplot(snr_H1~snr_L1|as.character(round(frequency*4)/4), coincidences[ExtraStats$Comment=="" & pmin(coincidences$snr_L1, coincidences$snr_H1)>CoincMinSNRCutoff,,drop=FALSE]))
dev.off()

start_png("coinc_snr_coincidences_sky_unmarked_largesnr_by_f0_overlayed.start_png", width=1024, height=1024)
X<-FullData[ExtraStats$Comment==""  & pmin(coincidences$snr_L1, coincidences$snr_H1)>CoincMinSNRCutoff,,drop=FALSE]
X[,"fgroup"]<-round(X$frequency*4)/4
plot(-X$ra, -X$dec, type="n")
Groups<-unique(X$fgroup)
cols<-rainbow(length(Groups), v=0.75)
for(i in 1:length(Groups)) {
	Y<-X[X$fgroup==Groups[[i]],,drop=FALSE]
	points(-Y$ra, -Y$dec, col=cols[[i]], pch=19)
	}
dev.off()

close(LOG)


