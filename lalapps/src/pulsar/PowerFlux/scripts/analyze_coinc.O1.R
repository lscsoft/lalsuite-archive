#!/usr/bin/env Rscript

require("lattice")
require("RMySQL")
require("parallel")

p<-function(...) {
	return(paste(..., sep=""))
	}

SkipLabel<- -1
#SkipLabel<-165600000

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

create_index<-function() {
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ADD COLUMN (ECL_X DOUBLE, ECL_Y DOUBLE, F0INDEX INTEGER NOT NULL)"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " MODIFY COLUMN `set` VARCHAR(255)"))
	dbGetQuery(con, p("UPDATE ", CoincTableName, " SET ECL_X=cos(ra)*cos(`dec`), ECL_Y=sin(ra)*cos(`dec`)*", ecliptic_y_axis[2], "+sin(`dec`)*", ecliptic_y_axis[3], ", F0INDEX=ROUND(FREQUENCY/", FrequencyTolerance, ")"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ADD INDEX F0IDX(`label`, F0INDEX)"))
	dbGetQuery(con, p("ALTER TABLE ", CoincTableName, " ORDER BY `label`, F0INDEX"))
	}

waterfall<-function(X, dfrequency=WaterFallFrequencyTolerance, dspindown=WaterFallSpindownTolerance, decl=WaterFallRelEclTolerance, nominal.f0=400.0) {
	Y<-X[order(X[,"snr"], decreasing=TRUE),,drop=FALSE]
	Y[,"Group"]<- -1
	
	group<-Y[,"Group"]
	snr<-Y[,"snr"]
	frequency<-Y[,"frequency"]
	spindown<-Y[,"spindown"]
	ECL_X<-Y[,"ECL_X"]
	ECL_Y<-Y[,"ECL_Y"]
	
	for(i in 1:(dim(Y)[1])) {
		g<-Y[i, "Group"]
		if(g<0)g<-i
		
		F<-(snr < snr[i]) & (group<0) & (abs(frequency-frequency[i])< dfrequency) & (abs(spindown-spindown[i])<dspindown) & (sqrt((ECL_X-ECL_X[i])^2+(ECL_Y-ECL_Y[i])^2)<decl*(nominal.f0/frequency[i]))
		group[F]<-g
		}
	Y[,"Group"]<-group
	return(Y)
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

estimated_fstat_timebase<-function(h0, ra, dec, iota, psi, M.H1, M.L1,  twoFcutoff=60, H1.uptime=0.75, L1.uptime=0.75, sft.coherence=sft_length) {
	# LHO position
	gamma_LHO<-171.8*pi/180
	lambda_LHO<-46.45*pi/180

	A2_LHO<- ifo_A2(ra, dec, iota, psi, gamma_LHO, lambda_LHO)
	
	# LLO position
	gamma_LLO<-243.0*pi/180
	lambda_LLO<- 30.56*pi/180

	A2_LLO<- ifo_A2(ra, dec, iota, psi, gamma_LLO, lambda_LLO)

	Sh.H1<- 0.5*M.H1^2*sft.coherence
	Sh.L1<- 0.5*M.L1^2*sft.coherence

	timebase<- twoFcutoff/(A2_LHO*h0^2*H1.uptime/Sh.H1+A2_LLO*h0^2*L1.uptime/Sh.L1)
	return(timebase)
	}

output.dir<-p("report/coinc_", Segment)
	
dir.create(output.dir)

LOG<-file(p(output.dir, "/coinc.log"), open="w")
cat(file=LOG, "Using sft coherence length of", sft_length, "\n")


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

Lines_skyband<-dbGetQuery(con, p("SELECT skyband FROM `", BandInfoTableName, "` WHERE skyband_name='Lines' LIMIT 1"))[1,1]

strong_lines<-dbGetQuery(con, p("SELECT dataset, strength, line_bin, flag FROM `", LineInfoTableName, "` WHERE strength>", StrongLine))

start_png<-function(name, width=600, height=600) {
	png(p(output.dir, "/", name), width=width, height=height)
	}

start_png("strong_lines.png", width=1200, height=600)
print(xyplot(I(strength*LineStrengthToSNR)~I(line_bin*1.0/sft_length)|dataset, strong_lines, xlab="Frequency", ylab="SNR"))
dev.off()


high_bands<- dbGetQuery(con, p("SELECT F0INDEX, `label`, MAX(snr) as max_snr, COUNT(*) as `count` FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", Segment, "_all' AND snr>7 GROUP BY F0INDEX, `label`"))

cat("Found", dim(high_bands)[1], "indices with SNR>7\n")
#ExcludeBands<-c(108.9, 193.4, 393.4)
ExcludeBands<-c()
for(f0 in ExcludeBands) {
	F<- abs(high_bands[,"F0INDEX"]*FrequencyTolerance-f0)<0.25
	high_bands<-high_bands[!F,,drop=FALSE]
	}
cat(dim(high_bands)[1], "remain after excluding", p(ExcludeBands, collapse=" "), "bands\n")
high_bands<-high_bands[order(high_bands[,"label"], high_bands[,"F0INDEX"]),,drop=FALSE]

#high_bands[,"snr_cutoff_all"]<-pmax(high_bands[,"max_snr"]/2.0, 7.0)
#high_bands[,"snr_cutoff_ifo"]<-pmax(high_bands[,"max_snr"]/4.0, 5.0)

high_bands[,"snr_cutoff_all"]<-7
high_bands[,"snr_cutoff_ifo"]<-5

cat(file=LOG, sum(high_bands[,"snr_cutoff_all"]>7.0), "indices with elevated snr level for _all\n")
cat(file=LOG, sum(high_bands[,"snr_cutoff_ifo"]>5.0), "indices with elevated snr level for ifos\n")

start_png("coinc_max_snr.png", width=600, height=600)
print(xyplot(pmin(max_snr, 50)~I(F0INDEX*FrequencyTolerance), high_bands))
dev.off()

start_png("coinc_count_start.png", width=600, height=600)
print(xyplot(count~I(F0INDEX*FrequencyTolerance), high_bands))
dev.off()


Labels<-unique(as.integer(as.character(high_bands[, "label"])))

L<-list()
for(i in Labels) {
	L[[i+1]]<-NA
	}

dir.create(CoincData)


fd<-NULL
fd_name<-NULL

for(i in 1:dim(high_bands)[1]) {
	fidx<-high_bands[i, "F0INDEX"]
	label<-as.integer(as.character(high_bands[i, "label"]))

	if(label<SkipLabel) {
		cat("Skipping", i, "\n")
		next
		}

	a_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_all' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_all"], "  AND F0INDEX=", fidx, " AND `label`=", label, " AND frequency_bin>=20 AND frequency_bin<=480")

	b_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LHO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND F0INDEX IN (", fidx, ", ", fidx+1, ", ", fidx-1, ") AND `label`=", label, " AND frequency_bin>=20 AND frequency_bin<=480")

	c_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LLO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND F0INDEX IN (", fidx, ", ", fidx+1, ", ", fidx-1, ") AND `label`=", label, " AND frequency_bin>=20 AND frequency_bin<=480")

# 	a_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_all' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_all"], " AND F0INDEX=", fidx, " AND `label`=", label)
# 
# 	b_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LHO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND `label`=", label)
# 
# 	c_table<-p("SELECT * FROM `", CoincTableName, "` WHERE `set`='", Segment, "_LLO' AND kind='snr' AND snr>", high_bands[i, "snr_cutoff_ifo"], " AND `label`=", label)

#	coincidences<-dbGetQueryCoinc(con, p("SELECT * FROM (", a_table, ") a, (", b_table, ") b, (", c_table, ") c WHERE a.snr>b.snr AND a.snr>c.snr AND abs(a.frequency-b.frequency)<", FrequencyTolerance, " AND abs(a.frequency-c.frequency)<", FrequencyTolerance, " AND abs(a.`dec`-b.`dec`)<", DecTolerance, " AND abs(a.`dec`-c.`dec`)<", DecTolerance, " AND abs(a.spindown-b.spindown)<", SpindownTolerance, " AND abs(a.spindown-c.spindown)<", SpindownTolerance))
	
	coincidences<-dbGetQueryCoinc(con, p("SELECT * FROM (", a_table, ") a, (", b_table, ") b, (", c_table, ") c WHERE a.snr>b.snr AND a.snr>c.snr",
		" AND a.pi=b.pi AND a.pi=c.pi",
		" AND abs(a.frequency-b.frequency)<", FrequencyTolerance, 
		" AND abs(a.frequency-c.frequency)<", FrequencyTolerance, 
		" AND abs(a.spindown-b.spindown)<", SpindownTolerance, 
		" AND abs(a.spindown-c.spindown)<", SpindownTolerance))
		
		
#		" AND power(a.ecl_x-b.ecl_x, 2)+power(a.ecl_y-b.ecl_y, 2)<", LocationTolerance(fidx*FrequencyTolerance)^2,
#		" AND power(a.ecl_x-c.ecl_x, 2)+power(a.ecl_y-c.ecl_y, 2)<", LocationTolerance(fidx*FrequencyTolerance)^2, 

        if(length(coincidences)<1 || dim(coincidences)[1]<1)next

	#if(dim(coincidences)[1]<1)next
	cat(i, "\n")
	cat("Found", dim(coincidences)[1], "raw matches in band", fidx*FrequencyTolerance, "for instance", label, "\n")

#	names(coincidences)<-p(names(coincidences), rep(c("", "_H1", "_L1"), each=dim(coincidences)[2]/3))
#	FDist<- ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_H1"], coincidences[,"dec_H1"])<LocationTolerance(coincidences[,"frequency"]) &
#		ecliptic_dist(coincidences[,"ra"], coincidences[,"dec"], coincidences[,"ra_L1"], coincidences[,"dec_L1"])<LocationTolerance(coincidences[,"frequency"])

	#if(sum(is.na(FDist))>0)browser()

	#FDist[is.na(FDist)]<-FALSE
#	coincidences<-coincidences[FDist,,drop=FALSE]
	#if(dim(coincidences)[1]<1)next

	cat("Found", dim(coincidences)[1], "matches in band", fidx*FrequencyTolerance, "for instance", label, "\n")

	filename<-L[[label+1]]
	if(length(filename)==1 && is.na(filename)) {
	        filename<-p(CoincData, "/inj_", Segment, "_", label, ".csv")
		L[[label+1]]<-filename
		if(!is.null(fd))close(fd)
		fd<-file(filename, open="w")
		fd_name<-filename
		write.table(coincidences, fd, col.names=TRUE, row.names=FALSE, append=FALSE)
		}
		else {
		if(fd_name!=filename) {
			if(!is.null(fd))close(fd)
			fd<-file(filename, open="w")
			fd_name<-filename
			}
		write.table(coincidences, fd, col.names=FALSE, row.names=FALSE, append=TRUE)
		}

	}
if(!is.null(fd))close(fd)
fd<-NULL

Primary<-list()
for(i in Labels) {

	filename<-p(CoincData, "/inj_", Segment, "_", i, ".csv")

	A<-NULL
       	try({A<-read.table(filename, header=TRUE)}, silent=TRUE)

	if(is.null(A) || (dim(A)[1]<1))next

	cat(i-1, "\n")
	
	B<-waterfall(A)
	B<-B[B[,"Group"]<0,,drop=FALSE]
	
	B[,"Instance"]<-i
	B[,"SNR.H1"]<-B[,"snr_H1"]
	B[,"SNR.L1"]<-B[,"snr_L1"]
	B[,"SNR.H1L1"]<-B[,"snr"]
	B[,"f0"]<-B[,"frequency"]
	B[,"fdot"]<-B[,"spindown"]
	
	B[,"fdist"]<-pmax(abs(B$frequency-B$frequency_H1), abs(B$frequency-B$frequency_H1), na.rm=TRUE)
	B[,"sdist"]<-pmax(abs(B$spindown-B$spindown_H1), abs(B$spindown-B$spindown_H1), na.rm=TRUE)
	B[,"dist"]<-pmax(ecliptic_dist(B$ra, B$dec, B$ra_H1, B$dec_H1), ecliptic_dist(B$ra, B$dec, B$ra_L1, B$dec_L1), na.rm=TRUE)
	
	for(col in c("Max.SNR.H1L1", "Max.SNR.H1", "Max.SNR.L1", "H1L1.outliers", "H1.outliers", "L1.outliers", "H1.snr", "L1.snr", "H1L1.snr", "line.f0", "Group")) {
		B[,col]<-0
		}
	B[,"Primary"]<-1
	B[,"line.comment"]<-""
	B[,"Group"]<-1:(dim(B)[1])



	Primary[[length(Primary)+1]]<-B
	}

rm(L)

PrimaryMatches<-do.call(rbind, Primary)

PrimaryMatches<-PrimaryMatches[order(-PrimaryMatches[,"snr"], pmin(PrimaryMatches[,"snr_H1"], PrimaryMatches[,"snr_L1"])),,drop=FALSE]

PrimaryMatches[,"line_id"]<-p(PrimaryMatches[,"Instance"], "_", 1:dim(PrimaryMatches)[1])
cat(file=LOG, "Found", dim(PrimaryMatches)[1], "outliers\n")

if(!is.null(KnownLines)) {
	line.list<-read.table(KnownLines, header=TRUE, sep="\t")
	}


ExtraStats<-data.frame(Idx=0, Comment="", SNR.H1L1=PrimaryMatches$snr,
		SNR.H1=PrimaryMatches$snr_H1,
		SNR.L1=PrimaryMatches$snr_L1,
		Max.SNR.H1L1=0,
		Max.SNR.H1=0,
		Max.SNR.L1=0,
		f0=PrimaryMatches$frequency,
		fdot=PrimaryMatches$spindown,
		ra=PrimaryMatches$ra,
		dec=PrimaryMatches$dec,
		fdist=pmax(abs(PrimaryMatches$frequency-PrimaryMatches$frequency_H1), abs(PrimaryMatches$frequency-PrimaryMatches$frequency_L1), na.rm=TRUE),
		sdist=pmax(abs(PrimaryMatches$spindown-PrimaryMatches$spindown_H1), abs(PrimaryMatches$spindown-PrimaryMatches$spindown_L1), na.rm=TRUE),
		dist=pmax(ecliptic_dist(PrimaryMatches$ra, PrimaryMatches$dec, PrimaryMatches$ra_H1, PrimaryMatches$dec_H1), ecliptic_dist(PrimaryMatches$ra, PrimaryMatches$dec, PrimaryMatches$ra_L1, PrimaryMatches$dec_L1), na.rm=TRUE),
		line.f0=0,
		line.comment="",
		OC.H1L1=0,
		OC.H1=0,
		OC.L1=0,
		FStat60.timebase=estimated_fstat_timebase(h0=PrimaryMatches[, "snr"]*PrimaryMatches[, "S"], ra=PrimaryMatches[, "ra"], dec=PrimaryMatches[, "dec"], iota=PrimaryMatches[, "iota"], psi=PrimaryMatches[, "psi"], M.H1=PrimaryMatches[, "M_H1"], M.L1=PrimaryMatches[, "M_L1"], twoFcutoff=60)/(24*3600),
		Group=0,
		Primary=1
		)

ExtraStats[,"Idx"]<-1:dim(ExtraStats)[1]

ExtraStats[,"Comment"]<-""
ExtraStats[,"line.comment"]<-""
LineDF<-ExtraStats[, "f0"]*1e-4+abs(ExtraStats[, "fdot"])*0.5*run_timebase

for(i in 1:dim(line.list)[1]) {
	if(line.list[i, "type"]<1000) {
		F<- (ExtraStats[,"f0"]>line.list[i, "f0"]-line.list[i,"wl"]-LineDF) & (ExtraStats[,"f0"]<line.list[i, "f0"]+line.list[i,"wr"]+LineDF)
		} else {
		F<- (ExtraStats[,"f0"]>line.list[i, "f0"]-line.list[i,"wl"]) & (ExtraStats[,"f0"]<line.list[i, "f0"]+line.list[i,"wr"])
		}
        ExtraStats[F, "Comment"]<- p(ExtraStats[F, "Comment"], line.list[i, "ifo"], ":", line.list[i, "comment"], " ")
        }

#F<-PrimaryMatches[,"skyband"] %in% Lines_bands
#ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "Lines ")

F<- abs( (PrimaryMatches$frequency+1.25) %% 60.0 ) < 2.5
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "60 Hz ")

# F<- abs( (PrimaryMatches$frequency+0.05) %% 16.0 ) < 0.1
# ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "16 Hz ")

F<- PrimaryMatches$max_m1_neg>0.44
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m1 ")

F<- PrimaryMatches$min_m4<1.95
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "NonGauss_m4 ")

F<- PrimaryMatches$weight_loss_fraction_H1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_H1 ")

F<- PrimaryMatches$weight_loss_fraction_L1>0.05
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "weight_loss_L1 ")

F<- PrimaryMatches$skyband==Lines_skyband | PrimaryMatches$skyband_H1==Lines_skyband | PrimaryMatches$skyband_L1==Lines_skyband
ExtraStats[F,"Comment"]<-p(ExtraStats[F,"Comment"], "lines ")

# write.table(data.frame(ExtraStats[,"Idx",drop=FALSE], PrimaryMatches[,c("label", "F0INDEX", "frequency", "spindown", "ra", "dec", "set", "set_H1", "set_L1"), drop=FALSE]), p(output.dir, "/primary_matches0.csv"), col.names=TRUE, row.names=FALSE, sep="\t")
# 
# dbGetQuery(con, p("CREATE TEMPORARY TABLE tmp_matches0 (Idx INTEGER, `label` VARCHAR(300), F0INDEX INTEGER, frequency DOUBLE, spindown DOUBLE, ra DOUBLE, `dec` DOUBLE, `set` VARCHAR(255), `set_H1` VARCHAR(255), `set_L1` VARCHAR(255))"))
# 
# dbGetQuery(con, p("LOAD DATA LOCAL INFILE '", p(output.dir, "/primary_matches0.csv"), "' INTO TABLE tmp_matches0 FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (Idx, `label`, F0INDEX, frequency, spindown, ra, `dec`, `set`, `set_H1`, `set_L1`)"))
# 
# dbGetQuery(con, p("ALTER TABLE tmp_matches0 ADD INDEX (`Idx`, `label`, `F0INDEX`)"))
# 
# L<-list()
# for(i in (-MaxF0INDEXTolerance):MaxF0INDEXTolerance) {
# 	cat(i, "\n")
# 	L[[length(L)+1]]<-dbGetQuery(con, p("SELECT b.Idx, b.`label` as `label`, MAX(IF(a.`set`=b.`set`, a.snr, NULL)) as `Max.SNR.H1L1`, MAX(IF(a.`set`=b.`set_H1`, a.snr, NULL)) as `Max.SNR.H1`, MAX(IF(a.`set`=b.`set_L1`, a.snr, NULL)) as `Max.SNR.L1` FROM `", CoincTableName, "` as a RIGHT JOIN tmp_matches0 as b ON a.kind='snr' AND (sin(a.`dec`)*sin(b.`dec`)+ cos(a.`dec`)*cos(b.`dec`)*cos(a.ra-b.ra)>", cos(MaxLocationTolerance), ") AND a.label=b.label AND a.F0INDEX=b.F0INDEX+", i, " AND abs(a.frequency-b.frequency)<", MaxFrequencyTolerance, " AND abs(a.spindown-b.spindown)<", MaxSpindownTolerance, " GROUP BY b.`Idx`"))
# 	}
# 
# MaxSNRL<-do.call(rbind, L)
# MaxSNR<-aggregate(MaxSNRL[,c("Max.SNR.H1L1", "Max.SNR.H1", "Max.SNR.L1"),drop=FALSE], MaxSNRL[, c("Idx", "label"),drop=FALSE], max, na.rm=TRUE)


MaxSNR_H1L1<-rep(0.0, dim(PrimaryMatches)[1])
MaxSNR_H1<-rep(0.0, dim(PrimaryMatches)[1])
MaxSNR_L1<-rep(0.0, dim(PrimaryMatches)[1])
line.comment<-rep("", dim(PrimaryMatches)[1])
line.f0<-rep(0.0, dim(PrimaryMatches)[1])

for(i in 1:dim(ExtraStats)[1]) {

	#if(ExtraStats[i, "Comment"]!="" && ExtraStats[i, "Primary"]<1)next
	if(ExtraStats[i, "Primary"]<1)next
	if(i %% 100==0) cat(".")
	F0INDEX_list<- p((PrimaryMatches[i, "F0INDEX"]-MaxF0INDEXTolerance):(PrimaryMatches[i, "F0INDEX"]+MaxF0INDEXTolerance), collapse=", ")
	ra<-PrimaryMatches[i, "ra"]
	dec<-PrimaryMatches[i, "dec"]

# 	ExtraStats[i, "Max.SNR.H1L1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", PrimaryMatches[i, 'set'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND label='", PrimaryMatches[i, "label"], "' AND F0INDEX IN (", F0INDEX_list, ") AND abs(frequency-", PrimaryMatches[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", PrimaryMatches[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]
# 
# 	ExtraStats[i, "Max.SNR.H1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", PrimaryMatches[i, 'set_H1'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND label='", PrimaryMatches[i, "label"], "' AND F0INDEX IN (", F0INDEX_list, ") AND abs(frequency-", PrimaryMatches[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", PrimaryMatches[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]
# 
# 	ExtraStats[i, "Max.SNR.L1"]<-dbGetQuery(con, p("SELECT MAX(snr) FROM `", CoincTableName, "` WHERE kind='snr' AND `set`='", PrimaryMatches[i, 'set_L1'], "' AND (sin(`dec`)*", sin(dec), " + cos(`dec`)*", cos(dec), "*cos(ra-", ra, ")>", cos(MaxLocationTolerance), ") AND label='", PrimaryMatches[i, "label"], "' AND F0INDEX IN (", F0INDEX_list, ") AND abs(frequency-", PrimaryMatches[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", PrimaryMatches[i, "spindown"], ")<", MaxSpindownTolerance))[1,1]

	MaxSNR<-dbGetQuery(con, p("SELECT `set`, MAX(snr) as `Max.SNR` FROM `", CoincTableName, "` WHERE kind='snr' AND `set` IN ('", PrimaryMatches[i, 'set'], "', '", PrimaryMatches[i, 'set_H1'], "', '", PrimaryMatches[i, 'set_L1'], "') AND (sin(`dec`)*(", sin(dec), ") + cos(`dec`)*(", cos(dec), ")*cos(ra-(", ra, "))>", cos(MaxLocationTolerance), ") AND label='", PrimaryMatches[i, "label"], "' AND F0INDEX IN (", F0INDEX_list, ") AND abs(frequency-", PrimaryMatches[i, "frequency"], ")<", MaxFrequencyTolerance, " AND abs(spindown-", PrimaryMatches[i, "spindown"], ")<", MaxSpindownTolerance, " GROUP BY `set`"))
	
	#ExtraStats[i, "Max.SNR.H1L1"]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set'], "Max.SNR"][1]
	#ExtraStats[i, "Max.SNR.H1"]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set_H1'], "Max.SNR"][1]
	#ExtraStats[i, "Max.SNR.L1"]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set_L1'], "Max.SNR"][1]
	MaxSNR_H1L1[i]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set'], "Max.SNR"][1]
	MaxSNR_H1[i]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set_H1'], "Max.SNR"][1]
	MaxSNR_L1[i]<-MaxSNR[MaxSNR[,"set"]==PrimaryMatches[i, 'set_L1'], "Max.SNR"][1]
	
# 	ExtraStats[i, "OC.H1L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", PrimaryMatches[i, 'label'], "' AND `set`='", PrimaryMatches[i, 'set'], "'AND `skyband`!=", Lines_skyband))[1,1]
# 
# 	ExtraStats[i, "OC.H1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", PrimaryMatches[i, 'label'], "' AND `set`='", PrimaryMatches[i, 'set_H1'], "' AND `skyband`!=", Lines_skyband))[1,1]
# 
# 	ExtraStats[i, "OC.L1"]<-dbGetQuery(con, p("SELECT SUM(IF(`index`>=5, `count`, 0)) as OC3  FROM `", OutlierHistogramTableName, "` WHERE `label`='", PrimaryMatches[i, 'label'], "' AND `set`='", PrimaryMatches[i, 'set_L1'], "'AND `skyband`!=", Lines_skyband))[1,1]

	bin_width<-sft_length*(PrimaryMatches[i, "frequency"]*1e-4+abs(PrimaryMatches[i, "spindown"]*0.5*run_timebase))

	joint_bin_width<-max(bin_width, VeryStrongLineBinWidth)

	lines<-dbGetQuery(con, p("SELECT dataset, strength, line_bin, flag FROM `", LineInfoTableName, "` WHERE ((strength>", ExtraStats[i, "SNR.H1L1"]*OutlierSNRToLineStrength, " AND line_bin>", PrimaryMatches[i, "frequency"]*sft_length-bin_width, 
		" AND line_bin<", PrimaryMatches[i, "frequency"]*sft_length+bin_width, ") OR (strength>", VeryStrongLine, " AND line_bin>",  PrimaryMatches[i, "frequency"]*sft_length-VeryStrongLineBinWidth, 
		" AND line_bin<", PrimaryMatches[i, "frequency"]*sft_length+VeryStrongLineBinWidth, ")) AND line_bin>", PrimaryMatches[i, "frequency"]*sft_length-joint_bin_width, 
		" AND ( line_bin<", PrimaryMatches[i, "frequency"]*sft_length+joint_bin_width, ")"))

	if(is.null(lines) || (dim(lines)[1]<1)) {
		#ExtraStats[i, "line.f0"]<- 0
		#ExtraStats[i, "line.comment"]<-""
		line.f0[i]<- 0.0
		line.comment<- ""
		} else {
		X<-aggregate(lines[,'strength'], lines[,'dataset', drop=FALSE], max, na.rm=TRUE)
		F<-X[,'x']>StrongLine
		F[is.na(F)]<-FALSE
		X<-X[F,,drop=FALSE]
		
# 		ExtraStats[i, "line.comment"]<-p(X[,'dataset'], "=", X[,'x']*LineStrengthToSNR, collapse=" ")
# 		ExtraStats[i, "line.f0"]<- lines[which.min(abs(lines[,'line_bin']-PrimaryMatches[i, "frequency"]*sft_length)), "line_bin"]/sft_length
		line.comment[i]<-p(X[,'dataset'], "=", X[,'x']*LineStrengthToSNR, collapse=" ")
		line.f0[i]<- lines[which.min(abs(lines[,'line_bin']-PrimaryMatches[i, "frequency"]*sft_length)), "line_bin"]/sft_length
		}
	}
cat("\n")

ExtraStats[,"Max.SNR.H1L1"]<-MaxSNR_H1L1
ExtraStats[,"Max.SNR.H1"]<-MaxSNR_H1
ExtraStats[,"Max.SNR.L1"]<-MaxSNR_L1
ExtraStats[,"line.f0"]<-line.f0
ExtraStats[,"line.comment"]<-line.comment

FullData<-data.frame(ExtraStats, PrimaryMatches)

write.table(FullData, p(output.dir, "/primary_matches.csv"), sep="\t", col.names=TRUE, row.names=FALSE)


close(LOG)
