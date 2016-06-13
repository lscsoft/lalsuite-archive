source("params.R")

require("lattice")
require("RMySQL")

p<-function(...) {
	return(paste(..., sep=""))
	}

fn<-p(prefix, "stat", suffix)
fnout<-p(fn, ".sed")

NAStrings<-c("NA", "NaN", "NAN")

vn<-function(field, i, position) {
	return(data[,p(var_prefix, field, ".", i, ".", position)])
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}


cat("Loading data from ", fn, "\n", sep="")
# Cheat - get a small sample first and use it to find out types and such.
header<-read.table(pipe(p("head --lines=10001 ", fn)), header=TRUE, na.strings=NAStrings, sep=" ")
cat("Found header with", dim(header)[2], "columns\n")

Types<-lapply(header, class)
Types<-lapply(Types, function(x)gsub("logical", "numeric", x))

BandNames<-grep("grid_points.*2$", colnames(header), value=TRUE)
NBands<- length(BandNames)

BandData<-data.frame(i=1:(NBands/2), idx=NA, Name=NA, ref_idx=NA, ref_Name=NA)
k<-1
for(i in 1:NBands) {
	if(regexpr("_ref$", as.character(header[1, p("grid_points.", i-1, ".2")]))>=0)next
	BandData[k, 'idx']<-i-1
	BandData[k, 'Name']<- as.character(header[1, p("grid_points.", i-1, ".2")])
	for(m in 1:NBands) {
		if(p(header[1, p("grid_points.", i-1, ".2")], "_ref")==header[1, p("grid_points.", m-1, ".2")]) {
			BandData[k, 'ref_idx']<-m-1
			BandData[k, 'ref_Name']<- as.character(header[1, p("grid_points.", m-1, ".2")])
			break
			}
		}
	k<-k+1
	}

FieldsUsed<-c(
	band="band.3",
	spindown="spindown.2",
	hist_residuals_min="hist_residuals.3",
	hist_residuals_max="hist_residuals.4",
	median="median.3",
	sky_max_dx_ra="max_dx.1",
	sky_max_dx_dec="max_dx.2",
	polarization="max_dx.3",
	sky_max_dx_snr="max_dx.4",
	sky_max_dx_strain="max_dx.5",
	sky_max_dx_freq="max_dx.7",
	masked="masked.2",
	cputime="cputime.2"
	)

for(skyband in 0:(NBands-1)) {
	L<-c(	max_dx_sky_band=p("max_dx_band.", skyband, ".1"), 
		max_dx_sky_band_name=p("max_dx_band.", skyband, ".2"), 
		max_dx_snr=p("max_dx_band.", skyband, ".3"), 
		max_dx_pol=p("max_dx_band.", skyband, ".4"), 
		max_dx_freq=p("max_dx_band.", skyband, ".5"), 
		max_dx_ra=p("max_dx_band.", skyband, ".6"), 
		max_dx_dec=p("max_dx_band.", skyband, ".7"), 
		max_dx_point_index=p("max_dx_band.", skyband, ".8"), 
		grid_points_sky_band_name=p("grid_points.", skyband, ".2"), 
		p("ks_hist.", skyband, ".3"), 
		p("max_band.", skyband, ".5"),
		p("max_band.", skyband, ".6"),
		p("max_high_ul_band.", skyband, ".2"),
		p("max_ratio.", skyband, ".3")
		)
	F<-names(L)==""
	names(L)<-p(names(L), "_", skyband)
	L[F]<-""
	FieldsUsed<-c(FieldsUsed, L)
	}

CreateQuery<-p("CREATE TABLE ", DataSet, "(Line INTEGER AUTO_INCREMENT")
LoadQuery<-p("LOAD DATA LOCAL INFILE '", fnout, "' INTO TABLE ", DataSet, " FIELDS TERMINATED BY ' ' OPTIONALLY ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (")

# "

for(i in 1:length(Types)) {
	Col<-names(Types)[i]

	F<- FieldsUsed %in% Col
	if(sum(F)>0) {
		L<-FieldsUsed[F]
		Name <- names(L)[1]
		#cat(Col, Name, "\n")
		if(Name=="")Name<-Col
		} else {
		Name<-Col
		}

	Name<-gsub("\\.", "_", Name)

	CreateQuery<-p(CreateQuery, ", ", Name, " ", switch(Types[[i]], integer="INTEGER", factor="TEXT", numeric="DOUBLE", "NA"))
	LoadQuery<-p(LoadQuery, ifelse(i==1, "", ", "), Name)
	}
CreateQuery<-p(CreateQuery, p(", PRIMARY KEY (Line))"))
cat(CreateQuery, "\n")
LoadQuery<-p(LoadQuery, ")")

cat("Preprocessing input data (NaN -> NULL)\n")
system(p("sed 's/[nN][aA][nN]/NULL/g' stats/stat.dat > ", fnout))

cat("Connecting to the database\n")
con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

cat("Dropping table", DataSet, "\n")
try(dbGetQuery(con, p("DROP TABLE ", DataSet)), silent=TRUE)

cat("Declaring table", DataSet, "\n")
dbGetQuery(con, CreateQuery)

cat("Loading table", DataSet, "\n")
dbGetQuery(con, LoadQuery)

cat("Warnings:\n")
print(dbGetQuery(con, "SHOW WARNINGS"))

#Types<-lapply(Types, function(x) {
#	if(x=="factor")return("character")
#	return(x)
#	})
#FieldsUnused<-setdiff(names(Types), FieldsUsed)
#Types[FieldsUnused]<-NULL



