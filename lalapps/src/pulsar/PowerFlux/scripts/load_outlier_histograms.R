require("lattice")
require("RMySQL")


source("params.R")

con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

p<-function(...) {
	return(paste(..., sep=""))
	}

fn<-p(prefix, "outlier_histogram", suffix)
fnout<-p(fn)

NAStrings<-c("NA", "NaN", "NAN")

vn<-function(field, i, position) {
	return(data[,p(var_prefix, field, ".", i, ".", position)])
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}



cat("Loading data from ", fn, "\n", sep="")
# Cheat - get a small sample first and use it to find out types and such.
header<-read.table(pipe(p("head --lines=10001 ", fn)), header=TRUE, na.strings=NAStrings, sep="\t")
cat("Found header with", dim(header)[2], "columns\n")
#colnames(header)<-unlist(col_names[1, 2:dim(col_names)[2]])

Types<-lapply(header, class)
Types<-lapply(Types, function(x)gsub("logical", "numeric", x))


CreateQuery<-p("CREATE TABLE ", OutlierHistogramTableName, "(Line INTEGER AUTO_INCREMENT")
LoadQuery<-p("LOAD DATA LOCAL INFILE '", fnout, "' INTO TABLE ", OutlierHistogramTableName, " FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (")

# "

for(i in 1:length(Types)) {
	Col<-names(Types)[i]

	Name<-Col

	Name<-gsub("\\.", "_", Name)

	CreateQuery<-p(CreateQuery, ", `", Name, "` ", switch(Types[[i]], integer="INTEGER", factor="VARCHAR(10000)", numeric="DOUBLE", "NA"))
	LoadQuery<-p(LoadQuery, ifelse(i==1, "`", ", `"), Name, "`")
	}
CreateQuery<-p(CreateQuery, p(", PRIMARY KEY (Line))"))
cat(CreateQuery, "\n")
LoadQuery<-p(LoadQuery, ")")

#cat("Preprocessing input data (NaN -> NULL)\n")
#system(p("sed 's/[nN][aA][nN]/NULL/g' stats/stat.dat > ", fnout))

cat("Connecting to the database\n")
con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

cat("Dropping table", OutlierHistogramTableName, "\n")
try(dbGetQuery(con, p("DROP TABLE ", OutlierHistogramTableName)), silent=TRUE)

cat("Declaring table", OutlierHistogramTableName, "\n")
dbGetQuery(con, CreateQuery)

cat("Loading table", OutlierHistogramTableName, "\n")
dbGetQuery(con, LoadQuery)


cat("Adding index to table", OutlierHistogramTableName, "\n")
dbGetQuery(con, p("ALTER TABLE `", OutlierHistogramTableName, "` ADD INDEX (`label`)"))


cat("Warnings:\n")
print(dbGetQuery(con, "SHOW WARNINGS"))

#Types<-lapply(Types, function(x) {
#	if(x=="factor")return("character")
#	return(x)
#	})
#FieldsUnused<-setdiff(names(Types), FieldsUsed)
#Types[FieldsUnused]<-NULL



