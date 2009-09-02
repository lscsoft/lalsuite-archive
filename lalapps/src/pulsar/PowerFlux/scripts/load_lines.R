require("lattice")
require("RMySQL")


source("params.R")

con<-dbConnect(dbDriver("MySQL"), user="volodya", password="", dbname="volodya")

p<-function(...) {
	return(paste(..., sep=""))
	}

fn<-p(prefix, "lines", suffix)
fnout<-p(fn)

NAStrings<-c("NA", "NaN", "NAN")

vn<-function(field, i, position) {
	return(data[,p(var_prefix, field, ".", i, ".", position)])
	}

ofn<-function(name) {
	return(file(paste(output_dir, "/", name, sep=""), open="wt"))
	}
cat("Loading firstbin indices\n")
firstbin<-read.table(p(prefix, "firstbin", suffix), header=FALSE)
names(firstbin)<-c("path", "", "firstbin")
firstbin[,"instance"]<-gsub("output/([^/]*)/.*", "\\1", as.character(firstbin[,"path"]))
write.table(firstbin[,c("instance", "firstbin"),drop=FALSE], "firstbin.txt.tmp", sep="\t", row.names=FALSE, col.names=FALSE)
instance.type<-class(read.table("firstbin.txt.tmp", sep="\t", header=FALSE)[,1])



cat("Loading data from ", fn, "\n", sep="")
# Cheat - get a small sample first and use it to find out types and such.
header<-read.table(pipe(p("head --lines=10001 ", fn)), header=TRUE, na.strings=NAStrings, sep="\t")
cat("Found header with", dim(header)[2], "columns\n")
#colnames(header)<-unlist(col_names[1, 2:dim(col_names)[2]])

Types<-lapply(header, class)
Types<-lapply(Types, function(x)gsub("logical", "numeric", x))


CreateQuery<-p("CREATE TABLE ", LineInfoTableName, "(Line INTEGER AUTO_INCREMENT")
LoadQuery<-p("LOAD DATA LOCAL INFILE '", fnout, "' INTO TABLE ", LineInfoTableName, " FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' LINES TERMINATED BY '\n' IGNORE 1 LINES (")

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

cat("Dropping table", LineInfoTableName, "\n")
try(dbGetQuery(con, p("DROP TABLE ", LineInfoTableName)), silent=TRUE)

cat("Declaring table", LineInfoTableName, "\n")
dbGetQuery(con, CreateQuery)

cat("Loading table", LineInfoTableName, "\n")
dbGetQuery(con, LoadQuery)

cat("Uploading firstbin table\n")
dbGetQuery(con, p("CREATE TEMPORARY TABLE firstbin_tmp (instance ", switch(instance.type, integer="INTEGER", factor="VARCHAR(10000)", numeric="DOUBLE", "NA") ,", firstbin INTEGER, INDEX I(instance))"))
dbGetQuery(con, "LOAD DATA LOCAL INFILE 'firstbin.txt.tmp' INTO TABLE firstbin_tmp FIELDS TERMINATED BY '\t' OPTIONALLY ENCLOSED BY '\"' LINES TERMINATED BY '\n' (instance, firstbin)")
#"

cat("Computing line bins\n")
dbGetQuery(con, p("ALTER TABLE `", LineInfoTableName, "` ADD COLUMN line_bin INTEGER"))
dbGetQuery(con, p("UPDATE `", LineInfoTableName, "` a JOIN firstbin_tmp b USING (instance) SET line_bin=b.firstbin+bin"))


cat("Adding index to table", LineInfoTableName, "\n")
dbGetQuery(con, p("ALTER TABLE `", LineInfoTableName, "` ADD INDEX (line_bin)"))


cat("Warnings:\n")
print(dbGetQuery(con, "SHOW WARNINGS"))

#Types<-lapply(Types, function(x) {
#	if(x=="factor")return("character")
#	return(x)
#	})
#FieldsUnused<-setdiff(names(Types), FieldsUsed)
#Types[FieldsUnused]<-NULL



