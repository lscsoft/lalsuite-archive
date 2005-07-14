prefix <- "stats/"
suffix <- ".H2.dat"
var_prefix <- ""
output_dir <- "report/"
Bands <- c(0,1,2,6,7,8)
DxCap <- 50.0
# 0.07 is nice cutoff for KS test 
ksLarge <- 0.07
# 7 is rather large detection strength
DxLarge <- 7
# High value for residuals:
ResLarge <- 1.5


h0ULtitle<-"h0 UL for S4 H2"
circULtitle <- "Circular UL for S4 H2"
DXtitle<- "Maximum detection strength for S4 H2"

#
# Function to make hardcopy plots
#

start_plot<-function(name) {
	png(filename=paste(output_dir, "/", name, ".png", sep=""), width=800, height=800)
	}

