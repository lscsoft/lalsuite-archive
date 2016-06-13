#
# Example of parameter file for single spindown analysis
#
prefix <- "stats/"
suffix <- ".H1.dat"
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

#
# Function to make hardcopy plots
#

start_plot<-function(name) {
	png(filename=paste(output_dir, "/", name, ".png", sep=""), width=800, height=800)
	}

#
# Alternative function to make PDF plots
#
