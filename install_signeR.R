
# load bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(ask=FALSE)

if(!require("devtools")) { biocLite("devtools", ask=FALSE) } # attempt to install devtools

if(require("devtools")) {
	# devtools loaded
	devtools::install_github("rvalieris/signeR")
} else {
	# something wrong with devtools, attempt manual installation
	biocLite(c("BiocGenerics","Biostrings","BSgenome","class",
		"GenomicRanges","nloptr","NMF","plot3D","R.methodsS3",
		"R.oo","tensorA","VariantAnnotation","Rcpp","RcppArmadillo"))
	install.packages("https://github.com/rvalieris/signeR/releases/download/0.5.0/signeR_0.5.0.tar.gz")
}

