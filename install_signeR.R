#
# this helper script will install signeR and its dependencies
#

local({
	target_R <- "3.2.2"
	if(getRversion() < target_R) {
		warning(paste0("R version detected (",getRversion(),") is below the recommended (",target_R,"), consider upgrading R."))
	}

	# dependencies
	deps <- c("R.oo", "VariantAnnotation", "RcppArmadillo", "nloptr", "NMF", "class", "tensorA", "devtools")

	# load bioconductor
	source("http://bioconductor.org/biocLite.R")

	# attempt to install missing packages
	for(p in deps[!(deps %in% installed.packages())]) {
		biocLite(p, suppressUpdates=TRUE)
		if(length(find.package(p, quiet=TRUE)) == 0) { # installation failed
			stop(paste0("failed to install package: ",p,"."))
		}
	}

	# install the latest signeR package
	devtools::install_url("https://github.com/rvalieris/signeR/releases/download/v0.5.0/signeR_0.5.0.tar.gz")

	# update bioconductor
	biocLite()
})
