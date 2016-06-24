#
# this helper script will install signeR and its dependencies
# run it on the R console:
# source("http://rvalieris.github.io/signeR/install_signeR.R")
#

local({
	targetR <- "3.2"
	targetBioc <- "3.2"
	deps <- c(
		"R.oo", "VariantAnnotation", "RcppArmadillo",
		"nloptr", "NMF", "class", "tensorA", "devtools"
	)

	# check R version
	if(getRversion() < targetR) {
		stop(paste0("R version (",getRversion(),") ",
		"is below the recommended (",targetR,"), ",
		"see http://bioconductor.org/install/ for help."))
	}

	# load bioconductor
	source("http://bioconductor.org/biocLite.R")
	library(BiocInstaller)

	# check bioconductor version
	if(biocVersion() < targetBioc) {
		stop(paste0("Bionconductor version (",biocVersion(),") ",
		"is below the recommended (",targetBioc,"), ",
		"see http://bioconductor.org/install/ for help."))
	}

	# attempt to install missing packages
	for(p in deps[!(deps %in% installed.packages())]) {
		biocLite(p, suppressUpdates=TRUE)
		if(length(find.package(p, quiet=TRUE)) == 0) {
			stop(paste0("failed to install package: ",p,"."))
		}
	}

	# get latest release url
	latest_release <- jsonlite::fromJSON("https://api.github.com/repos/rvalieris/signeR/releases/latest")

	# install the latest signeR package
	devtools::install_url(latest_release$assets$browser_download_url)

	# update bioconductor
	biocLite()
})
