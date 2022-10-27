.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath(
    prefix = "assets",
    directoryPath = system.file(
      "www/assets",
      package = "signeR"
    )
  )
}

.onUnload <- function(libname, pkgname) {
  shiny::removeResourcePath("assets")
}

.get_cache <- function() {
    cache <- tools::R_user_dir("signerflow", which = "cache")
    BiocFileCache::BiocFileCache(cache, ask = FALSE)
}

download_data_file <- function(type, tumor, verbose = FALSE) {
    rootURL <- "https://gitlab.com/lbcb/signer-data/-/raw/main/tcga/"

    url <- paste0(
      rootURL, type, "/", tumor, "-signatures.RData"
    )

    bfc <- .get_cache()
    rid <- bfcquery(bfc, paste0(type, "-", tumor), "rname")$rid
    if (!length(rid)) {
     if (verbose)
         message("Downloading TCGA data")
     rid <- names(bfcadd(
        bfc, paste0(type, "-", tumor), url
       ))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid)))
    bfcdownload(bfc, rid, ask = FALSE)

    data <- bfcrpath(bfc, rids = rid)
    return(data)
}

download_opp_file <- function(build, verbose = FALSE) {
    rootURL <- "https://gitlab.com/lbcb/signer-data/-/raw/main/"

    url <- paste0(
      rootURL, "/", build, "/opp.refgene.",build,".txt"
    )

    bfc <- .get_cache()
    rid <- bfcquery(bfc, build, "rname")$rid
    if (!length(rid)) {
     if (verbose)
         message("Downloading genome opportunity data")
     rid <- names(bfcadd(
        bfc, build, url
       ))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid)))
    bfcdownload(bfc, rid, ask = FALSE)

    data <- bfcrpath(bfc, rids = rid)
    return(data)
}

download_clinical_file <- function(tumor, verbose = FALSE) {
    rootURL <- "https://gitlab.com/lbcb/signer-data/-/raw/main/tcga/"

    url <- paste0(
      rootURL, "clinical/", tumor, ".tsv.gz"
    )

    bfc <- .get_cache()
    rid <- bfcquery(bfc, paste0(tumor,"-clinical"), "rname")$rid
    if (!length(rid)) {
     if (verbose)
         message("Downloading TCGA clinical data")
     rid <- names(bfcadd(
        bfc, paste0(tumor,"-clinical"), url
       ))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid)))
    bfcdownload(bfc, rid, ask = FALSE)

    data <- bfcrpath(bfc, rids = rid)
    return(data)
}