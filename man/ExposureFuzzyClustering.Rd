\name{FuzzyClustExp}
\docType{methods}
\alias{FuzzyClustExp}
\alias{FuzzyClustExp,SignExp,numeric-method}

\title{Fuzzy Clustering of exposure data}

\description{
    FuzzyClustExp : Performs fuzzy C-means clustering of samples, based on 
    exposures. The number of clusters is defined by optimizing the PBMF index 
    of obtained clustering.
}

\usage{
\S4method{FuzzyClustExp}{SignExp,numeric}(signexp_obj, Clim,
                method.dist="euclidean", method.clust="fcm", relative=FALSE, 
                m=2, plot_to_file=FALSE, file="FuzzyClustExp.pdf",colored=TRUE)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.}
    \item{Clim}{number of groups por a vector with mimum and maximum accepted 
    number of groups. The algorithm will maximize the PBMF-index whithin this 
    range.}
    \item{method.dist}{used distance metric}
    \item{method.clust}{clustering method. Either "fcm", default, for 
    fuzzy C-means or  "km" for k-means.}
    \item{relative}{Whether to normalize exposures of each sample so that they
    sum up to one. Default is FALSE, thus clustering samples by the absolute
    contributions of signatures to mutation counts. Otherwise, clustering will 
    be based on relative contributions.}
    \item{m}{Expoent used in PBMF-index}
    \item{plot_to_file}{Whether to save a heatmap of results to the file 
    parameter. Default is FALSE.}
    \item{file}{Output file to export a heatmap with the levels of pertinence 
    of samples to found groups.}
    \item{colored}{Whether plots will be in color or B&W. Default is TRUE.}
}

\value{
A list with the following items:
Meanfuzzy=Meanfuzzy,
                        AllFuzzy=Fuzzy[[1]],
                        Centroids=Fuzzy[[2]]
\item{Meanfuzzy}{Final clustering: mean levels of pertinence of samples to 
                 found groups.}
\item{AllFuzzy}{All levels of pertinence of samples to 
                 found groups in repeated runs of the clustering algorithm.}
\item{Centroids}{All centroids of found groups in repeated runs of the 
                 clustering algorithm.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    muts <- t(rmultinom(30, 500, rep(1:6,each=16))) #30 samples
    signatures <- signeR(M=muts, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}
# Limits to number of goups:
cl <- c(2,4)

FuzClust <- FuzzyClustExp(signatures$SignExposures, Clim = cl)

# see also
vignette(package="signeR")
}


