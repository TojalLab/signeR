\name{HClustExp}
\docType{methods}
\alias{HClustExp}
\alias{HClustExp,SignExp-method}
\alias{HClustExp,SignExp,numeric-method}

\title{Hierarchical Clustering of exposure data}

\description{
    HClustExp: Performs hierarchical clustering of samples, based on 
    exposures.
    }

\usage{
\S4method{HClustExp}{SignExp,numeric}(signexp_obj, Med_exp=NA,
        max_instances=200, method.dist="euclidean", method.hclust="average", 
        use.cor=FALSE, relative=FALSE, plot_to_file=FALSE, 
        file="HClustExp_dendrogram.pdf", colored=TRUE)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.}
    \item{Med_exp}{optional matrix with (median) exposures.}
    \item{max_instances}{Maximum number of the exposure matrix instances to be 
    analyzed. If the number of available E instances is bigger then this 
    parameter, a subset of those will be randomly selected for analysis.}
    \item{method.dist}{used distance metric}
    \item{method.hclust}{clustering method.}
    \item{use.cor}{used in pv.distance}
    \item{relative}{Whether to normalize exposures of each sample so that they
    sum up to one. Default is FALSE, thus clustering samples by the absolute
    contributions of signatures to mutation counts. Otherwise, clustering will 
    be based on relative contributions.}
    \item{plot_to_file}{Whether to save a heatmap of results to the file 
    parameter. Default is FALSE.}
    \item{file}{Output file to export a heatmap with the levels of pertinence 
    of samples to found groups.}
    \item{colored}{Whether plots will be in color or B&W. Default is TRUE.}

}

\value{
A pvclust object, as described in package pvclust.
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    muts <- t(rmultinom(30, 500, rep(1:6,each=16))) #30 samples
    signatures <- signeR(M=muts, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}

HClust <- HClustExp(signatures$SignExposures)

# see also
vignette(package="signeR")
}


