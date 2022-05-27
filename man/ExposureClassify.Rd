\name{ExposureClassify}
\docType{methods}
\alias{ExposureClassify}
\alias{ExposureClassify,SignExp,character-method}

\title{Classify samples by exposure levels}

\description{
    Assign unlabeled samples to previously defined groups.
}

\usage{
\S4method{ExposureClassify}{SignExp,character}(signexp_obj, labels, method="knn",
    k=3, weights=NA, plot_to_file=FALSE, file="Classification_barplot.pdf",
    colors=NA_character_, min_agree=0.75, ...)
}

\arguments{
    \item{signexp_obj}{A SignExp object returned by signeR function.}
    \item{labels}{Sample labels. Every sample labeled as NA will be classified
according to its mutational profile and the profiles of labeled samples.}
    \item{method}{Classification algorithm used. Default is k-Nearest
Neighbors (kNN). Any other algorithm may be used, as long as it is customized
to satisfy the following conditions:\cr
        Input: a matrix of labeled samples, with one sample per line and one
feature per column;
                a matrix of unlabeled samples to classify, with the same
structure;
                an array of labels, with one entry for each labeled
sample.\cr
        Output: an array of assigned labels, one for each unlabeled sample.
}
    \item{k}{Number of nearest neighbors considered for classification,
used only if method="kNN". Default is 3.}
    \item{weights}{Vector of weights applied to the signatures when performing
    classification. Default is NA, which leads all the signatures to have
    weight=1.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
Default is FALSE.}
    \item{file}{File that will be generated with classification graphic
output.}
    \item{colors}{Array of color names, one for each sample class. Colors will
be recycled if the length of this array is less than the number of classes.}
    \item{min_agree}{Minimum frequency of agreement among individual
classifications. Samples showing a frequency of agreement below this value are
considered as "undefined". Default is 0.75.}
    \item{...}{additional parameters for classification algorithm (defined by
"method" above).}
}

\value{
A list with the following items:
\item{class}{The assigned classes for each unlabeled sample.}
\item{freq}{Classification agreement for each unlabeled sample: the relative
frequency of assignment of each sample to the group specified in "class".}
\item{allfreqs}{Matrix with one column for each unlabeled sample and one row
for each group label. Contains the assignment frequencies of each sample to
each group.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    mut <- t(cbind(rmultinom(5, 500, 5:100),rmultinom(5, 500, 100:5)))
    signatures <- signeR(M=mut, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, main_eval=10, main_burn=10)
}

my_labels <- c("a","a","a","a",NA,"b","b","b","b",NA)
Class <- ExposureClassify(signatures$SignExposures, labels=my_labels)

# see also
vignette(package="signeR")
}

