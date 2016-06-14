\name{Classify}
\alias{Classify.SignExp}
\alias{Classify}

\title{Classify unknown samples}

\description{
    Classify: Assign unknown samples to previously defined groups.
}

\usage{
\S3method{Classify}{SignExp}(this, labels, method="knn", k=3, plotfile="Classification_barplot.pdf", ...)
}

\arguments{
    \item{this}{a SignExp object returned by signeR function.}
    \item{labels}{sample labels. Every sample labeled as NA will be classified according to its mutational profile and the profiles of labeled samples.}
    \item{method}{classification algorithm used. Default is k-Nearest Neighbors (kNN). }
    \item{k}{number of nearest neighbors considered for classification, used only if method="kNN". Default is 3.}
    \item{plotfile}{file that will be generated with classification graphic output.}
    \item{...}{additional parameters for classification algorithm (defined by "method" above).}
}

\value{
A list with the following items:
\item{class}{The assigned classes for each unlabeled sample.}
\item{freq}{Classification agreement for each unlabeled sample: the relative frequency of assignment of each sample to the group specified in "class".}
\item{allfreqs}{Matrix with one column for each unlabeled sample and one row for each group label. Contains the assignment frequencies of each sample to each group.}
}

\examples{
# assuming signatures is the return value of signeR()

\dontrun{
    my_labels <- c("a","a",NA,"b","b",NA)
    Class <- Classify(signatures$SignExposures, labels=my_labels, plotfile="Sample_classification.pdf")
}
# see also
vignette(package="signeR")
}




