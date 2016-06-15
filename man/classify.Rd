\name{Classify}
\alias{Classify.SignExp}
\alias{Classify}

\title{Classify unknown samples}

\description{
    Classify: Assign unknown samples to previously defined groups.
}

\usage{
\S3method{Classify}{SignExp}(this, labels, method="knn", k=3, 
    plotfile="Classification_barplot.pdf", colors=NA, ...)
}

\arguments{
    \item{this}{a SignExp object returned by signeR function.}
    \item{labels}{sample labels. Every sample labeled as NA will be classified according to its mutational profile and the profiles of labeled samples.}
    \item{method}{Classification algorithm used. Default is k-Nearest Neighbors (kNN). 
Any other algorithm may be used, as long as it is customized to satisfy the following conditions:\cr
        Input: a matrix of labeled samples, with one sample per line and one feature per column;
                  a matrix of unlabeled samples to classify, with the same structure;
                  an array of labels, with one entry for each labeled sample.\cr
        Output: an array of assigned labels, one for each unlabeled sample.
}
    \item{k}{number of nearest neighbors considered for classification, used only if method="kNN". Default is 3.}
    \item{plotfile}{file that will be generated with classification graphic output.}
    \item{colors}{Array of color names, one for each sample class. Colors will be recycled if the length of this array is less than the number of classes.}
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




