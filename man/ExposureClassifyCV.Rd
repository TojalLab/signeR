\name{ExposureClassifyCV}
\docType{methods}
\alias{ExposureClassifyCV}
\alias{ExposureClassifyCV,SignExp,character-method}
\alias{ExposureClassifyCV,ANY,character-method}

\title{k-fold cross-validation of sample classifcation by exposure levels}

\description{
    Splits labeled samples in k groups (deafult k=8), keeping the proportion of 
    classes stable among groups. Classify samples in each group according to the 
    k-1 remaining ones. Gather results and evaluate global classification 
    performance.
}

\usage{
\S4method{ExposureClassifyCV}{SignExp,character}(signexp_obj, labels, method="knn",
    max_instances=200, k=3, weights=NA, plot_to_file=FALSE, 
    file="Classification_CV_barplot.pdf", colors=NA_character_, 
    min_agree=0.75, fold=8, ...)
}

\arguments{
    \item{signexp_obj}{A SignExp object returned by signeR function.}
    \item{labels}{Sample labels. Unlabeled samples (NA labels) will be ignored.}
    \item{method}{Classification algorithm used. Default is k-Nearest
Neighbors (kNN). Any other algorithm may be used, as long as it is customized
to satisfy the following conditions:\cr
        Input: a matrix of labeled samples, with one sample per line and one
feature per column;
           a matrix of unlabeled samples to classify, with the same structure;
          an array of labels, with one entry for each labeled sample.\cr
        Output: an array of assigned labels, one for each unlabeled sample.}
    \item{max_instances}{Maximum number of the exposure matrix instances to be 
    analyzed. If the number of available E instances is bigger than this 
    parameter, a subset of those will be randomly selected for analysis.}
    \item{k}{Number of nearest neighbors considered for classification,
used only if method="kNN". Default is 3.}
    \item{weights}{Vector of weights applied to the signatures when performing
    classification. Default is NA, which leads all the signatures to have
    weight=1.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
Default is FALSE.}
    \item{file}{File that will be generated with cross validation graphic
output.}
    \item{colors}{Array of color names, one for each sample class. Colors will
be recycled if the length of this array is less than the number of classes.}
    \item{min_agree}{Minimum frequency of agreement among individual
classifications. Samples showing a frequency of agreement below this value are
considered as "undefined". Default is 0.75.}
    \item{fold}{Number of subsets in which labeled samples will be split}
    \item{...}{additional parameters for classification algorithm (defined by
"method" above).}
}

\value{
A list with the following items:
\item{confusion_matrix}{Contingency table of attributed sample classes against 
original labels.}
\item{class}{The assigned classes for each sample.}
\item{freq}{Classification agreement for each sample: the relative
frequency of assignment of each sample to the group specified in "class".}
\item{allfreqs}{Matrix with one column for each sample and one row
for each class label. Contains the assignment frequencies of each sample to
each class.}
\item{probs}{As above, a matrix with samples in columns and class 
labels in rows. Contains the average probability, among repeated exposure 
classifications, of each sample belonging to each class.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    mut <- t(cbind(rmultinom(5, 500, 5:100),rmultinom(5, 500, 100:5)))
    rownames(mut)<-paste("sample",letters[1:10],sep="_")
    signatures <- signeR(M=mut, Mheader=TRUE, nsig=2, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}

my_labels <- c("a","a","a","a","a","b","b","b","b","b")
ClassCV <- ExposureClassifyCV(signatures$SignExposures, labels=my_labels,fold=5)

# see also
vignette(package="signeR")
}

