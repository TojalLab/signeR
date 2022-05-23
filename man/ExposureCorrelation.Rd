\name{ExposureCorrelation}
\docType{methods}
\alias{ExposureCorrelation}
\alias{ExposureCorrelation,SignExp,character-method}

\title{Exposure correlation analysis (given a known sample feature)}

\description{
    ExposureCorrelation : Identify signatures which are significantly 
    correlated with a provided (numeric) sample feature.
}

\usage{
\S4method{ExposureCorrelation}{SignExp,numeric}(signexp_obj, feature,
    method=spearman, quant=0.5, cutoff_pvalue=0.05,
    plot_to_file=FALSE, file="ExposureCorr.pdf",
    colored=TRUE, ...)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.}
    \item{feature}{numeric feature associated to each sample, such as
     age, weight or the expression of a gene.}
    \item{method}{a character string indicating which correlation 
    coefficient should be used for the test. Options are "pearson", 
    "kendall", or "spearman" (default).}
    \item{quant}{the p-values quantile which, after log-transform, will 
    be used for selecting significantly correlated signatures. Deafult 
    is 0.5, which means the median p-value will be considered.}
    \item{cutoff_pvalue}{threshold for p-values quantile for signatures to be
    considered as showing significant correlation.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
    Default is FALSE.}
    \item{file}{Output file to export p-values boxplot and scatterplots 
    showing the correlations of exposures and the provided feature.}
    \item{colored}{Boolean variable, if TRUE p-values boxplots of 
    significantly correlated signatures will be colored in green, 
    cutoff line will be colored in red and line segments showing the 
    transformed p-value quantile used for significance evaluation will be 
    colored in blue. Otherwise the plot will be black & white.}
    \item{...}{additional parameters for test algorithm defined by the method
    parameter.}
}

\value{
A list with the following items:
\item{Significance}{boolean array with one entry for each signature, indicating
whether it shows significant correlation with the provided feature.}
\item{Correlation_quantiles}{vector of correlation quantiles, with one entry 
for each signature.}
\item{Pvalues_quantiles}{vector of p-values quantiles used for significance 
evaluation.}
\item{Correlations}{matrix containing all computed correlations, with one row 
for each signature.}
\item{Pvalues}{matrix containing all computed p-values, with one row for each
signature.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    muts <- t(rmultinom(30, 500, rep(1:6,each=16))) #30 samples
    signatures <- signeR(M=muts, Mheader=FALSE, nsig=1, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}
# feature vector, with one value for each sample
my_feature <- rnorm(30,100,20)

Exp_corr <- ExposureCorrelation(signatures$SignExposures,feature=my_feature)

# see also
vignette(package="signeR")
}


