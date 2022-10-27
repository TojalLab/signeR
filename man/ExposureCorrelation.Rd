\name{ExposureCorrelation}
\docType{methods}
\alias{ExposureCorrelation}
\alias{ExposureCorrelation,SignExp,numeric-method}
\alias{ExposureCorrelation,matrix,numeric-method}

\title{Exposure correlation analysis (given a known sample feature)}

\description{
    ExposureCorrelation : Identify signatures which are significantly 
    correlated with a provided (numeric) sample feature.
}

\usage{
\S4method{ExposureCorrelation}{SignExp,numeric}(Exposures, feature,
    method="spearman", max_instances=200, cutoff_pvalue=0.05, quant=0.5,
    plot_to_file=FALSE, file="ExposureCorrelation_plot.pdf",
    colors=TRUE,...)
}

\arguments{
    \item{Exposures}{a SignExp object returned by signeR function or a matrix of 
    exposures (with signatures in rows and a column for each sample).}
    \item{feature}{numeric feature associated to each sample, such as
     age, weight or the expression of a gene.}
    \item{method}{a character string indicating which correlation 
    coefficient should be used for the test. Options are "pearson", 
    "kendall", or "spearman" (default).}
    \item{max_instances}{Maximum number of the exposure matrix instances to be 
    analyzed. If the number of available E instances is bigger then this 
    parameter, a subset of those will be randomly selected for analysis.}
    \item{cutoff_pvalue}{threshold for p-values quantile for signatures to be
    considered as showing significant correlation.}
    \item{quant}{the p-values quantile which, after log-transform, will 
    be used for selecting significantly correlated signatures. Deafult 
    is 0.5, which means the median p-value will be considered.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
    Default is FALSE.}
    \item{file}{Output file to export p-values boxplot and scatterplots 
    showing the correlations of exposures and the provided feature.}
    \item{colors}{Boolean variable, if TRUE p-values boxplots of 
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
    signatures <- signeR(M=muts, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}
# feature vector, with one value for each sample
my_feature <- rnorm(30,100,20)+signatures$SignExposures@Exp[1,,1]

Exp_corr <- ExposureCorrelation(signatures$SignExposures,feature=my_feature)

# see also
vignette(package="signeR")
}


