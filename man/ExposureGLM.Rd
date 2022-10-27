\name{ExposureGLM}
\docType{methods}
\alias{ExposureGLM}
\alias{ExposureGLM,SignExp,numeric-method}
\alias{ExposureGLM,matrix,numeric-method}

\title{Exposure Generalized Linear Model}

\description{
     Fits a GLM to exposure data, with a given sample feature as the target 
     of the model.
}

\usage{
\S4method{ExposureGLM}{SignExp,numeric}(Exposures, feature, max_instances=200,
     cutoff_pvalue=0.05,quant=0.5, plot_to_file=FALSE,
     file="ExposureGLM_plot.pdf",colors=TRUE,...)
}
\arguments{
    \item{Exposures}{A SignExp object returned by signeR function or a matrix of 
    exposures (with signatures in rows and a column for each sample).}
    \item{feature}{numeric feature associated to each sample, such as
     age, weight or the expression of a gene.}
    \item{max_instances}{Maximum number of the exposure matrix instances to be 
    analyzed. If the number of available E instances is bigger then this 
    parameter, a subset of those will be randomly selected for analysis.}
    \item{cutoff_pvalue}{threshold for p-values quantile for signatures to be
    considered as significant on the model.}
    \item{quant}{p-values quantile used to evaluate if signatures are 
    significant. Default is 0.5, meaning that median p-values are adopted.}
    \item{plot_to_file}{Whether to save plots to the file parameter.
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
      whether it shows a significant contribution to the model.}
\item{Stats}{matrix of model statistics, with one line for each signature.}
\item{Pvalues}{vector of p-values used for significance evaluation.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    muts <- t(rmultinom(30, 500, rep(1:6,each=16))) #30 samples
    signatures <- signeR(M=muts, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, 
    main_eval=10, main_burn=10)
}

my_feature <- rnorm(30,100,20)+signatures$SignExposures@Exp[1,,1]
EGlm <- ExposureGLM(signatures$SignExposures, feature=my_feature)

# see also
vignette(package="signeR")
}

