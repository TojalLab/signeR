\name{ExposureSurvival}
\docType{methods}
\alias{ExposureSurvival}
\alias{ExposureSurvival,SignExp,character-method}

\title{Exposure survival analysis}

\description{
    ExposureSurvival: Given survival data, identify signatures which are 
    significantly related with differences in hazards.
}

\usage{
\S4method{ExposureSurvival}{SignExp,Surv}(signexp_obj, surv,
    method=logrank, quant=0.5, cutoff_pvalue=0.05, cutoff_hr=NA,
    plot_to_file=FALSE, file="ExposureSurvival_plot.pdf",
    colors=TRUE, ...)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.}
    \item{surv}{a Surv object from package survival or a matrix with columns 
    "time" and "status" (the last indicates whether 1:an event occurred or 
    0:there was a loss of follow up).}
    \item{method}{a character string indicating which approach should be used 
    for the test. Options are "logrank" (default) or "cox" (fit a Cox 
    proportional hazards model to data).}
    \item{quant}{the quantile of p-values and hazard ratios which will 
    be used for selecting survival significant signatures. Deafult 
    is 0.5, which means the median p-value and hazard ratio will be considered.}
    \item{cutoff_pvalue}{threshold for p-values quantile for signatures to be
    considered as significant.}
    \item{cutoff_hr}{threshold for hazard ratio quantile for signatures to be
    considered as significant.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
    Default is FALSE.}
    \item{file}{Output file to export p-values boxplots and Kaplan-Meier curves.}
    \item{colored}{Boolean variable, if TRUE p-values boxplots of 
    significant signatures will be colored in green, 
    cutoff line will be colored in red and line segments showing the 
    transformed p-value quantile used for significance evaluation will be 
    colored in blue. Otherwise the plot will be black & white.}
    \item{...}{additional parameters for test algorithm defined by the method
    parameter.}
}

\value{
A list with the following items:
\item{Significance}{boolean array with one entry for each signature, indicating
whether its levels of exposure are significant to survival.}
\item{Correlation_quantiles}{vector of correlation quantiles, with one entry 
for each signature.}
\item{pvalues}{vector of p-values used for significance evaluation.}
\item{limits}{vector containing one cut value for the exposures of each 
signature, such that splitting the samples according to this value leads to 
maximal differences in survival among generated groups.}
\item{Groups}{matrix containing one line for each signature, definig a division 
of the samples in two groups according to its exposures, such that survival 
differences between the groups are maximal.}
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
my_surv <- Surv(rnorm(30,730,100),sample(c(0:1),30,replace=T))

Exp_corr <- ExposureSurvival(signatures$SignExposures, surv = my_surv)

# see also
vignette(package="signeR")
}


