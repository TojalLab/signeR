\name{ExposureSurvModel}
\docType{methods}
\alias{ExposureSurvModel}
\alias{ExposureSurvModel,SignExp,character-method}
\alias{ExposureSurvModel,SignExp,Surv-method}

\title{Exposure Cox model}

\description{
    ExposureSurvModel: Given survival data, fits a multivariate Cox proportional 
    hazards model to exposure data.
}

\usage{
\S4method{ExposureSurvModel}{SignExp,Surv}(Exposures, surv, addata,
    quant=0.5, cutoff_pvalue=0.05, cutoff_hr=NA, plot_to_file=FALSE, 
    file="ExposureSurvival_plot.pdf", colors=TRUE, ...)
}

\arguments{
    \item{Exposures}{A SignExp object returned by signeR function or a matrix of 
    exposures (with signatures in rows and a column for each sample).}
    \item{surv}{a Surv object from package survival or a matrix with columns 
    "time" and "status" (the last indicates whether 1:an event occurred or 
    0:there was a loss of follow up).}
    \item{addata}{a data frame  withh additional data (one sample per row) that 
    will be used in the Cox model along with exposure data.}
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
    \item{colors}{Boolean variable, if TRUE p-values boxplots of 
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
\item{Stats}{data frame containing hazard ratios and pvalues for signatures 
(one per line) on fitted Cox models.}
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
library(survival)
my_surv <- Surv(rnorm(30,730,100),sample(c(0:1), 30, replace = TRUE))

Exp_corr <- ExposureSurvModel(signatures$SignExposures, surv = my_surv)

# see also
vignette(package="signeR")
}


