\name{DiffExp}
\alias{DiffExp.SignExp}
\alias{DiffExp}

\title{Differential Exposure Analysis}

\description{
    DiffExp : Identify signatures with significantly different activities among sample groups.
}

\usage{
\S3method{DiffExp}{SignExp}(this, labels, method='kw', contrast="all", quant=0.5, cutoff=0.05, plotfile="Diffexp_boxplot.pdf", pal="color", ...)
}

\arguments{
    \item{this}{a SignExp object returned by signeR function.}
    \item{labels}{sample labels used to define sample groups.}
    \item{method}{algorithm used to compare each signature exposures among sample groups. Default is "kw", which leads to the use of Kruskal-Wallis Rank Sum Test.}
    \item{contrast}{defines which sample groups will be considered in the analysis. Default is "all", which leads the algorithm to evaluates the null hypothesis of exposure levels beeing constant in all groups. Instead, if this parameter contains a list of group labels, the algorithm will evaluate the null hypothesis of exposure levels beeing constant among those groups.}
    \item{quant}{the p-values quantile which, after log-transform, will be used as DES (Differental Exposure Score). Deafult is 0.5, which means the median log-transformed p-value will be considered as DES.}
    \item{cutoff}{threshold for p-values quantile for signatures to be considered as showing differential exposure.}
    \item{plotfile}{Output file to export p-values boxplot.}
    \item{pal}{Color palette.}
    \item{...}{additional parameters for test algorithm defined by the method parameter.}
}

\value{
A list with the following items:
\item{Pvquant}{boolean array with one entry for each signature, indicating whether it shows differential exposure.}
\item{Pvalues}{matrix containing all computed p-values, with one row for each signature.}
\item{MostExposed}{for each differentially exposed signature, this array contains the label of the group where it showed higher levels of exposure. Contains NA for signatures not showing differential exposure.}
}

\examples{
# assuming signatures is the return value of signeR()

\dontrun{
    diff_exposure <- DiffExp(signatures$SignExposures, labels=my_labels, plotfile="Diffexp_boxplot.pdf")
}

# see also
vignette(package="signeR")
}

