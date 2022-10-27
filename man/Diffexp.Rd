\name{DiffExp}
\docType{methods}
\alias{DiffExp}
\alias{DiffExp,SignExp,character-method}

\title{Differential Exposure Analysis}

\description{
    DiffExp : Identify signatures with significantly different activities
among sample groups.
}

\usage{
\S4method{DiffExp}{SignExp,character}(signexp_obj, labels, max_instances=200, 
    method=kruskal.test, contrast="all", quant=0.5, cutoff=0.05,
    p.adj= "BH",plot_to_file=FALSE, file="Diffexp_boxplot.pdf",
    colored=TRUE, relative = FALSE, ...)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.}
    \item{labels}{sample labels used to define sample groups.}
    \item{max_instances}{Maximum number of the exposure matrix instances to be 
    analyzed. If the number of available E instances is bigger then this 
    parameter, a subset of those will be randomly selected for analysis.}
    \item{method}{algorithm used to compare each signature exposures among
sample groups. Default is kruskal.test, which leads to the use of
Kruskal-Wallis Rank Sum Test.}
    \item{contrast}{defines which sample groups will be considered in the
analysis. Default is "all", which leads the algorithm to evaluates the null
hypothesis of exposure levels beeing constant in all groups. Instead, if this
parameter contains a list of group labels, the algorithm will evaluate the
null hypothesis of exposure levels beeing constant among those groups.}
    \item{quant}{the p-values quantile which, after log-transform, will be
used as DES (Differental Exposure Score). Deafult is 0.5, which means the
median log-transformed p-value will be considered as DES.}
    \item{p.adj}{correction method for p-values adjust at the post-hoc tests
    performed when there are more then two group labels. See p.adjust for
    options.}
    \item{cutoff}{threshold for p-values quantile for signatures to be
considered as showing differential exposure.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
Default is FALSE.}
    \item{file}{Output file to export p-values boxplot.}
    \item{colored}{Boolean variable, if TRUE boxplots of differentially
exposed signatures will be colored in green, cutoff line will be colored in
red and line segments showing the transformed p-value quantile used for DE
evaluation will be colored in blue. Otherwise the plot will be black & white.}
    \item{relative}{Whether tests should be performed on absolute or relative
    signature contributions to each sample mutations.
    Default is FALSE (absolute contributions will be tested).}
    \item{...}{additional parameters for test algorithm defined by the method
parameter.}
}

\value{
A list with the following items:
\item{Pvquant}{boolean array with one entry for each signature, indicating
whether it shows differential exposure.}
\item{Pvalues}{matrix containing all computed p-values, with one row for each
signature.}
\item{MostExposed}{for each differentially exposed signature, this array
contains the label of the group where it showed higher levels of exposure.
Contains NA for signatures not showing differential exposure.}
\item{Differences}{List of matrices, exported only when there are more than two
groups in the analysis and any signature is found to be differentially active.
Each matrix corresponds to one of the highlighted signatures and show the
results of comparisons among groups, with the significant ones marked as TRUE.}
}

\examples{
# assuming signatures is the return value of signeR()
\dontshow{
    # mock data
    mut <- t(cbind(rmultinom(2, 500, rep(1:6,each=16)),
                   rmultinom(2, 500, rep(6:1,each=16))))
    signatures <- signeR(M=mut, Mheader=FALSE, nsig=2, EMit_lim=10, 
                         EM_eval=10, main_eval=10, main_burn=10)
}
# labels vector, one for each sample
my_labels <- c("a","a","b","b")

diff_exposure <- DiffExp(signatures$SignExposures,labels=my_labels)

# see also
vignette(package="signeR")
}

