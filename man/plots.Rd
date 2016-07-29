\name{plots}
\docType{methods}
\alias{SignPlot}
\alias{SignHeat}
\alias{Paths}
\alias{ExposureBoxplot}
\alias{BICboxplot}
\alias{SignPlot,SignExp-method}
\alias{SignHeat,SignExp-method}
\alias{Paths,SignExp-method}
\alias{ExposureBoxplot,SignExp-method}

\title{signeR plot functions}

\description{
    SignPlot: Plot the mutational signatures in a barchart, with error bars 
        according to the variation of individual entries along generated
        Gibbs samples.\cr

    SignHeat: Plot the mutations signatures in a heatmap.\cr

    ExposureBoxplot: Boxplot of exposure values, showing their variation
        along generated Gibbs samples.\cr

    Paths   : Plot the convergence of the Gibbs sampler for signatures and
        exposures on separate charts.\cr

    BICboxplot : Plot the measured values of the Bayesian Information
        Criterion (BICs) for tested model dimensions.
}

\usage{
\S4method{SignPlot}{SignExp}(signexp_obj, plot_to_file=FALSE,
    file="Signature_plot.pdf", pal="bcr1", threshold=0, plots_per_page=4,
    gap=1, reord=NA, ...)
\S4method{SignHeat}{SignExp}(signexp_obj, plot_to_file=FALSE,
    file="Signature_heatmap.pdf", nbins=20, pal="roh", ...)
\S4method{ExposureBoxplot}{SignExp}(signexp_obj, plot_to_file=FALSE, 
    file="Exposure_boxplot.pdf", col='tan2', threshold=0,
    plots_per_page=4, ...)
\S4method{Paths}{SignExp}(signexp_obj, plot_to_file=FALSE,
    file_suffix="plot.pdf", plots_per_page=4, ...)
BICboxplot(signeRout, plot_to_file=FALSE, file="Model_selection_BICs.pdf")
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.
e.g.: sig$SignExposures}
    \item{signeRout}{the value returned by the signeR function.}
    \item{plot_to_file}{Whether to save the plot to the file parameter.
Default is FALSE.}
    \item{file}{Output pdf file of the plots.}
    \item{pal}{Color palette used. Options are: "brew","lba","bcr1",
"bcr2","bw","roh".}
    \item{threshold}{entries below this value will be rounded to 0. Default
is 0 (all entries are kept).}
    \item{plots_per_page}{How many plots in a single page, default is 4.}
    \item{gap}{Distance between consecutive bars on the plot.}
    \item{reord}{Order of signatures for plotting. Should be a permutation of
1:nsig, where nsig is the number of signatures. By default, signatures are
ordered by the total exposure, in decreasing order.}
    \item{nbins}{The range of signature entries is divided in this number of
bins for plotting, each bin corresponding to a different color.}
    \item{file_suffix}{The suffix of the output file.}
    \item{col}{Single color name for boxplots.}
    \item{...}{.}
}

\value{
The plot result is saved in the file defined by the file argument.
}

\examples{
\dontshow{
    # mock data
    mut <- rmultinom(96, 100, 1:4)
    signatures <- signeR(M=mut, Mheader=FALSE, EMit_lim=10, nsig=1)
}

# each plot function need the SignExposures object
# which is part of the result of the signeR() call
SignPlot(signatures$SignExposures)
Paths(signatures$SignExposures)
# etc ..

# BICboxplot needs the returned value itself
BICboxplot(signatures)

# see also
vignette(package="signeR")
}

