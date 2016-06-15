\name{plots}
\alias{SignPlot.SignExp}
\alias{SignHeat.SignExp}
\alias{Paths.SignExp}
\alias{ExposureBoxplot.SignExp}
\alias{SignPlot}
\alias{SignHeat}
\alias{Paths}
\alias{ExposureBoxplot}
\alias{BICboxplot}

\title{signeR plot functions}

\description{
    SignPlot: Plot the mutational signatures in a barchart, with error bars 
        according to the variation of individual entries along generated Gibbs samples.\cr

    SignHeat: Plot the mutations signatures in a heatmap.\cr

    ExposureBoxplot: Boxplot of exposure values, showing their variation along generated Gibbs samples.\cr

    Paths   : Plot the convergence of the Gibbs sampler for signatures and exposures on separate charts.\cr

    BICboxplot : Plot the measured values of the Bayesian Information Criterion (BICs) for tested model dimensions.
}

\usage{
\S3method{SignPlot}{SignExp}(this, plotfile="Signature_plot.pdf", pal="bcr1",
    threshold=0, plots_per_page=4, gap=1, reord=NA, ...)
\S3method{SignHeat}{SignExp}(this, plotfile="Signature_heatmap.pdf", nbins=20, ...)
\S3method{ExposureBoxplot}{SignExp}(this, plotfile="Exposure_boxplot.pdf", col='tan2',
    threshold=0, plots_per_page=4, ...)
\S3method{Paths}{SignExp}(this, file_suffix="plot.pdf", plots_per_page=4, ...)
BICboxplot(signeRout, plotfile="Model_selection_BICs.pdf")
}

\arguments{
    \item{this}{a SignExp object returned by signeR function. e.g.: sig$SignExposures}
    \item{signeRout}{the value returned by the signeR function.}
    \item{plotfile}{Output pdf file of the plots.}
    \item{pal}{Color palette used. Options are: "brew","lba","bcr1","bcr2","bw".}
    \item{threshold}{entries below this value will be rounded to 0. Default is 0 (all entries are kept).}
    \item{plots_per_page}{How many plots in a single page, default is 4.}
    \item{gap}{Distance between consecutive bars on the plot.}
    \item{reord}{Order of signatures for plotting. Should be a permutation of 1:nsig, where nsig is the number of signatures. By default, signatures are ordered by the total exposure, in decreasing order.}
    \item{nbins}{The range of signature entries is divided in this number of bins for plotting, each bin corresponding to a different color.}
    \item{file_suffix}{The suffix of the output file.}
    \item{col}{Single color name for boxplots.}
    \item{...}{.}
}

\value{
The plot result is saved in the file defined by the file argument.
}

\examples{
# assuming signatures is the return value of signeR()

\dontrun{
    SignPlot(signatures$SignExposures)
    Paths(signatures$SignExposures)
}

vignette(package="signeR")
}

