\name{signeR-package}
\alias{signeR-package}
\docType{package}
\title{
\packageTitle{signeR}
}
\description{
\packageDescription{signeR}
}
\details{
signeR package focus on the characterization and analysis of mutational processes. Its functionalities can be divided in three steps. Firstly, it provides tools to process VCF files and generate matrices of SNV mutation counts and mutational opportunities, both divided according to a 3bp context (mutation site and its neighboring bases). Secondly, the main part of the package takes those matrices as input and applies a Bayesian approach to estimate the number of underlying signatures and their mutational profiles. Thirdly, the package provides tools to correlate the activities of those signatures with other relevant information, e.g. clinical data, in order to infer conclusions about the analyzed genome samples, which can be useful for clinical applications.
}
\author{
Rodrigo Drummond, Rafael Rosales, Renan Valieris, Israel Tojal da Silva

Maintainer: \packageMaintainer{signeR}
}
\references{
here be references
}
\examples{
vignette(package="signeR")
}
\keyword{ package }
