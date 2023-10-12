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
signeR package focuses on the characterization and analysis of mutational
processes. Its functionalities can be divided into three steps. Firstly, it
provides tools to process VCF files and generate matrices of SNV mutation
counts and mutational opportunities, both divided according to a 3bp context
(mutation site and its neighboring bases). Secondly, the main part of the
package takes those matrices as input and applies a Bayesian approach to
estimate the number of underlying signatures and their mutational profiles.
Thirdly, the package provides tools to correlate the activities of those
signatures with other relevant information, e.g. clinical data, to
infer conclusions about the analyzed genome samples, which can be useful for
clinical applications.
}
\author{
Rodrigo Drummond, Rafael Rosales, Renan Valieris, Israel Tojal da Silva

Maintainer: \packageMaintainer{signeR}
}
\references{
    This work has been submitted to Bioinformatics under the title
    "signeR: An empirical Bayesian approach to mutational signature discovery".

    L. B. Alexandrov, S. Nik-Zainal, D. C. Wedge, P. J. Campbell, and M.
    R. Stratton. Deciphering Signatures of Mutational Processes Operative in
    Human Cancer. Cell Reports, 3(1):246-259, Jan. 2013.
    doi:10.1016/j.celrep.2012.12.008.

    A. Fischer, C. J. Illingworth, P. J. Campbell, and V. Mustonen.
    EMu: probabilistic inference of mutational processes and their
    localization in the cancer genome. Genome biology, 14(4):R39, Apr. 2013.
    doi:10.1186/gb-2013-14-4-r39.
}
\examples{
vignette(package="signeR")
}
\keyword{ package }
