\name{signeR}
\alias{signeR}
\title{signeR}
\description{
    Generates the signatures.
}
\usage{
    signeR(M, Mheader = TRUE, samples = "rows", Opport = NA,
        Oppheader = FALSE, nsig = NA, nlim = c(NA, NA),
        try_all = FALSE, ap = NA, bp = NA, ae = NA, be = NA,
        lp = NA, le = NA, var.ap = 10, var.ae = 10,
        testing_burn = 1000, testing_eval = 1000, EM_eval = 100,
        main_burn = 10000, main_eval = 2000, start = "lee",
        estimate_hyper = FALSE, EMit_lim=100)
}
\arguments{
    \item{M}{mutation counts matrix of samples x features.}
    \item{Mheader}{if M have colnames defined use TRUE, if FALSE a default
order will be assumed.}
    \item{samples}{if the samples are row-wise or column-wise in M, default
is "row".}
    \item{Opport}{context count matrix of samples x features in the target
genome or region.}
    \item{Oppheader}{if Opport have header defined.}
    \item{nsig}{number of signatures, which can be provided or estimated by
the algorithm.}
    \item{nlim}{define a interval to search for the optimal number of
signatures.}
    \item{try_all}{if true, all possible values for nsig will be tested}
    \item{ap}{shape parameter of the gamma distribution used to generate the
entries of a matrix of rate parameters of the gamma distributions which
generate signatures.}
    \item{bp}{rate parameter of the gamma distribution used to generate the
entries of a matrix of rate parameters of the gamma distributions which
generate signatures.}
    \item{ae}{shape parameter of the gamma distribution used to generate the
entries of a matrix of rate parameters of the gamma distributions which
generate exposures.}
    \item{be}{rate parameter of the gamma distribution used to generate the
entries of a matrix of rate parameters of the gamma distributions which
generate exposures.}
    \item{lp}{parameter of the exponential distribution used to generate the
entries of a matrix of shape parameters of the gamma distributions which
generate signatures.}
    \item{le}{parameter of the exponential distribution used to generate the
entries of a matrix of shape parameters of the gamma distributions which
generate exposures.}
    \item{var.ap}{variance of the gamma distribution used to generate
proposals for shape parameters of signatures}
    \item{var.ae}{variance of the gamma distribution used to generate
proposals for shape parameters of exposures}
    \item{testing_burn}{number of burning iterations of the Gibbs sampler used
to estimate the number of signatures in data. Corresponds to R0 at Algorithm 1
on signeR paper.}
    \item{testing_eval}{number of iterations of the Gibbs sampler used to
estimate the number of signatures in data. Corresponds to R2 at Algorithm 1
on signeR paper.}
    \item{EM_eval}{number of samples generated at each iteration of the EM
algorithm. Corresponds to R1 at Algorithm 1 on signeR paper.}
    \item{main_burn}{number of burning iterations of the final Gibbs sampler.}
    \item{main_eval}{number of iterations of the final Gibbs sampler.}
    \item{start}{NMF algorithm used to generate initial values for signatures
and exposures,options: "brunet","KL","lee","Frobenius","offset","nsNMF",
"ls-nmf","pe-nmf","siNMF","snmf/r" or "snmf/l".}
    \item{estimate_hyper}{if TRUE, algorithm estimates optimal values of
ap,bp,ae,be,lp,le. Start values can still be provided.}
    \item{EMit_lim}{limit of EM iterations for the estimation of
hyper-hyperparameters ap,bp,ae,be,lp,le. Default is 100. Corresponds to U at
Algorithm 1 on signeR paper.}
}
\value{
    signeR output is a list with the following items:\cr
    \item{Nsign}{selected number of signatures.}
    \item{tested_n}{array containing the numbers of signatures tested by the
algorithm.}
    \item{Test_BICs}{list of measured BIC values when testing different
numbers of signatures.}
    \item{Phat}{Estimated signatures, median of P samples.}
    \item{Ehat}{Estimated exposures, median of E samples.}
    \item{SignExposures}{SignExp object which contain the set of samples for
the model parameters.}
    \item{Bics}{measured BIC values on the final run of the sampler.}
    \item{HyperParam}{evolution of estimated hyperparameters when testing
different numbers of signatures.}
}
\examples{
vignette(package="signeR")
}
