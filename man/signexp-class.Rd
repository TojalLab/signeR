\name{SignExp}
\docType{class}
\alias{SignExp}

\title{SignExp class}

\description{
Keep samples for signature and exposure matrices.
}

\value{
Object fields:
\item{.Sign}{tensor of signature matrix samples.}
\item{.Exp}{tensor of exposure matrix samples.}
\item{.sigSums}{Signature sums for each sample, organized by row. Normalizing factors.}
\item{.samples}{Genome sample IDs.}
\item{.mutations}{mutation names.}
\item{.normalized}{boolean variable, indicating whether Sign tensor has been normalized.}
}

