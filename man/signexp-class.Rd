\name{SignExp}
\docType{class}
\alias{SignExp}
\alias{SignExp-class}

\title{SignExp class}

\description{
Keep samples for signature and exposure matrices.
}

\value{
Object fields:
\item{@Sign}{array of signature matrix samples.}
\item{@Exp}{array of exposure matrix samples.}
\item{@sigSums}{Signature sums for each sample, organized by row.
Normalizing factors.}
\item{@samples}{Genome sample IDs.}
\item{@mutations}{mutation names.}
\item{@normalized}{boolean variable, indicating whether Sign array has been
normalized.}
}
