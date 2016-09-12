\name{methods}
\docType{methods}
\alias{setSamples}
\alias{setMutations}
\alias{Normalize}
\alias{Reorder}
\alias{Average_sign}
\alias{Median_sign}
\alias{Average_exp}
\alias{Median_exp}
\alias{setSamples,SignExp-method}
\alias{setMutations,SignExp-method}
\alias{Normalize,SignExp-method}
\alias{Reorder,SignExp,numeric-method}
\alias{Average_sign,SignExp-method}
\alias{Median_sign,SignExp-method}
\alias{Average_exp,SignExp-method}
\alias{Median_exp,SignExp-method}

\title{SignExp class methods}

\description{
    setSamples: Define sample names for a SignExp object, according to the 
        "names" argument.\cr
        
    setMutations: Define mutation names for a SignExp object, according to the 
        "mutations" argument.\cr
        
    Normalize: Normalize a SignExp object so that the entries of each 
        signature sum up to one.\cr
        
    Reorder: Change the order of the signatures in a SignExp object. New 
    signature order wil be defined by the "ord" argument.\cr
    
    Average_sign: Exports an aproximation of the signatures obtained by the 
        averages of the samples for the signature matrix P.\cr
        
    Median_sign: Exports an aproximation of the signatures obtained by the 
        medians of the samples for signature matrix P.\cr
        
    Average_exp: Exports an aproximation of the exposures obtained by the 
        averages of the samples for exposure matrix E.\cr
        
    Median_exp: Exports an aproximation of the exposures obtained by the 
        medians of the samples for exposure matrix E.
}

\usage{
\S4method{setSamples}{SignExp}(signexp_obj, names)
\S4method{setMutations}{SignExp}(signexp_obj, mutations)
\S4method{Normalize}{SignExp}(signexp_obj)
\S4method{Reorder}{SignExp,numeric}(signexp_obj, ord)
\S4method{Average_sign}{SignExp}(signexp_obj, normalize=TRUE)
\S4method{Median_sign}{SignExp}(signexp_obj, normalize=TRUE)
\S4method{Average_exp}{SignExp}(signexp_obj, normalize=TRUE)
\S4method{Median_exp}{SignExp}(signexp_obj, normalize=TRUE)
}

\arguments{
    \item{signexp_obj}{a SignExp object returned by signeR function.
e.g.: sig$SignExposures}
    \item{names}{Vector of sample names.}
    \item{mutations}{Vector of mutations, e.g. "C>A:TCG".}
    \item{normalize}{Whether the signatures should be normalized before 
    extracting approximations. Default is TRUE.}
    \item{ord}{Vector with the new signature order.}
}

\value{ setSamples, setMutations, Normalize and Reorder return a modified 
SignExp object.\cr 
Average_sign, Median_sign, Average_exp and Median_exp return a
matrix with the corresponding approximation.
}

\examples{
\dontshow{
    # mock data
    mut <- rmultinom(96, 50, 1:4)
    signatures <- signeR(M=mut, Mheader=FALSE, nsig=2, EMit_lim=10, EM_eval=10, main_eval=10, main_burn=10)
}

# each function needs the SignExposures object
# which is part of the result of the signeR() call
signexp  <- Normalize(signatures$SignExposures)
signexp  <- Reorder(signatures$SignExposures,ord=c(2,1))
matrix_p <- Median_sign(signatures$SignExposures)
# etc ...

# see also
vignette(package="signeR")
}

