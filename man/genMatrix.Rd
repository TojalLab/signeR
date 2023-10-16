\name{generateMatrix}
\alias{genCountMatrixFromVcf}
\alias{genCountMatrixFromMAF}
\alias{genOpportunityFromGenome}

\title{count matrix and opportunity matrix generators}

\description{
    genCountMatrixFromVcf    : generate a count matrix from a VCF file.\cr
    genCountMatrixFromMAF    : generate a count matrix from an MAF file.\cr
    genOpportunityFromGenome : generate an opportunity matrix from a target
regions set.
}

\usage{
    genCountMatrixFromVcf(bsgenome, vcfobj)
    genCountMatrixFromMAF(bsgenome, maf_file)
    genOpportunityFromGenome(bsgenome, target_regions, nsamples=1)
}

\arguments{
    \item{bsgenome}{A BSgenome object, equivalent to the genome used for
the variant call.}
    \item{vcfobj}{A VCF object. See VCF-class from the VariantAnnotation
package.}
    \item{maf_file}{Path to a MAF file.}
    \item{target_regions}{A GRanges object, describing the target region
analyzed by the variant caller.}
    \item{nsamples}{Number of samples to generate the matrix, should be the
same number as rows of the count matrix.}
}

\value{
A matrix of samples x (96 features).\cr
Each feature is an SNV change with a 3bp context.
}

\examples{
library(rtracklayer)
library(VariantAnnotation)

# input files, variant call and target
vcf_file <- system.file("extdata","example.vcf", package="signeR")
bed_file <- system.file("extdata","example.bed", package="signeR")

# BSgenome, will depend on your variant call
library(BSgenome.Hsapiens.UCSC.hg19)

vcfobj <- readVcf(vcf_file, "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con=bed_file, format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
    target_regions, nsamples=nrow(mut))


# see also
vignette(package="signeR")
}
