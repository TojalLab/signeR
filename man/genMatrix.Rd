\name{generateMatrix}
\alias{genCountMatrixFromVcf}
\alias{genOpportunityFromGenome}

\title{count matrix and opportinity matrix generators}

\description{
    genCountMatrixFromVcf    : generate count matrix from a VCF file.\cr
    genOpportunityFromGenome : generate opportunity matrix from a target regions set.
}

\usage{
    genCountMatrixFromVcf(bsgenome, vcfobj)
    genOpportunityFromGenome(bsgenome, target_regions, nsamples=1)
}

\arguments{
    \item{bsgenome}{A BSgenome object, equivalent to the genome used for the variant call.}
    \item{vcfobj}{A VCF object. See VCF-class from the VariantAnnotation package.}
    \item{target_regions}{A GRanges object, describing the target region analyzed by the variant caller.}
    \item{nsamples}{Number of samples to generate the matrix, should be the same number as rows of the count matrix.}
}

\value{
A matrix of samples x (96 features).\cr
Each feature is a SNV change with a 3bp context.
}

\examples{
\dontrun{
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(rtracklayer)
    library(VariantAnnotation)

    vcfobj <- readVcf("/path/to/a/file.vcf", "hg19")
    mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)
    
    target_regions <- import(con="/path/to/a/target.bed", format="bed")
    opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19, target_regions, nsamples=nrows(mut))
}

# see also
vignette(package="signeR")
}
