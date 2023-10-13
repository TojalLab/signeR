
change_triplet <- c(
    "C>A:ACA","C>A:ACC","C>A:ACG","C>A:ACT","C>A:CCA","C>A:CCC","C>A:CCG",
    "C>A:CCT","C>A:GCA","C>A:GCC","C>A:GCG","C>A:GCT","C>A:TCA","C>A:TCC",
    "C>A:TCG","C>A:TCT","C>G:ACA","C>G:ACC","C>G:ACG","C>G:ACT","C>G:CCA",
    "C>G:CCC","C>G:CCG","C>G:CCT","C>G:GCA","C>G:GCC","C>G:GCG","C>G:GCT",
    "C>G:TCA","C>G:TCC","C>G:TCG","C>G:TCT","C>T:ACA","C>T:ACC","C>T:ACG",
    "C>T:ACT","C>T:CCA","C>T:CCC","C>T:CCG","C>T:CCT","C>T:GCA","C>T:GCC",
    "C>T:GCG","C>T:GCT","C>T:TCA","C>T:TCC","C>T:TCG","C>T:TCT","T>A:ATA",
    "T>A:ATC","T>A:ATG","T>A:ATT","T>A:CTA","T>A:CTC","T>A:CTG","T>A:CTT",
    "T>A:GTA","T>A:GTC","T>A:GTG","T>A:GTT","T>A:TTA","T>A:TTC","T>A:TTG",
    "T>A:TTT","T>C:ATA","T>C:ATC","T>C:ATG","T>C:ATT","T>C:CTA","T>C:CTC",
    "T>C:CTG","T>C:CTT","T>C:GTA","T>C:GTC","T>C:GTG","T>C:GTT","T>C:TTA",
    "T>C:TTC","T>C:TTG","T>C:TTT","T>G:ATA","T>G:ATC","T>G:ATG","T>G:ATT",
    "T>G:CTA","T>G:CTC","T>G:CTG","T>G:CTT","T>G:GTA","T>G:GTC","T>G:GTG",
    "T>G:GTT","T>G:TTA","T>G:TTC","T>G:TTG","T>G:TTT"
)

maf_req_columns <- c("Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2",
  "Variant_Type","Chromosome","Start_Position","End_Position","Strand",
  "Tumor_Sample_Barcode")

revcomp <- function(s) {
    return(as.character(reverseComplement(DNAString(s))))
}

isBadGeno <- function(genoIndex) {
    if(isGenoIndexRef(genoIndex)){return(FALSE)}
    if("." %in% genoIndex){return(TRUE)}
    if(length(genoIndex) > 2){return(TRUE)}
    return(FALSE)
}

isGenoIndexRef <- function(genoIndex) {
    if(length(genoIndex) == 0){return(TRUE)}
    return(length(genoIndex) == 1
        && (genoIndex[[1]] == "." || genoIndex[[1]] == "0"))
}

getFirstGenoAltIndex <- function(genoIndex) {
    return(Filter(function(x)x!="0",genoIndex)[1])
}

genOpportunityFromGenome <- function(bsgenome, target_regions, nsamples=1) {

    # make sure there are no overlaps
    target_regions <- reduce(target_regions, drop.empty.ranges=TRUE)

    # count the kmers in the region
    kmers <- colSums(oligonucleotideFrequency(
        getSeq(bsgenome, target_regions),3))

    # turn back into matrix
    c0 <- names(kmers)
    dim(kmers) <- c(1, length(kmers))
    colnames(kmers) <- c0

    uk <- sapply(colnames(kmers),
        function(x){y<-substr(x,2,2);if(y=="C"||y=="T"){x}else{revcomp(x)}})

    m <- matrix(0, nrow=1, ncol=length(change_triplet))
    colnames(m) <- change_triplet
    for(i in 1:length(kmers)) {
        v <- kmers[i]
        k <- uk[[i]]
        for(j in grep(k,change_triplet)) {
            m[1,j] = v + m[1,j]
        }
    }

    # replicate the matrix over the number of desired samples
    m <- do.call("rbind",rep(list(m),nsamples))
    return(m)
}

genCountMatrixFromVcf <- function(bsgenome, vcfobj) {

    # keep only SNVs
    vcfobj <- vcfobj[isSNV(vcfobj),]

    # convert alt() to DNAString if needed
    # https://github.com/TojalLab/signeR/issues/17
    alt(vcfobj) <- DNAStringSetList(alt(vcfobj))

    # remove mutations on chrom edges
    # but not if the chrom is circular
    si = BiocGenerics::as.data.frame(seqinfo(bsgenome))
    si = si[(is.na(si$isCircular) | si$isCircular==FALSE),]
    si = GRanges(seqinfo(bsgenome)[rownames(si)])
    if(length(subsetByOverlaps(vcfobj, si, type='start')) > 0 ||
        length(subsetByOverlaps(vcfobj, si, type='end')) > 0) {
        warning("Mutations on the edge of chromosomes will be ignored.")
        vcfobj = subsetByOverlaps(vcfobj, si, type='start', invert=T)
        vcfobj = subsetByOverlaps(vcfobj, si, type='end', invert=T)
    }

    contexts <- getSeq(bsgenome, resize(granges(vcfobj), 3, fix="center"))
    alts <- alt(vcfobj)
    refs <- ref(vcfobj)

    gtsMat <- geno(vcfobj)$GT
    gtsMat <- structure(lapply(strsplit(gtsMat,"[/|]",fixed=FALSE), unique), dim=dim(gtsMat))
    sample_names <- colnames(vcfobj)

    count_matrix <- matrix(0,
        nrow=length(sample_names), ncol=length(change_triplet))

    rownames(count_matrix) <- sample_names
    colnames(count_matrix) <- change_triplet

    for(i in 1:nrow(gtsMat)) {
        rb <- refs[[i]]
        cc <- contexts[i]
        if(rb == DNAString("C") || rb == DNAString("T")) {
            cts <- sapply(alts[i], function(ab){paste0(rb,">",ab,":",cc)})
        } else {
            cts <- sapply(alts[i], function(ab){
                paste0(reverseComplement(rb),">",
                    reverseComplement(ab),":",
                    reverseComplement(cc))})
        }
        for(j in 1:length(sample_names)) {
            gi <- gtsMat[[i,j]]
            if(isGenoIndexRef(gi)) {next}
            if(isBadGeno(gi)){
                warning(paste0("Genotype of Line ",i," sample ",
                    sample_names[j]," is not supported, it will be skipped."))
                next
            }
            ai <- as.integer(getFirstGenoAltIndex(gi))
            ct <- cts[[ai]]
            cx <- match(ct, change_triplet)
            if(is.na(cx)) {
                warning(paste0("ct = ",ct))
                next
            }
            count_matrix[j,cx] = 1 + count_matrix[j,cx]
        }
    }
    return(count_matrix)
}

genCountMatrixFromMAF <- function(bsgenome, maf_file) {
    
    maf <- readr::read_tsv(maf_file, col_types='')

    # assert all required columns are present
    if(!all(maf_req_columns %in% colnames(maf))) {
        dc <- paste(collapse = ',',setdiff(maf_req_columns, colnames(maf)))
        stop(paste("Error: Missing columns in input file: ",dc))
    }

    # keep only SNVs
    maf <- dplyr::filter(maf, Variant_Type == 'SNP')

    maf <- dplyr::mutate(maf,
        ref=Reference_Allele,
        alt=ifelse(Tumor_Seq_Allele2=="", Tumor_Seq_Allele1, Tumor_Seq_Allele2))
    maf <- dplyr::filter(maf, ref!=alt)

    # assert all chrom from maf are present on the genome
    # attempt to add or remove "chr" if the names dont match
    si <- seqinfo(bsgenome)
    if(!all(unique(maf$Chromosome) %in% si@seqnames)) {
        if(all(paste0("chr",unique(maf$Chromosome)) %in% si@seqnames)) {
            warning("Warning: sequence names from MAF don't match genome sequence names, trying to continue by adding the chr prefix.")
            maf <- dplyr::mutate(maf, Chromosome=paste0("chr",Chromosome))
        } else if(all(gsub("^chr","",unique(maf$Chromosome)) %in% si@seqnames)) {
            warning("Warning: sequence names from MAF don't match genome sequence names, trying to continue by removing the chr prefix.")
            maf <- dplyr::mutate(maf, Chromosome=gsub("^chr","",Chromosome))
        } else {
            stop("Error: sequence names from MAF don't match genome sequence names, please check the genome parameter.")
        }
    }

    mut_ranges <- GRanges(maf$Chromosome, IRanges(maf$Start_Position, maf$End_Position), maf$Strand)
    contexts <- getSeq(bsgenome, resize(mut_ranges, 3, fix='center'))
    sample_names <- sort(unique(maf$Tumor_Sample_Barcode))

    refs <- DNAStringSet(maf$ref)
    alts <- DNAStringSet(maf$alt)
    count_matrix <- matrix(0, ncol=length(change_triplet), nrow=length(sample_names))
    rownames(count_matrix) = sample_names
    colnames(count_matrix) = change_triplet
    for(i in 1:nrow(maf)) {
        rb <- DNAString(maf[[i, "Reference_Allele"]])
        cc <- contexts[i]
   
        if(rb == DNAString("C") || rb == DNAString("T")) {
            ct <- sapply(alts[i], function(ab){paste0(rb,">",ab,":",cc)})
   
        } else {
            ct <- sapply(alts[i], function(ab) {
                paste0(reverseComplement(rb),">",
                       reverseComplement(ab),":",
                       reverseComplement(cc))})
        }
   
        j <- match(maf[[i, "Tumor_Sample_Barcode"]], sample_names)
        cx <- match(ct, change_triplet)
        count_matrix[j, cx] = 1 + count_matrix[j,cx]
    }

    return(count_matrix)
}
