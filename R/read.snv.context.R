read.snv.context<-function(stvet){
    triplets<-c()
    mutations<-c()
    for (st in stvet){
        if (grepl("[ACGT]{3,3}",st,ignore.case=TRUE,perl=TRUE)){
            strt<-regexpr("[ACGT]{3,3}",st)[[1]]
            trivet<-strsplit(st,"")[[1]][strt:(strt+2)]
            triplet<-paste(trivet,collapse="")
            st<-sub(triplet,"",st,ignore.case=TRUE)
            st<-sub(trivet[2],"",st,ignore.case=TRUE)
            strt<-regexpr("[ACGT]",st,ignore.case=TRUE)[[1]]
            newbase<-strsplit(st,"")[[1]][strt]
            mut<-paste(trivet[2],newbase,sep=">")
        }else if(grepl("[ACTG]>[ACTG]",st,ignore.case=TRUE,perl=TRUE)){
            strt<-regexpr("[ACTG]>[ACTG]",st,ignore.case=TRUE)[[1]]
            mutvet<-strsplit(st,"")[[1]][strt:(strt+2)]
            mut<-paste(mutvet,collapse="")
            st<-sub(mut,"",st,ignore.case=TRUE)
            splitvet<-strsplit(st,"")[[1]]
            splitvet<-splitvet[splitvet %in% c("A","C","G","T","a","c","g","t")]
            firstbase<-splitvet[1]
            lastbase<-splitvet[2]
            triplet<-paste(c(firstbase,mutvet[1],lastbase),collapse="")
        }else if (grepl("[ACTG]\\.[ACTG]",st,ignore.case=TRUE,perl=TRUE)){
            strt<-regexpr("[ACTG]\\.[ACTG]",st,ignore.case=TRUE)[[1]]
            trivet<-strsplit(st,"")[[1]][strt:(strt+2)]
            st<-sub(paste(trivet,collapse=""),"",st,ignore.case=TRUE)
            splitvet<-strsplit(st,"")[[1]]
            splitvet<-splitvet[splitvet %in% c("A","C","G","T","a","c","g","t")]
            oldbase<-splitvet[1]
            newbase<-splitvet[2]
            trivet[2]<-oldbase
            triplet<-paste(trivet,collapse="")
            mut<-paste(oldbase,newbase,sep=">")
        }else{
            triplet<-st
            mut<-""
        }
        triplets<-c(triplets,triplet)
        mutations<-c(mutations,mut)
    }
    return(list(triplets,mutations))
}
