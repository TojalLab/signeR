Comp_labels<-function(Equival){ 
    # Find group labels indicating classes equivalence. 
    # Equival is an array with results for all significant signatures, 
    # where TRUE indicates classes significantly different.
    nclasses<-dim(Equival)[2]+1
    Allcomp<-list()
    for (i in 1:dim(Equival)[3]){
        Thiseq <- Equival[,,i]
        multicomp<-vector("list", nclasses)
        labelmat<-matrix(1,nclasses,1)
        for(cl in 1:(nclasses-1)){
            vtest<-as.vector(Thiseq[cl:(nclasses-1),cl])
            if(sum(vtest)>0){
                differ<-cl+which(vtest)
                for (d in differ){
                    update<-apply(labelmat,2,function(v){v[cl]==1 & v[d]==1})
                    newcols<-labelmat[,update,drop=FALSE]
                    labelmat[d,update]<-0
                    newcols[cl,]<-0
                    labelmat<-cbind(labelmat,newcols)
                    keep<-apply(labelmat,2,function(lab){
                        sum(t(labelmat) %*% lab == sum(lab))==1
                    }) 
                    labelmat<-labelmat[,keep]    
                }
            }
        }
        for(cl in 1:nclasses){
            multicomp[[cl]]<-letters[which(labelmat[cl,]==1)]
        }
        Allcomp[[i]]<-multicomp
    }
    return(Allcomp)
}
