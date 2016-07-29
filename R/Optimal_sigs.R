Optimal_sigs<-function(testfun,liminf,limsup,step,oldpoints=c(),oldvalues=c(),
    oldextra=list(),left_eval=NA){
    points<-seq(liminf,limsup,by=step)
    if(step>1){
        if (points[length(points)]<limsup) points<-c(points,limsup)
        maxpos<-1
        values<-rep(0,length(points))
        extra_out<-list()
        downhill<-FALSE
        k=0
        while(!downhill & k<length(points)){
            k=k+1
            if (points[k] %in% oldpoints){
                thisval<-oldvalues[oldpoints==points[k]]
                thisextra<-oldextra[[which(oldpoints==points[k])]]
            }else{
                if (is.na(left_eval)){
                    cat(paste("Evaluating ",points[k]," signatures. \n",sep=""))
                }else{
                    if ( left_eval > 1 ){
                        cat(paste("Evaluating ",points[k]," signatures (",
                            left_eval," evaluations left). \n",sep=""))
                    }else{
                        cat(paste("Evaluating ",points[k],
                            " signatures (last evaluation). \n",sep=""))
                    }
                    left_eval <- left_eval - 1
                }
                eval<-testfun(points[k])#ebayesNMF(...,n=points[k])
                thisval<-eval[[1]]
                thisextra<-eval[-1]
            }
            values[k]<-thisval
            extra_out[[k]]<-thisextra
            if (thisval>=values[maxpos]){
                maxpos<-k
            }else if ( k>(maxpos+1) ){
                downhill<-TRUE
            }
        }
        points<-points[1:k]
        values<-values[1:k]
        extra_out<-extra_out[1:k]
        maxvals<-which(values==max(values))
        newinf<-points[max(min(maxvals)-1,1)]
        newsup<-points[min(max(maxvals)+1,length(points))]
        newstep<-floor(step/2)
        cat(paste("Refining search for the number of signatures ranging from ",
            newinf," to ",newsup,", please be patient.\n",sep=""))
        Os<-Optimal_sigs(testfun,liminf=newinf,limsup=newsup,step=newstep,
            oldpoints=points,oldvalues=values,oldextra=extra_out,
            left_eval=2*floor(log2(newstep)+1))
        n<-Os[[1]]
        old_small<- points < min(Os[[2]])
        old_big  <- points > max(Os[[2]])
        points<-c(points[old_small],Os[[2]],points[old_big])
        values<-c(values[old_small],Os[[3]],values[old_big])
        extra_out<-c(extra_out[old_small],Os[[4]],extra_out[old_big])
    }else{
        values<-rep(0,length(points))
        extra_out<-list()
        for(k in 1:length(points)){
            if (points[k] %in% oldpoints){
                thisval<-oldvalues[oldpoints==points[k]]
                thisextra<-oldextra[[which(oldpoints==points[k])]]
            }else{
                if (is.na(left_eval)){
                    cat(paste("Evaluating ",points[k]," signatures. \n",sep=""))
                }else{
                    if ( left_eval > 1 ){
                        cat(paste("Evaluating ",points[k]," signatures (",
                            left_eval," evaluations left). \n",sep=""))
                    }else{
                        cat(paste("Evaluating ",points[k],
                            " signatures (last evaluation). \n",sep=""))
                    }
                    left_eval <- left_eval - 1
                }
                eval<-testfun(points[k])#ebayesNMF(...,n=points[k])
                thisval<-eval[[1]]
                thisextra<-eval[-1]
            }
            values[k]<-thisval
            extra_out[[k]]<-thisextra
        }
        n<-points[which(values==max(values))[1]]
    }
    return(list(n,points,values,extra_out))
}
