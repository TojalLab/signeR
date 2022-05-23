Optimal_sigs<-function(testfun,liminf,limsup,step,
  significance=FALSE,pcrit=0.05,parplan="multisession"){
  plan(parplan)
  points<-seq(liminf,limsup,by=1)
  values<-listenv()
  for(k in 1:length(points)) values[[k]]<-NA
  names(values)<-as.character(points)
  to_eval<-seq(liminf,limsup,by=step)
  if (to_eval[length(to_eval)]<limsup) to_eval<-c(to_eval,limsup)
  evalued<-c()
  ncores<-detectCores()-1
  downhill<-FALSE
  nmax<-c()
  if(length(to_eval)>ncores){
    eval_now<-to_eval[1:ncores]
    eval_latter<-to_eval[-c(1:ncores)]
  }else{
    eval_now<-to_eval
    eval_latter<-c()
  } 
  for(k in eval_now) values[[as.character(k)]] %<-% {testfun(k)}
  running<-eval_now
  while(!downhill & length(running)>0){ #large steps, looking for downhill descent
    still_running<-running
    for(k in still_running){
      if(resolved(futureOf(values[[as.character(k)]]))){
        thisbics<-(values[[as.character(k)]][[1]])
        thismedian<-median(thisbics)
        #compare k and nmax
        if(length(nmax)>0){
          if(significance){
            testeq<-sapply(nmax,function(n){
              maxbics<-values[[as.character(n)]][[1]]
              maxmed<-median(maxbics)
              kt<-kruskal.test(list(thisbics,maxbics))
              if(kt$p.value<=pcrit){
                if(thismedian>maxmed){ 
                  testval=1 
                }else{ 
                  testval=-1 
                }
              }else testval = 0
              return(testval)
            })
          }else{
            testeq<-sapply(nmax,function(n){
              maxbics<-values[[as.character(n)]][[1]]
              maxmed<-median(maxbics)
              if(thismedian>maxmed){ 
                testval=1 
              }else if(thismedian<maxmed){ 
                testval=-1 
              }else testval = 0
              return(testval)
            })
          }
          nmax<-nmax[testeq %in% c(0,-1)] #keep those n not worse than k
          if( length(nmax)==0 | min(testeq)>=0 ) nmax<-c(nmax,k) #keep k if no other n is better
        }else{
          nmax<-k
        }### end of comparison k and nmax
        running<-running[!running==k]
        evalued<-c(evalued,k)
        if(length(eval_latter)>0){ 
          #start testing another number of signatures:
          nextrun<-eval_latter[1]
          values[[as.character(nextrun)]] %<-% {testfun(nextrun)}
          running<-c(running,nextrun)
          eval_latter<-eval_latter[-1]
        }
      }
    }
    #check downhill
    if(length(nmax)>0 & length(evalued)>=3){ 
      if(sum(evalued>max(nmax))>=2) {downhill<-TRUE}
    }
  }
  while(step>1){ #refining steps
    step<-ceiling(step/2)        
    left_eval=2*floor(log2(step)+1)
    to_eval<-unique(sort(c(nmax-step,nmax+step)))
    to_eval<-to_eval[to_eval>=liminf & to_eval<=limsup]
    for(k in running){ #check which ks are already finished
      if(resolved(futureOf(values[[as.character(k)]]))){
        running<-running[!running==k]
        evalued<-c(evalued,k)
      }
    }
    if(length(to_eval)>(ncores-length(running))){
      eval_now<-to_eval[ncores-length(running)]
      eval_latter<-to_eval[-c(1:(ncores-length(running)))]
    }else{
      eval_now<-to_eval
      eval_latter<-c()
    } 
    cat(paste(left_eval," evaluations left:\n",sep=""))
    for(k in eval_now) values[[as.character(k)]] %<-% {testfun(k)}
    running<-c(running,eval_now)
    while(length(eval_now)>0){
      still_running<-eval_now
      for(k in still_running){
        if(resolved(futureOf(values[[as.character(k)]]))){
          thisbics<-(values[[as.character(k)]][[1]])
          thismedian<-median(thisbics)
          #compare k and nmax
          if(significance){
            testeq<-sapply(nmax,function(n){
              maxbics<-values[[as.character(n)]][[1]]
              maxmed<-median(maxbics)
              kt<-kruskal.test(list(thisbics,maxbics))
              if(kt$p.value<=pcrit){
                if(thismedian>maxmed){ 
                  testval=1 
                }else{ 
                  testval=-1 
                }
              }else testval = 0
              return(testval)
            })
          }else{
            testeq<-sapply(nmax,function(n){
              maxbics<-values[[as.character(n)]][[1]]
              maxmed<-median(maxbics)
              if(thismedian>maxmed){ 
                testval=1 
              }else if(thismedian<maxmed){ 
                testval=-1 
              }else testval = 0
              return(testval)
            })
          }
          nmax<-nmax[testeq %in% c(0,-1)] #keep those n better than k
          if( length(nmax)==0 | min(testeq)>=0) nmax<-c(nmax,k) #keep k if no other n is better
          ### end of comparison k and nmax
          eval_now<-eval_now[!eval_now==k]
          running<-running[!running==k]
          evalued<-c(evalued,k)
          if(length(eval_latter)>0){ 
            #start testing another number of signatures:
            nextrun<-eval_latter[1]
            values[[as.character(nextrun)]] %<-% {testfun(nextrun)}
            eval_now<-c(eval_now,nextrun)
            eval_latter<-eval_latter[-1]
          }
        }
      }
    }
  }
  nopt<- min(nmax) #in case of ties, less signatures are better 
  for(k in running){ #check which ks are already finished
    if(resolved(futureOf(values[[as.character(k)]]))){
      evalued<-c(evalued,k)
    }
  }
  cond<-points %in% evalued
  return(list(nopt,points[cond],as.list(values[cond])))
}

