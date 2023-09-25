#### S4 Class SignExp: pair of tensors with signatures and exposures ####
setClass("SignExp",
    slots = c(Sign="array",
        Exp="array",
        samples="character",
        mutations="character",
        signames="character",
        sigSums="matrix",
        normalized="logical",
        Psummary="array",
        Esummary="array",
        Eoutliers="list"),
    prototype = list(Sign=array(NA,dim=c('i'=1,'n'=1,'k'=1)),
                     Exp=array(NA,dim=c('n'=1,'j'=1,'k'=1)),
        samples=NA_character_,
        mutations=NA_character_,
        signames=NA_character_,
        sigSums=matrix(NA_real_,1,1),
        normalized=FALSE,
        Psummary=array(NA,dim=c('i'=1,'n'=1,'q'=6)),
        Esummary=array(NA,dim=c('n'=1,'j'=1,'q'=6)),
        Eoutliers=list())
)

SignExpConstructor<-function(Ps=NA,Es=NA,samplenames=NA,mutnames=NA,signames=NA){
    if(!(is.na(Ps[1]) | is.na(Es[1]))){
        if(!(is.array(Ps) & is.array(Es))){
            stop("Signatures and exposures must be arrays.")
        }
        dp<-dim(Ps) #[i,n,r]
        de<-dim(Es) #[n,j,r]
        if(!(dp[[3]]==de[[3]] & dp[[2]]==de[[1]])){
            stop("Signatures and exposures are not compatible")
        }
        n<-dp[[2]]
        sSums<-apply(Ps,c(2,3),sum)
        if(all(is.na(samplenames))){
            samplenames<-paste("sample",1:de[[2]],sep="_")
        }
        if(all(is.na(mutnames))){
            bases=c("A","C","G","T")
            mutnames<-paste(rep(bases,each=4,6),rep(c("C","T"),each=48),
                rep(bases,24),">",
                c(rep(bases[-2],each=16),rep(bases[-4],each=16)),
                sep="")
        }
        if(all(is.na(signames))){
            signames<-paste("S",1:dp[[2]],sep="")
        }
        SE<-new("SignExp",Sign=Ps,Exp=Es,samples=samplenames,
            mutations=mutnames,signames=signames,
            sigSums=sSums,normalized=FALSE,
            Psummary=array(NA,dim=c('i'=1,'n'=1,'q'=6)),
            Esummary=array(NA,dim=c('n'=1,'j'=1,'q'=6)),
            Eoutliers=list())
        SE<-Normalize(SE)
    }else{
        SE<-new("SignExp")
    }
    return(SE)
}

setGeneric("setSamples",
    def=function(signexp_obj,names){
        standardGeneric("setSamples")
    }
)
setMethod("setSamples",signature(signexp_obj="SignExp",names="ANY"),
    function(signexp_obj,names){
        attr(signexp_obj,"samples")<-names
        return(signexp_obj)
    }
)

setGeneric("setMutations",
    def=function(signexp_obj,mutations){
        standardGeneric("setMutations")
    }
)
setMethod("setMutations",signature(signexp_obj="SignExp",mutations="ANY"),
    function(signexp_obj,mutations){
        attr(signexp_obj,"mutations")<-mutations
        return(signexp_obj)
    }
)


setGeneric("setSignames",
    def=function(signexp_obj,signames){
        standardGeneric("setSignames")
    }
)
setMethod("setSignames",signature(signexp_obj="SignExp",signames="ANY"),
    function(signexp_obj,signames){
        attr(signexp_obj,"signames")<-signames
        return(signexp_obj)
    }
)

setGeneric("Normalize",
    def=function(signexp_obj){
        standardGeneric("Normalize")
    }
)
setMethod("Normalize",signature("SignExp"), function(signexp_obj){
    #Normalize signatures.
    if(!signexp_obj@normalized){
        Ps<-signexp_obj@Sign
        Es<-signexp_obj@Exp
        dp<-dim(Ps)
        de<-dim(Es)
        i<-dp[[1]]
        n<-dp[[2]]
        j<-de[[2]]
        r<-de[[3]]
        for(k in 1:r){
            P<-Ps[,,k]
            E<-Es[,,k]
            vm<-signexp_obj@sigSums[,k,drop=TRUE]
            signexp_obj@Sign[,,k] <- t(t(matrix(as.vector(P),i,n))/vm)
            signexp_obj@Exp[,,k] <- matrix(as.vector(E),n,j)*vm
        }
        signexp_obj@normalized=TRUE
        Psum<-array(NA,dim=c('i'=i,'n'=n,'q'=6))
        Esum<-array(NA,dim=c('n'=n,'j'=j,'q'=5))
        Eout<-list()
        for(k in 1:n){
            P<-signexp_obj@Sign[,k,,drop=TRUE]
            Psum[,k,]<-t(apply(P,1,quantile,c(0.05,0.25,0.50,0.75,0.95,1)))
            E<- signexp_obj@Exp[k,,,drop=TRUE]
            if (j==1){
              v<-as.vector(E)
              bp<-boxplot(v,plot=FALSE)
              bpstats<-list(list(bp$stats,bp$out))
            }else{
              bpstats<-apply(E,1,function(v){
                bp<-boxplot(v,plot=FALSE)
                return(list(bp$stats,bp$out))
              })
            }
            Esum[k,,]<-t(matrix(unlist(lapply(bpstats,function(l){
                l[[1]]})),5,j))
            Eout[[k]]<-lapply(bpstats,function(l){l[[2]]})
        }
        signexp_obj@Psummary <- Psum
        signexp_obj@Esummary <- Esum
        signexp_obj@Eoutliers <- Eout
    }
    return(signexp_obj)
})

setGeneric("Reorder_signatures",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_signatures")
    }
)
setMethod("Reorder_signatures",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change signatures order.
        if(length(ord)==dim(signexp_obj@Sign)[[2]]){
            signexp_obj@Sign<-signexp_obj@Sign[,ord,]
            signexp_obj@Exp<-signexp_obj@Exp[ord,,]
            if(!all(is.na(signexp_obj@sigSums))){
                signexp_obj@sigSums<-signexp_obj@sigSums[ord,]
            }
            if(signexp_obj@normalized){
                signexp_obj@Psummary <- signexp_obj@Psummary[,ord,]
                signexp_obj@Esummary <- signexp_obj@Esummary[ord,,]
                signexp_obj@Eoutliers <- signexp_obj@Eoutliers[ord]
            }
            return(signexp_obj)
        }else{
            stop(paste("'ord' needs to be a vector of length",
                "equal to the number of signatures.",sep=" "))
        }
    }
)

setGeneric("Reorder_samples",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_samples")
    }
)
setMethod("Reorder_samples",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change sample order or take subsets.
        if(!length(ord)==dim(signexp_obj@Sign)[[1]]){
            warning("Reorder_samples will generate a new SignExp object",
                " with the sample subset enumerated in 'ord'.\n")
        }
        signexp_obj@Exp<-signexp_obj@Exp[,ord,,drop=FALSE]
        if(!all(is.na(signexp_obj@samples))){
            signexp_obj@samples<-signexp_obj@samples[ord]
        }
        if(signexp_obj@normalized){
            signexp_obj@Esummary <- signexp_obj@Esummary[,ord,,drop=FALSE]
            for(k in 1:length(signexp_obj@Eoutliers)){
                signexp_obj@Eoutliers[[k]] <- signexp_obj@Eoutliers[[k]][ord]
            }
        }
        return(signexp_obj)
    }
)

setGeneric("Reorder_mutations",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_mutations")
    }
)
setMethod("Reorder_mutations",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change mutation order.
        if(length(ord)==dim(signexp_obj@Sign)[[1]]){
            signexp_obj@Sign<-signexp_obj@Sign[ord,,]
            if(!all(is.na(signexp_obj@mutations))){
                signexp_obj@mutations<-signexp_obj@mutations[ord]
            }
            if(signexp_obj@normalized){
                signexp_obj@Psummary <- signexp_obj@Psummary[ord,,]
            }
            return(signexp_obj)
        }else{
            stop(paste("'ord' needs to be a vector of length",
                "equal to the number of mutations.",sep=" "))
        }
    }
)

setGeneric("Average_sign",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Average_sign")
    }
)
setMethod("Average_sign",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        Ps<-signexp_obj@Sign #[i,n,r]
        Phat<-apply(Ps,c(1,2),mean)
        rownames(Phat)<-signexp_obj@mutations
        colnames(Phat)<-signexp_obj@signames
        return(Phat)
    }
)

setGeneric("Median_sign",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Median_sign")
    }
)
setMethod("Median_sign",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        i<-dp[[1]]; n<-dp[[2]]; r<-dp[[3]]
        Phat<-signexp_obj@Psummary[,,3,drop=TRUE]
        if(n==1) Phat<-matrix(as.vector(Phat),i,n)
        rownames(Phat)<-signexp_obj@mutations
        colnames(Phat)<-signexp_obj@signames
        return(Phat)
    }
)

setGeneric("Average_exp",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Average_exp")
    }
)
setMethod("Average_exp",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        Es<-signexp_obj@Exp #[n,j,r]
        Ehat<-apply(Es,c(1,2),mean)
        rownames(Ehat)<-signexp_obj@signames
        colnames(Ehat)<-signexp_obj@samples
        return(Ehat)
    }
)

setGeneric("Median_exp",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Median_exp")
    }
)
setMethod("Median_exp",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        de <- dim(signexp_obj@Exp) #[n,j,r]
        n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
        Ehat<-signexp_obj@Esummary[,,3,drop=TRUE]
        if(n==1 | j==1) Ehat<-matrix(as.vector(Ehat),n,j)
        rownames(Ehat)<-signexp_obj@signames
        colnames(Ehat)<-signexp_obj@samples
        return(Ehat)
    }
)

#Update signatures and expositions
setGeneric("AddSamples",
           def=function(Signatures, originalCounts,originalOpp, NewCounts, 
                        NewOpp, updateSigs, updateSigNum, BICsignificance, 
                        critical_p, ap, bp, ae, be, lp, le, var.ap, var.ae, 
                        burn_it, eval_it, main_burn, main_eval, estimate_hyper, 
                        EMit_lim, EM_eval, parallel, ncore, samplerGap, 
                        min_samples_per_core){
             standardGeneric("AddSamples")
           }
)

setMethod("AddSamples",signature(Signatures="ANY", originalCounts="ANY", 
                                 originalOpp="ANY", NewCounts="matrix", 
                                 NewOpp="ANY", updateSigs="logical",
                                 updateSigNum="logical", 
                                 BICsignificance="logical", 
                                 critical_p = "numeric", 
                                 ap="ANY", bp="ANY", ae="ANY", be="ANY", 
                                 lp="ANY", le="ANY", 
                                 var.ap="ANY",var.ae="ANY", 
                                 burn_it="numeric", eval_it="numeric", 
                                 main_burn="numeric", main_eval="numeric",
                                 estimate_hyper="logical", EMit_lim="numeric", 
                                 EM_eval="numeric", parallel="logical",
                                 ncore="ANY", samplerGap="numeric", 
                                 min_samples_per_core="numeric"),
          function(Signatures, originalCounts,originalOpp, NewCounts, NewOpp, 
                   updateSigs=FALSE, updateSigNum=FALSE, BICsignificance=FALSE, 
                   critical_p = 0.05, ap=NA,bp=NA,ae=NA,be=NA,lp=NA,le=NA,
                   var.ap=10, var.ae=10, burn_it=1000, eval_it=1000, 
                   main_burn=10000, main_eval=2000, estimate_hyper=FALSE,
                   EMit_lim=100, EM_eval=100, parallel=TRUE,ncore=NA,
                   samplerGap=10, min_samples_per_core=200){
            if(is.list(Signatures)){
              Nsign<-Signatures@Nsign
              tested_n<-Signatures@tested_n
              HH<-Signatures@Hyper_paths
              Signatures<-Signatures@SignExposures
              BestHH<-HH[[which(tested_n==Nsign)]]
              best_hyperparam<-BestHH[NROW(BestHH),]
              if(is.na(ap)) ap<-best_hyperparam$ap
              if(is.na(bp)) bp<-best_hyperparam$bp
              if(is.na(ae)) ae<-best_hyperparam$ae
              if(is.na(be)) be<-best_hyperparam$be
              if(is.na(lp)) lp<-best_hyperparam$lp
              if(is.na(le)) le<-best_hyperparam$le
              eh<-estimate_hyper
            }else if(class(Signatures)=="SignExp"){
              if(is.na(ap)) ap<-1
              if(is.na(bp)) bp<-1
              if(is.na(ae)) ae<-1
              if(is.na(be)) be<-1
              if(is.na(lp)) lp<-1
              if(is.na(le)) le<-1
              eh<-TRUE
            }else{ 
              stop("Signatures must be a signeR output or a SignExp object.\n")
            }
            if(!Signatures@normalized) Signatures<-Normalize(Signatures)
            dp <- dim(Signatures@Sign) #[i,n,r]
            de <- dim(Signatures@Exp) #[n,j,r]
            i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
            Ps <- Signatures@Sign
            Es <- Signatures@Exp
            P<-Ps[,,r]
            E<-Es[,,r]
            dn<-dim(NewCounts)
            newsamp<-dn[[2]]
            fixP<-!updateSigs
            keep_param<-TRUE ############# FALSE would make no difference
            #Estimate new expositions based on P
            NewExp<-matrix(0,n,newsamp)
            for(k in 1:newsamp){
              found<-NewCounts[,k]
              opp<-NewOpp[,k]
              RSS0<-function(vet){
                expected<-Ps%*%matrix(vet,n,1) * opp
                return( sum( expected+lgamma(found+1)-found*log(expected) ) )
                
              }
              Fit0 <- nloptr(x0=rep(1,n), eval_f = RSS0, lb=rep(0,n), 
                             opts = list(algorithm = "NLOPT_LN_SBPLX", 
                                         xtol_rel=1e-100, xtol_abs=1e-100, 
                                         maxeval = 1e6))
              NewExp[,k]<-as.vector(Fit0$solution)
            }
            #proceed with sampler, keeping or changing P according to updateSigs
            if(!is.na(originalCounts[1])){
              M<-cbind(originalCounts,NewCounts)
              E<-cbind(E,NewExp)
            }else{
              M<-NewCounts
              E<-NewExp
            }            
            if(!is.na(originalOpp[1])){
              W<-cbind(originalOpp,NewOpp)
            }else{
              W<-NewOpp
            }
            #Initial guess for Z
            Z  <- array(0,dim=c('i'=i,'j'=j,'n'=n))
            PE <- P %*% E
            Fi <- array(0,dim=c('i'=i,'j'=j,'n'=n))
            for (m in 1:n){
              Fi[,,'n'=m] <- (P[,m,drop=FALSE] %*% E[m,,drop=FALSE])/PE
            }
            for (s in 1:i){
              Z['i'=s,,] <- t(sapply(1:j,function(g){
                rmultinom(1, size = M[s,g], prob = as.vector(Fi['i'=s,'j'=g,]))
              }))
            }
            rm(PE,Fi)
            if (fixP){ 
              n<-NCOL(P)
              Ap <- matrix(0,i,n)
              Bp <- matrix(0,i,n) 
              apvet<-NA; bpvet<-NA; lpvet<-NA 
            }else{
              Ap <- matrix(rexp(i*n,rate=lp),i,n)
              Bp <- matrix(rgamma(i*n,shape=ap,rate=bp),i,n)
              apvet<-ap; bpvet<-bp; lpvet<-lp 
            }
            Ae <- matrix(rexp(n*j,rate=le),n,j)
            Be <- matrix(rgamma(n*j,shape=ae,rate=be),n,j)
            aevet<-ae; bevet<-be; levet<-le
            if(estimate_hyper){#Hyperparameter optimization
              burn_HO<-burn_it
              upgrade<-100
              it<-0
              cat("EM algorithm:\n")
              progbar <- txtProgressBar(style=3)
              while(upgrade>0.05 & it < EMit_lim){
                SamplesHO <-GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,
                                            ap,bp,ae,be,lp,le,
                                            var.ap,var.ae,burn=burn_HO,
                                            eval=EM_eval, Pfixed=fixP,
                                            Zfixed=FALSE,Thetafixed=FALSE,
                                            Afixed=FALSE,keep_par=TRUE)
                Zs <- SamplesHO[[1]]
                Ps <- SamplesHO[[3]]
                Es <- SamplesHO[[4]]
                Aps <- SamplesHO[[5]]
                Bps <- SamplesHO[[6]]
                Aes <- SamplesHO[[7]]
                Bes <- SamplesHO[[8]]
                Z <- Zs[,,,as.numeric(dim(Zs)[4])]
                dim(Z)<-dim(Zs)[1:3]
                E <- Es[,,as.numeric(dim(Es)[3])]
                dim(E)<-dim(Es)[1:2]
                Ae <- Aes[,,as.numeric(dim(Aes)[3])]
                dim(Ae)<-dim(Aes)[1:2]
                Be <- Bes[,,as.numeric(dim(Bes)[3])]
                dim(Be)<-dim(Bes)[1:2]
                aeold <- ae; beold <- be; leold <- le
                New_param<-EstimateParameters(Aes,Bes,aeold,beold)
                ae<-New_param[[1]]
                be<-New_param[[2]]
                le<-New_param[[3]]
                aevet<-c(aevet,ae)
                bevet<-c(bevet,be)
                levet<-c(levet,le)
                if(fixP){
                  upgrade<-max(abs(c(ae-aeold,be-beold,le-leold)))
                }else{
                  P <- Ps[,,as.numeric(dim(Ps)[3])]
                  dim(P)<-dim(Ps)[1:2]
                  Ap <- Aps[,,as.numeric(dim(Aps)[3])]
                  dim(Ap)<-dim(Aps)[1:2]
                  Bp <- Bps[,,as.numeric(dim(Bps)[3])]
                  dim(Bp)<-dim(Bps)[1:2]
                  apold <- ap; bpold <- bp; lpold <- lp 
                  New_param<-EstimateParameters(Aps,Bps,apold,bpold)
                  ap<-New_param[[1]]
                  bp<-New_param[[2]]
                  lp<-New_param[[3]]
                  apvet<-c(apvet,ap)
                  bpvet<-c(bpvet,bp)
                  lpvet<-c(lpvet,lp)
                  upgrade<-max(abs(c(ap-apold,bp-bpold,ae-aeold,
                                     be-beold,lp-lpold,le-leold)))
                }                
                it <- it+1
                burn_HO<-200
                setTxtProgressBar(progbar, it/EMit_lim)
              }
              setTxtProgressBar(progbar, 1)
              cat("\n")
            }
            testedn<-c(n)
            HH<-list(data.frame('ap'=apvet,'bp'=bpvet,'ae'=aevet,
                                'be'=bevet,'lp'=lpvet,'le'=levet))
            #if updateSigNum, test sampler with one more signature or one less?.
            if(updateSigs & updateSigNum){
              #Sampler with original sigs
              if(n==1){
                cat("Running  Gibbs sampler for 1 signature...")
              }else{
                cat(paste0("Running  Gibbs sampler for ",n," signatures...",
                           collapse=""))
              }
              Samples <- GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,
                                         ap,bp,ae,be,lp,le,
                                         var.ap,var.ae,burn=burn_it,
                                         eval=eval_it, Pfixed=fixP,
                                         Zfixed=FALSE,Thetafixed=FALSE,
                                         Afixed=FALSE, keep_par=keep_param)
              if(keep_param){
                Zs <- Samples[[1]] 
                if (!fixP){
                  Aps <- Samples[[5]] 
                  Bps <- Samples[[6]] 
                }
                Aes <- Samples[[7]]
                Bes <- Samples[[8]]
              }
              Ps <- Samples[[3]]
              Es <- Samples[[4]]
              rm(Samples)
              # BICs #
              lgM<-lgamma(M+1)
              loglikes<-sapply(1:eval_it,function(r){
                PE <- (matrix(as.vector(Ps[,,'r'=r]),i,n) %*% 
                         matrix(as.vector(Es[,,'r'=r]),n,j))*W
                return(sum(M*log(PE)-lgM-PE))
              }) #loglike is the sum of log-densities of Poisson distributions
              Bics <- 2*loglikes -n*(i+j)*log(j)
              testBics<-list(Bics)
              update<-TRUE
              while(update){
                Expect<-Phat%*%Ehat
                Mres<-apply(M-Expect,c(1,2),function(x){max(x,0)})
                Extra_run<-eBayesNMF(Mres,W,n=1,ap,bp,ae,be,lp,le,
                                     var.ap,var.ae,burn_it=burn_it,
                                     eval_it=eval_it, estimate_hyper=FALSE)
                ExtraP<-Extra_run[[1]]
                ExtraE<-Extra_run[[2]]
                P2<-cbind(Phat,ExtraP)
                E2<-rbind(Ehat,ExtraE)
                i=NROW(P2); j=NCOL(E2); n2=NCOL(P2)
                #Initial guess for Z
                Z  <- array(0,dim=c('i'=i,'j'=j,'n'=n2))
                PE <- P2 %*% E2
                Fi <- array(0,dim=c('i'=i,'j'=j,'n'=n2))
                for (m in 1:n2){
                  Fi[,,'n'=m] <- (P2[,m,drop=FALSE] %*% E2[m,,drop=FALSE])/PE
                }
                for (s in 1:i){
                  Z['i'=s,,] <- t(sapply(1:j,function(g){
                    rmultinom(1,size=M[s,g],prob=as.vector(Fi['i'=s,'j'=g,]))
                  }))
                }
                rm(PE,Fi)
                Ap2 <- matrix(rexp(i*n2,rate=lp),i,n2)
                Bp2 <- matrix(rgamma(i*n2,shape=ap,rate=bp),i,n2)
                apvet<-ap; bpvet<-bp; lpvet<-lp 
                Ae2 <- matrix(rexp(n2*j,rate=le),n2,j)
                Be2 <- matrix(rgamma(n2*j,shape=ae,rate=be),n2,j)
                aevet<-ae; bevet<-be; levet<-le
                if(estimate_hyper){#Hyperparameter optimization
                  burn_HO<-burn_it
                  upgrade<-100
                  it<-0
                  cat("EM algorithm:\n")
                  progbar <- txtProgressBar(style=3)
                  while(upgrade>0.05 & it < EMit_lim){
                    SamplesHO <-GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,
                                                ap,bp,ae,be,lp,le,
                                                var.ap,var.ae,burn=burn_HO,
                                                eval=EM_eval, Pfixed=fixP,
                                                Zfixed=FALSE,Thetafixed=FALSE,
                                                Afixed=FALSE,keep_par=TRUE)
                    Zs <- SamplesHO[[1]]
                    Ps <- SamplesHO[[3]]
                    Es <- SamplesHO[[4]]
                    Aps <- SamplesHO[[5]]
                    Bps <- SamplesHO[[6]]
                    Aes <- SamplesHO[[7]]
                    Bes <- SamplesHO[[8]]
                    Z <- Zs[,,,as.numeric(dim(Zs)[4])]
                    dim(Z)<-dim(Zs)[1:3]
                    E <- Es[,,as.numeric(dim(Es)[3])]
                    dim(E)<-dim(Es)[1:2]
                    Ae <- Aes[,,as.numeric(dim(Aes)[3])]
                    dim(Ae)<-dim(Aes)[1:2]
                    Be <- Bes[,,as.numeric(dim(Bes)[3])]
                    dim(Be)<-dim(Bes)[1:2]
                    aeold <- ae; beold <- be; leold <- le
                    New_param<-EstimateParameters(Aes,Bes,aeold,beold)
                    ae<-New_param[[1]]
                    be<-New_param[[2]]
                    le<-New_param[[3]]
                    aevet<-c(aevet,ae)
                    bevet<-c(bevet,be)
                    levet<-c(levet,le)
                    if(fixP){
                      upgrade<-max(abs(c(ae-aeold,be-beold,le-leold)))
                    }else{
                      P <- Ps[,,as.numeric(dim(Ps)[3])]
                      dim(P)<-dim(Ps)[1:2]
                      Ap <- Aps[,,as.numeric(dim(Aps)[3])]
                      dim(Ap)<-dim(Aps)[1:2]
                      Bp <- Bps[,,as.numeric(dim(Bps)[3])]
                      dim(Bp)<-dim(Bps)[1:2]
                      apold <- ap; bpold <- bp; lpold <- lp 
                      New_param<-EstimateParameters(Aps,Bps,apold,bpold)
                      ap<-New_param[[1]]
                      bp<-New_param[[2]]
                      lp<-New_param[[3]]
                      apvet<-c(apvet,ap)
                      bpvet<-c(bpvet,bp)
                      lpvet<-c(lpvet,lp)
                      upgrade<-max(abs(c(ap-apold,bp-bpold,ae-aeold,
                                         be-beold,lp-lpold,le-leold)))
                    }                
                    it <- it+1
                    burn_HO<-200
                    setTxtProgressBar(progbar, it/EMit_lim)
                  }
                  setTxtProgressBar(progbar, 1)
                  cat("\n")
                }
                if(n2==1){
                  cat("Running  Gibbs sampler for 1 signature...")
                }else{
                  cat(paste0("Running  Gibbs sampler for ",n2," signatures...",
                             collapse=""))
                }
                Samples <- GibbsSamplerCpp(M,W,Z,P2,E2,Ap2,Bp2,Ae2,Be2,
                                           ap,bp,ae,be,lp,le,
                                           var.ap,var.ae,burn=burn_it,
                                           eval=eval_it, Pfixed=FALSE,
                                           Zfixed=FALSE,Thetafixed=FALSE,
                                           Afixed=FALSE, keep_par=keep_param)
                Ps2 <- Samples[[3]]
                Es2 <- Samples[[4]]
                rm(Samples)
                # BICs #
                lgM<-lgamma(M+1)
                loglikes<-sapply(1:eval_it,function(r){
                  PE <- (matrix(as.vector(Ps2[,,'r'=r]),i,n2) %*% 
                           matrix(as.vector(Es2[,,'r'=r]),n2,j))*W
                  return(sum(M*log(PE)-lgM-PE))
                }) #loglike is the sum of log-densities of Poisson distributions
                Bics2 <- 2*loglikes -n2*(i+j)*log(j)
                testedn<-c(testedn,n2)
                HH<-c(HH,list(data.frame('ap'=apvet,'bp'=bpvet,'ae'=aevet,
                                         'be'=bevet,'lp'=lpvet,'le'=levet)))
                testBics<-c(testBics,list(Bics2))
                if(BICsignificance){ #Compare Bics
                  kt<-kruskal.test(list(Bics,Bics2))
                  sig<-kt$p.value>=critical_p
                }else{
                  sig<-TRUE
                }
                update <- sig & median(Bics2)>median(Bics) 
                if(update){
                  n<-n2
                  P<-P2
                  E<-E2
                  Ap<-Ap2
                  Ae<-Ae2
                  Bp<-Bp2
                  Be<-Be2
                  #ap,bp,ae,be,le,lp are updated already.
                  SE<-SignExpConstructor(Ps,Es)
                  Phat<-Median_sign(SE)
                  Ehat<-Median_exp(SE)
                }
              }
            }
            ################# Final run #################
            if(n==1){
              cat("Running  Gibbs sampler for 1 signature...")
            }else{
              cat(paste0("Running  Gibbs sampler for ",n," signatures...",
                         collapse=""))
            }
            Samples <- GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,
                                       ap,bp,ae,be,lp,le,
                                       var.ap,var.ae,burn=main_burn,
                                       eval=main_eval, Pfixed=fixP,
                                       Zfixed=FALSE,Thetafixed=FALSE,
                                       Afixed=FALSE, keep_par=keep_param)
            if(keep_param){
              Zs <- Samples[[1]] 
              if (!fixP){
                Aps <- Samples[[5]] 
                Bps <- Samples[[6]] 
              }
              Aes <- Samples[[7]]
              Bes <- Samples[[8]]
            }
            Ps <- Samples[[3]]
            Es <- Samples[[4]]
            rm(Samples)
            # BICs #
            lgM<-lgamma(M+1)
            loglikes<-sapply(1:eval_it,function(r){
              PE <- (matrix(as.vector(Ps[,,'r'=r]),i,n) %*% 
                       matrix(as.vector(Es[,,'r'=r]),n,j))*W
              return(sum(M*log(PE)-lgM-PE))
            }) #loglike is the sum of log-densities of Poisson distributions
            Bics <- 2*loglikes -n*(i+j)*log(j)
            SE<-SignExpConstructor(Ps,Es)
            Phat<-Median_sign(SE)
            Ehat<-Median_exp(SE)
            if(!fixP){ #Ordering by total exposure.
              totalexp<-rowSums(Ehat)*colSums(Phat) 
              ord<-order(totalexp,decreasing=TRUE)
              Ps<-Ps[,ord,,drop=FALSE]
              Es<-Es[ord,,,drop=FALSE]
              SE<-SignExpConstructor(Ps,Es)
              Phat<-Median_sign(SE)
              Ehat<-Median_exp(SE)
            }
            if(!(updateSigs & updateSigNum)){
              testBics<-list(Bics)
            }
            # Output #
            if(is.list(Signatures)){
              result<-list(Nsign=n,
                           tested_n=testedn,
                           Test_BICs=testBics,
                           Phat=Phat,
                           Ehat=Ehat,
                           SignExposures=SE,
                           BICs=Bics,
                           Hyper_paths=HH)
            }else{
              result<-SE
            }
            return(result)
          }
)
