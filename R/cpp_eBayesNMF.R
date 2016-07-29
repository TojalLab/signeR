############################################################################
# Parameters:
# M: matrix of mutation counts, each column corresponds to one sample 
# W: matrix of mutation opportunities, each column corresponds to one sample
# n: number of signatures
# ap, bp: hyperparameters of rate parameters used to generate signatures
# ae, be: hyperparameters of rate parameters used to generate exposures
# lp: hyperparameters of shape parameters used to generate signatures
# le: hyperparameters of shape parameters used to generate exposures
# var.ap: variance of the gamma distribution used to generate proposals 
#         for shape parameters of signatures
# var.ae: variance of the gamma distribution used to generate proposals 
#         for shape parameters of exposures
# burn_it: number of burning iterations of the Gibbs sampler
# eval_it: number of iterations of the Gibbs sampler
# start: NMF algorithm used to generate initial values for signatures and 
#        exposures, options: "brunet","KL","lee","Frobenius","offset",
#                    "nsNMF","ls-nmf","pe-nmf","siNMF","snmf/r" or "snmf/l".
# estimate_hyper: if TRUE, algorithm estimates optimal values of 
#                 ap,bp,ae,be,lp,le. Start values are still required.
# EM_lim: limit of EM iterations for the estimation of ap,bp,ae,be,lp,le.
# keep_param: if TRUE, algorithm export hyperparameters
############################################################################

eBayesNMF<-function(M,W,n,ap,bp,ae,be,lp,le,
    var.ap=10,var.ae=10,
    burn_it=500,eval_it=1000, EM_eval_it=100,
    start='random',estimate_hyper=FALSE,
    EM_lim=100, keep_param=FALSE){
    i<-dim(M)[1]; j<-dim(M)[2]
    Bp <- matrix(rgamma(i*n,shape=ap,rate=bp),i,n)
    Be <- matrix(rgamma(n*j,shape=ae,rate=be),n,j)
    Ap <- matrix(rexp(i*n,rate=lp),i,n)
    Ae <- matrix(rexp(n*j,rate=lp),n,j)
    if (start=="random"){
        P = matrix(0,i,n); E = matrix(0,n,j)
        for(m in 1:n){
            for(s in 1:i){ P[s,m] = rgamma(1, shape=Ap[s,m]+1, rate=Bp[s,m] ) }
            for(g in 1:j){ E[m,g] = rgamma(1, shape=Ae[m,g]+1, rate=Be[m,g] ) }
        }
    }else{
        Mcor<-M
        if (min(colSums(Mcor))==0 | min(rowSums(Mcor))==0){Mcor[Mcor==0]<-0.01}
        Wcor<-W
        Wcor[Wcor==0] <- 0.01
        nnegfact = nmf(Mcor/Wcor,n,method=start) 
        #start can be the name of one NMF algorithm
        P = basis(nnegfact)
        E = coef(nnegfact)
        rm(nnegfact)
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
    apvet<-ap; bpvet<-bp; aevet<-ae; bevet<-be; lpvet<-lp; levet<-le
    if(estimate_hyper){#Hyperparameter optimization
        burn_HO<-burn_it
        upgrade<-100
        it<-0
        cat("EM algorithm:\n")
        progbar <- txtProgressBar(style=3)
        while(upgrade>0.05 & it < EM_lim){
            SamplesHO <-GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,ap,bp,ae,be,lp,le,
                var.ap,var.ae,burn=burn_HO,eval=EM_eval_it,
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
            P <- Ps[,,as.numeric(dim(Ps)[3])]
            dim(P)<-dim(Ps)[1:2]
            E <- Es[,,as.numeric(dim(Es)[3])]
            dim(E)<-dim(Es)[1:2]
            Ap <- Aps[,,as.numeric(dim(Aps)[3])]
            dim(Ap)<-dim(Aps)[1:2]
            Bp <- Bps[,,as.numeric(dim(Bps)[3])]
            dim(Bp)<-dim(Bps)[1:2]
            Ae <- Aes[,,as.numeric(dim(Aes)[3])]
            dim(Ae)<-dim(Aes)[1:2]
            Be <- Bes[,,as.numeric(dim(Bes)[3])]
            dim(Be)<-dim(Bes)[1:2]
            apold <- ap; bpold <- bp; aeold <- ae 
            beold <- be; lpold <- lp; leold <- le
            New_param<-EstimateParameters(Aps,Bps,Aes,Bes,
                apold,bpold,aeold,beold)
            ap<-New_param[[1]]
            bp<-New_param[[2]]
            ae<-New_param[[3]]
            be<-New_param[[4]]
            lp<-New_param[[5]]
            le<-New_param[[6]]
            apvet<-c(apvet,ap)
            bpvet<-c(bpvet,bp)
            aevet<-c(aevet,ae)
            bevet<-c(bevet,be)
            lpvet<-c(lpvet,lp)
            levet<-c(levet,le)
            upgrade<-max(abs(c(ap-apold,bp-bpold,ae-aeold,
                be-beold,lp-lpold,le-leold)))
            it <- it+1
            burn_HO<-200
            setTxtProgressBar(progbar, it/EM_lim)
        }
        setTxtProgressBar(progbar, 1)
        cat("\n")
    }
    #Sampler
    if(n==1){
        cat("Running  Gibbs sampler for 1 signature...")
    }else{
        cat(paste0("Running  Gibbs sampler for ",n," signatures...",
            collapse=""))
    }
    Samples <- GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,ap,bp,ae,be,lp,le,
        var.ap,var.ae,burn=burn_it,eval=eval_it,
        Zfixed=FALSE,Thetafixed=FALSE,Afixed=FALSE,
        keep_par=keep_param)
    if(keep_param){
        Zs <- Samples[[1]] 
        Aps <- Samples[[5]] 
        Bps <- Samples[[6]] 
        Aes <- Samples[[7]]
        Bes <- Samples[[8]]
    }
    Ps <- Samples[[3]]
    Es <- Samples[[4]]
    rm(Samples)
    Phat<-apply(Ps,c(1,2),median)
    Ehat<-apply(Es,c(1,2),median)
    # BICs #
    lgM<-lgamma(M+1)
    loglikes<-sapply(1:eval_it,function(r){
        PE <- (matrix(as.vector(Ps[,,'r'=r]),i,n) %*% 
                matrix(as.vector(Es[,,'r'=r]),n,j))*W
        return(sum(M*log(PE)-lgM-PE))
    }) #loglike is the sum of log-densities of Poisson distributions
    Bics <- 2*loglikes -n*(i+j)*log(j)
    totalexp<-rowSums(Ehat)*colSums(Phat) #Ordering by total exposure.
    ord<-order(totalexp,decreasing=TRUE)
    Ps<-Ps[,ord,,drop=FALSE]
    Phat<-Phat[,ord,drop=FALSE]
    Es<-Es[ord,,,drop=FALSE]
    Ehat<-Ehat[ord,,drop=FALSE]
    # Output #
    SE<-SignExpConstructor(Ps,Es)
    out_list<-list(Phat,Ehat,Bics,SE)
    if(keep_param){ out_list<-c(out_list,list(Aps,Bps,Aes,Bes)) }
    if(estimate_hyper){
        out_list[[9]]<-data.frame('ap'=apvet,'bp'=bpvet,'ae'=aevet,
            'be'=bevet,'lp'=lpvet,'le'=levet)
    }else{
        out_list[[9]]<-data.frame('ap'=ap,'bp'=bp,'ae'=ae,
            'be'=be,'lp'=lp,'le'=le)
    }
    cat("Done.\n")
    return(out_list)
}
