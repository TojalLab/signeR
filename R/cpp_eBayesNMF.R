library(tensorA)
library(NMF)
library(nloptr)
  
eBayesNMF<-function(M,W,n,ap,bp,ae,be,lp,le,
                    var.ap=10,var.ae=10,
                    burn_it=500,eval_it=1000, EM_eval_it=100,
                    start='random',estimate_hyper=FALSE,
                    EM_lim=100, keep_param=FALSE){
    ####################################################################################################################
    # Parameters:
    # M: matrix of mutation counts, each column corresponds to one sample (96 x nsamples)
    # W: matrix of mutation opportunities, each column corresponds to one sample (96 x nsamples)
    # n: number of signatures
    # ap, bp: hyperparameters of rate parameters used to generate signatures (gamma distributions)
    # ae, be: hyperparameters of rate parameters used to generate exposures (gamma distributions)
    # lp: hyperparameters of shape parameters used to generate signatures (exponential distributions)
    # le: hyperparameters of shape parameters used to generate exposures (exponential distributions)
    # var.ap: variance of the gamma distribution used to generate proposals for shape parameters of signatures
    # var.ae: variance of the gamma distribution used to generate proposals for shape parameters of exposures
    # burn_it: number of burning iterations of the Gibbs sampler
    # eval_it: number of iterations of the Gibbs sampler
    # start: NMF algorithm used to generate initial values for signatures and exposures, 
    #        options: "brunet","KL","lee","Frobenius","offset","nsNMF","ls-nmf","pe-nmf","siNMF","snmf/r" or "snmf/l".
    # estimate_hyper: if TRUE, algorithm estimates optimal values of ap,bp,ae,be,lp,le. Start values are still required.
    # EM_lim: limit of EM iterations for the estimation of ap,bp,ae,be,lp,le.
    # keep_param: if TRUE, algorithm export hyperparameters
    ####################################################################################################################
    #Start
    i<-dim(M)[1]
    j<-dim(M)[2]
    Bp <- matrix(rgamma(i*n,shape=ap,rate=bp),i,n)
    Be <- matrix(rgamma(n*j,shape=ae,rate=be),n,j)
    Ap <- matrix(rexp(i*n,rate=lp),i,n)
    Ae <- matrix(rexp(n*j,rate=lp),n,j)
    if (start=="random"){
      P = matrix(0,i,n)
      E = matrix(0,n,j)
      for(m in 1:n){
        for(s in 1:i){
          P[s,m] = rgamma(1, shape=Ap[s,m]+1, rate=Bp[s,m] )
        }
        for(g in 1:j){
          E[m,g] = rgamma(1, shape=Ae[m,g]+1, rate=Be[m,g] )
        }
      }
    }else{
      Mcor<-M
      if (min(colSums(Mcor))==0 | min(rowSums(Mcor))==0){
        Mcor[Mcor==0] <- 0.01
      }
      Wcor<-W
      Wcor[Wcor==0] <- 0.01
      nnegfact = nmf(Mcor/Wcor,n,method=start) 
      #start can be the name of one NMF algorithm ("brunet","KL","lee","Frobenius","offset","nsNMF","ls-nmf","pe-nmf","siNMF","snmf/r","snmf/l")
      P = basis(nnegfact)
      E = coef(nnegfact)
      rm(nnegfact)
    }
    #Initial guess for Z
    Z <- to.tensor(0,c('i'=i,'j'=j,'n'=n))
    PE <- P %*% E
    Fi <- to.tensor(0,c('i'=i,'j'=j,'n'=n))
    for (m in 1:n){
      Fi[,,'n'=m] <- (P[,m,drop=FALSE] %*% E[m,,drop=FALSE])/PE
    }
    for (s in 1:i){
      Z['i'=s,,] <- t(sapply(1:j,function(g){
          rmultinom(1, size = M[s,g], prob = as.vector(Fi['i'=s,'j'=g,]))
      }))
    }
    rm(PE,Fi)
    apvet<-ap
    bpvet<-bp
    aevet<-ae
    bevet<-be
    lpvet<-lp
    levet<-le
    if(estimate_hyper){#Hyperparameter optimization
      burn_HO<-burn_it
      upgrade<-100
      it<-0
      cat("EM algorithm:\n")
      progbar <- txtProgressBar(style=3)
      while(upgrade>0.05 & it < EM_lim){
        SamplesHO <- GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,ap,bp,ae,be,lp,le,var.ap,var.ae,
                                     burn=burn_HO,eval=EM_eval_it,Zfixed=FALSE,Thetafixed=FALSE,Afixed=FALSE,keep_par=TRUE)
        Zs <- to.tensor(SamplesHO[[1]])
        Ps <- to.tensor(SamplesHO[[3]])
        Es <- to.tensor(SamplesHO[[4]])
        Aps <- to.tensor(SamplesHO[[5]])
        Bps <- to.tensor(SamplesHO[[6]])
        Aes <- to.tensor(SamplesHO[[7]])
        Bes <- to.tensor(SamplesHO[[8]])
        Z <- Zs[,,,as.numeric(dim(Zs)[4])]
        dim(Z)<-dim(Zs)[1:3]
        P <- matrix(as.vector(Ps[,,as.numeric(dim(Ps)[3])]),dim(Ps)[1],dim(Ps)[2])
        E <- matrix(as.vector(Es[,,as.numeric(dim(Es)[3])]),dim(Es)[1],dim(Es)[2])
        Ap <- matrix(as.vector(Aps[,,as.numeric(dim(Aps)[3])]),dim(Aps)[1],dim(Aps)[2])
        Bp <- matrix(as.vector(Bps[,,as.numeric(dim(Bps)[3])]),dim(Bps)[1],dim(Bps)[2])
        Ae <- matrix(as.vector(Aes[,,as.numeric(dim(Aes)[3])]),dim(Aes)[1],dim(Aes)[2])
        Be <- matrix(as.vector(Bes[,,as.numeric(dim(Bes)[3])]),dim(Bes)[1],dim(Bes)[2])
        apold <- ap
        bpold <- bp
        aeold <- ae
        beold <- be
        lpold <- lp
        leold <- le
        New_param<-Param_Estimate(Aps,Bps,Aes,Bes,apold,bpold,aeold,beold)
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
        upgrade<-max(abs(c(ap-apold,bp-bpold,ae-aeold,be-beold,lp-lpold,le-leold)))
        it <- it+1
        burn_HO<-200
        setTxtProgressBar(progbar, it/EM_lim)
      }
      setTxtProgressBar(progbar, 1)
      cat("\n")
      #cat("Optimal Parameters are: ",ap,",",bp,",",ae,",",be,",",lp,",",le,".\n")
      #cat("Obtained after ",it," iterations.\n")
    }
    #Sampler
    if(n==1){
      cat("Running  Gibbs sampler for 1 signature...")
    }else{
      cat(paste0("Running  Gibbs sampler for ",n," signatures...",collapse=""))
    }
    Samples <- GibbsSamplerCpp(M,W,Z,P,E,Ap,Bp,Ae,Be,ap,bp,ae,be,lp,le,var.ap,var.ae,
                            burn=burn_it,eval=eval_it,Zfixed=FALSE,Thetafixed=FALSE,Afixed=FALSE,keep_par=keep_param)
               #Samples = list(Zs,Fis,Ps,Es,Aps,Bps,Aes,Bes) #,Zstar,Pstar,Estar,Apstar,Bpstar,Aestar,Bestar)
    if(keep_param){
      Zs <- to.tensor(Samples[[1]])
      Aps <- to.tensor(Samples[[5]])
      Bps <- to.tensor(Samples[[6]])
      Aes <- to.tensor(Samples[[7]])
      Bes <- to.tensor(Samples[[8]])
    }
    Ps <- to.tensor(Samples[[3]])
    Es <- to.tensor(Samples[[4]])
    #cat("GibbsSampler finished.\n")
    rm(Samples)
    Phat<-matrix(0,i,n)
    Ehat<-matrix(0,n,j)
    for(k in 1:n){
      Pcut<-matrix(as.vector(Ps[,k,]),i,eval_it)
      Phat[,k]<-apply(Pcut,1,median)
      Ecut<-matrix(as.vector(Es[k,,]),j,eval_it)
      Ehat[k,]<-apply(Ecut,1,median)
    }
    #cat('Estimators obtained.\n')
    ######### BICs calculation############
    lgM<-lgamma(M+1)
    loglikes<-sapply(1:eval_it,function(r){
      PE <- (matrix(as.vector(Ps[,,'r'=r]),i,n) %*% matrix(as.vector(Es[,,'r'=r]),n,j))*W
      return(sum(M*log(PE)-lgM-PE))
    })
    #loglike is the sum of log-densities of Poisson distributions
    Bics <- 2*loglikes -n*(i+j)*log(j)
    ############################
    #Ordering by total exposure.
    totalexp<-rowSums(Ehat)*colSums(Phat)
    ord<-order(totalexp,decreasing=TRUE)
    Ps<-Ps[,ord,,drop=FALSE]
    Phat<-Phat[,ord,drop=FALSE]
    Es<-Es[ord,,,drop=FALSE]
    Ehat<-Ehat[ord,,drop=FALSE]
    #Output
    SE<-SignExp(Ps,Es)
    out_list<-list(Phat,Ehat,Bics,SE)
    if(keep_param){
      out_list<-c(out_list,list(Aps,Bps,Aes,Bes))
    }
    if(estimate_hyper){
      out_list[[9]]<-data.frame('ap'=apvet,'bp'=bpvet,'ae'=aevet,'be'=bevet,'lp'=lpvet,'le'=levet)
    }else{
      out_list[[9]]<-data.frame('ap'=ap,'bp'=bp,'ae'=ae,'be'=be,'lp'=lp,'le'=le)
    }
    cat("Done.\n")
    return(out_list)
}

## support function

Param_Estimate<-function(Aps,Bps,Aes,Bes,apold,bpold,aeold,beold){
  lp<-1/mean(as.vector(Aps))
  le<-1/mean(as.vector(Aes))
  GNL<-function(pars,data){
     alpha <- pars[[1]]
     beta <- pars[[2]]
     return (-sum(dgamma(as.vector(data), shape = alpha, rate = beta, log = TRUE)))
  }
  Fitp <- nloptr(x0 = c(apold,bpold), eval_f = GNL, lb = c(0,0), data=Bps,
                opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e5))
  Fite <- nloptr(x0 = c(aeold,beold), eval_f = GNL, lb = c(0,0), data=Bes,
                 opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e5))
  ap <- Fitp$solution[[1]]
  bp <- Fitp$solution[[2]]
  ae <- Fite$solution[[1]]
  be <- Fite$solution[[2]]
  return(list(ap,bp,ae,be,lp,le))
}
