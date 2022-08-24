################################################################################
# Parameters:
#
# Input data
# M: matrix of mutation counts, each row corresponds to one sample 
#   (dimension = nsamples x 96). Can be a path to a file.
# Mheader: if M is to be read from a file, does it contains a header line? 
#          Otherwise, if M is a matrix, does it have mutations and samples as 
#          colnames and rownames?  
# samples: if "rows" M is supposed to contain one sample per row. 
#          Otherwise, one sample per column is expected.
# Opport: matrix of mutation opportunities, each column corresponds to 
#         one sample (96 x nsamples)
# Opportheader: if Opport is to be read from a file, does it contains 
#               a header line?
# P: matrix of signatures, can be an initial guess or a fixed matrix
# fixedP: if P matrix is provided, should it be used as an initial guess 
#   that will be altered by the sampler (fixedP = FALSE, default) or should 
#   it be held fixed along the algorithm (fixedP = TRUE)
#
# Numerical parameters of the model
# nsig: number of signatures, which can be provided or estimated 
#       by the algorithm
# nlim: range of possible values for the number of signatures
# try_all: if true, all possible values for n will be tested
# BICsignificance: should the number of signatures be selected by a simple 
#                  comparison between BIC medians or by a statistical criterion  
#                  of comparison? Default is FALSE 
# critical_p: critical p-value required to consider differences in BIC values
#             as significant. Default is 0.05. 
# ap, bp: hyperparameters of rate parameters used to generate signatures 
#        (gamma distributions)
# ae, be: hyperparameters of rate parameters used to generate exposures 
#        (gamma distributions)
# lp: hyperparameters of shape parameters used to generate signatures 
#     (exponential distributions)
# le: hyperparameters of shape parameters used to generate exposures 
#     (exponential distributions)
# varp.ap: variance of the gamma distribution used to generate proposals for 
#          shape parameters of signatures
# varp.ae: variance of the gamma distribution used to generate proposals for 
#          shape parameters of exposures
# Sampler parameters
# start: NMF algorithm used to generate initial values for signatures and 
#        exposures, options: "brunet","KL","lee","Frobenius","offset","nsNMF",
#                            "ls-nmf","pe-nmf","siNMF","snmf/r" or "snmf/l".
# testing_burn: number of burning iterations of the Gibbs sampler used to 
#               estimate the number of signatures in data. 
# testing_eval: number of iterations of the Gibbs sampler used to estimate 
#               the number of signatures in data. 
# main_burn: number of burning iterations of the final Gibbs sampler
# main_eval: number of iterations of the final Gibbs sampler
#
# Further options
# estimate_hyper: if TRUE, algorithm estimates optimal values of 
#                 ap,bp,ae,be,lp,le. Start values can still be provided.
# EMit_lim: limit of EM iterations for the estimation of ap,bp,ae,be,lp,le.
# EM_eval: number of samples generated at each EM iteration.
# parallelization: strategy for parallel computing. Defaults to "multisession" 
################################################################################

signeR<-function(M, Mheader=TRUE, samples = "rows", Opport=NA, 
    Oppheader=FALSE, P=NA, fixedP=FALSE, nsig=NA,nlim=c(NA,NA),
    try_all=FALSE, BICsignificance=FALSE, critical_p = 0.05,
    ap=NA,bp=NA,ae=NA,be=NA,lp=NA,le=NA, var.ap=10,var.ae=10,
    start='lee', testing_burn=1000, testing_eval=1000,
    main_burn=10000, main_eval=2000, estimate_hyper=FALSE, 
    EMit_lim=100, EM_eval=100, parallelization="multisession"){ 
    # Read data
    if(is.character(M)){ 
        Mread<-read.table(M, header=Mheader, check.names=FALSE)
    }else if (is.data.frame(M)){
        Mread<-M
    }else if (is.matrix(M)){
        Mread <-M
    }else stop("M should be a matrix, a data frame or a link to a file.")
    if(samples=="rows"){
        M<-t(as.matrix(Mread))
    }else{
        M<-as.matrix(Mread)
    }
    if(Mheader){
        samplenames<-colnames(M)
        mutnames<-rownames(M)
    }else{
        samplenames<-paste("sample",1:NCOL(M),sep="_")
        bases=c("A","C","G","T")
        mutnames<-paste(rep(bases,each=4,6),rep(c("C","T"),each=48),
            rep(bases,24),">",
            c(rep(bases[-2],each=16),rep(bases[-4],each=16)),sep="")
    }
    if(is.character(Opport)){
        Opportread <- read.table(Opport,header=Oppheader)
    }else if (is.data.frame(Opport)){
        Opportread <- Opport
    }else if (is.matrix(Opport)){
        Opportread <- Opport
    }else if(is.na(Opport)){
        Opportread<-M
        Opportread[,]<-1
    }else{ 
        stop("Opport should be a matrix, a data frame, a link to a file or NA.")
    }
    if(all(dim(Opportread)==dim(M))){
        W<-as.matrix(Opportread)
    }else if (all(dim(t(Opportread))==dim(M))){
        W<-t(as.matrix(Opportread))
    }else stop("Dimensions of M and Opport do not match.")
    #Initialize parameters
    if (is.na(P[1])){ fixedP <- FALSE }else{ nsig <- NCOL(P) }
    if (fixedP){
        if(sum(is.na(c(ae,be,le)))>0){
            eh<-TRUE
            if(is.na(ae)) ae<-1
            if(is.na(be)) be<-1
            if(is.na(le)) le<-1
        }else{
            eh<-estimate_hyper
        }
    }else{
        if(sum(is.na(c(ap,bp,ae,be,lp,le)))>0){
            eh<-TRUE
            if(is.na(ap)) ap<-1
            if(is.na(bp)) bp<-1
            if(is.na(ae)) ae<-1
            if(is.na(be)) be<-1
            if(is.na(lp)) lp<-1
            if(is.na(le)) le<-1
        }else{
            eh<-estimate_hyper
        }
    }
    if(is.na(nsig) & !fixedP ){
        #Search for best n
        if(is.na(nlim[1])){
            Nmin <- 1
        }else{
            Nmin <- nlim[1]
        }
        if(is.na(nlim[2])){
            Nmax <- min(dim(M))-1
        }else{
            Nmax <- nlim[2]
        }
        if(try_all){
            step0 <- 1
        }else{
            step0 <- 2^max((floor(log2(Nmax-Nmin+1))-2),0)
        }
        cat(paste("Evaluating models with the number of signatures",
            " ranging from ",Nmin," to ",Nmax,", please be patient.\n",sep=""))
        Ops<-Optimal_sigs(testfun=function(n){
            cat(paste("Evaluating",n,"signatures.\n",sep=" "))
            if(!eh){be<-be/(n^0.5)}
            ebNMF<-eBayesNMF(M,W,n,ap,bp,ae,be,lp,le,
                             var.ap,var.ae,P=NA,fixP=FALSE,
                start=start, burn_it=testing_burn,eval_it=testing_eval,
                estimate_hyper=eh, EM_lim=EMit_lim, EM_eval_it=EM_eval,
                keep_param=FALSE)
            bics<-ebNMF[[3]]
            HH<-ebNMF[[9]]
            rm(ebNMF)
            return(list(bics,HH))
        },
            liminf=Nmin, limsup=Nmax, step=step0,
            significance=BICsignificance, pcrit=critical_p, 
            parplan = parallelization)
        nopt<-Ops[[1]]
        evalued_n<-Ops[[2]]
        Test_BICs<-list()
        HH<-list()
        for (k in 1:length(evalued_n)){
            Test_BICs[[k]]<-Ops[[3]][[k]][[1]]
            HH[[k]]<-Ops[[3]][[k]][[2]]
        }
        cat(paste("The optimal number of signatures is ",nopt,".\n",sep=""))
        BestHH<-HH[[which(evalued_n==nopt)]]
        best_hyperparam<-BestHH[NROW(BestHH),]
        ap<-best_hyperparam$ap
        bp<-best_hyperparam$bp
        ae<-best_hyperparam$ae
        be<-best_hyperparam$be
        lp<-best_hyperparam$lp
        le<-best_hyperparam$le
        eh<-FALSE
    }else{
        if (fixedP){
            nopt<-NCOL(P)
        }else{
            nopt<-nsig
        }
        Ops<-list(NA,NA,NA,NA)
        Test_BICs<-NA
        HH<-NA
    }
    if(fixedP & !is.null(colnames(P)[1])){
        signames<-colnames(P)
    }else{
        signames<-paste("S",1:nopt,sep="")
    }
    #cat(paste("Running Gibbs sampler for ",nopt," signatures.\n",sep=""))
    Final_run<-eBayesNMF(M,W,n=nopt,ap,bp,ae,be,lp,le,
        var.ap,var.ae,P,fixP=fixedP,
        start=start,burn_it=main_burn,eval_it=main_eval, 
        estimate_hyper=eh,EM_lim=EMit_lim, EM_eval_it=EM_eval,
        keep_param=TRUE)
    SE<-Final_run[[4]]
    SE<-setSamples(SE,samplenames)
    SE<-setMutations(SE,mutnames)
    SE<-setSignames(SE,signames)
    result<-list(Nsign=nopt,
        tested_n=Ops[[2]],
        Test_BICs=Test_BICs,
        Phat=Final_run[[1]],
        Ehat=Final_run[[2]],
        SignExposures=SE,
        BICs=Final_run[[3]],
        Hyper_paths=HH)
    #cat("Done.\n")
    return(result)
}


