EstimateParameters<-function(Aps,Bps,Aes,Bes,apold,bpold,aeold,beold){
    lp<-1/mean(as.vector(Aps))
    le<-1/mean(as.vector(Aes))
    GNL<-function(pars,data){
        alpha <- pars[[1]]
        beta <- pars[[2]]
        return (-sum(dgamma(as.vector(data), shape=alpha, rate=beta, log=TRUE)))
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
