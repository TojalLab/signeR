EstimateParameters<-function(As,Bs,aold,bold){
    l<-1/mean(as.vector(As))
    GNL<-function(pars,data){
        alpha <- pars[[1]]
        beta <- pars[[2]]
        return (-sum(dgamma(as.vector(data), shape=alpha, rate=beta, log=TRUE)))
    }
    Fitp <- nloptr(x0 = c(aold,bold), eval_f = GNL, lb = c(0,0), data=Bs,
        opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e5))
    a <- Fitp$solution[[1]]
    b <- Fitp$solution[[2]]
    return(list(a,b,l))
}
