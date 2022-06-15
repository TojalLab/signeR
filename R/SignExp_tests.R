setGeneric("DiffExp",
    def=function(signexp_obj,labels, method="kruskal.test", contrast="all",
        quant=0.5, cutoff=0.05, p.adj="BH", plot_to_file=FALSE,
        file="Diffexp_boxplot.pdf",colored=TRUE,relative=FALSE,...){
        standardGeneric("DiffExp")
    }
)
setMethod("DiffExp",signature(signexp_obj="SignExp", labels="character",
    method="ANY", contrast="ANY", quant="ANY", cutoff="ANY", p.adj="ANY",
    plot_to_file="ANY", file="ANY", colored="ANY",relative="ANY"),
    function(signexp_obj, labels, method, contrast, quant, cutoff, plot_to_file,
        file, colored,relative,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        cl<-labels[!is.na(labels)]
        if(all(contrast=="all")){ classes <- as.vector(levels(as.factor(cl)))
        }else{ classes <- contrast  }
        used <- labels %in% classes
        used_labels<-as.factor(as.vector(labels[used]))
        if(colored){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
        }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
        nclasses<-length(classes)
        classes<-classes[classes %in% labels]
        if(length(classes)<nclasses){
            warning("There are labels in 'contrast' not found among samples.\n",
                "Comparison will be peformed between the following groups: ",
                paste(classes,collapse=", "),".\n")
            nclasses<-length(classes)
        }
        Exposures<-signexp_obj@Exp[,used,]
        if(relative) {
            Exposures<-future_apply(Exposures,c(1,3),function(v){v/sum(v)})
            Exposures<-aperm(Exposures,c(2,1,3))
        }
        Pval_and_max<-future_apply(Exposures,c(1,3),function(vet){
            switch(method,
                   kruskal.test = {
                       Test<-kruskal.test(vet,used_labels,...)
                   }, #Kruskal test
                   Logistic.reg = {
                       Test<-Logistic.reg(vet,used_labels,...)
                   }, #regressao logistica
                   AUC = {
                       Test<-AUC(vet,used_labels,...)
                   }#AUC
            )
            mean_exp<-sapply(classes,function(cl){
                mean(vet[used_labels==cl])
            })
            maxexp <- max(mean_exp)
            return(c(Test$p.value,which(mean_exp==maxexp)[1]))
        })
        Pval<-Pval_and_max[1,,]
        #Converting Pval_and_max[2,,] to classes
        MoreExp<-future_apply(Pval_and_max[2,,],c(1,2),function(cl){
            classes[cl]
            })

        if (!method=="AUC"){
            Lpval <- -1*log(Pval) #n x r
            Lpval[Lpval>750]<-750 ### avoid Inf
            rownames(Lpval)<-paste("S",1:n,sep="")
            y.min<-min(Lpval); y.max<-max(Lpval)
            lcut <- -1*log(cutoff)
            invquant<-1-quant
            Lpmed<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
            cor<-rep("black",n)
            cor[Lpmed>=lcut]<-col1
            signif<-which(Lpmed>=lcut)
            my.ylab <- "-log(pvalue)"
        }else{
            Lpval<-Pval
            rownames(Lpval)<-paste("S",1:n,sep="")
            y.min<-min(Lpval); y.max<-max(Lpval)
            lcut<-cutoff
            Lpmed<-apply(Lpval,1,quantile,quant,na.rm=TRUE)
            cor<-rep("black",n)
            cor[Lpmed>=lcut]<-col1
            signif<-which(Lpmed>=lcut)
            my.ylab <- "AUC"
        }
        bigexp<-rep(NA,n)
        boxnames<-paste("S",1:n,sep="")
        boxlines<-rep(0.5,n)
        signif<-which(Lpmed>=lcut)
        for(k in signif){
            morevet<-MoreExp[k,]
            freqcl<-sapply(classes, function(cl){ sum(morevet==cl) })
            maxfreq <- max(freqcl)
            bigexp[k]<-classes[which(freqcl==maxfreq)[1]]
            boxnames[k]<-paste(boxnames[k],bigexp[k],sep="\n")
            boxlines[k]<-1.5
        }
        multicompare <- length(signif)>0 & nclasses>2 & method=="kruskal.test"
        if(multicompare){ #Multiple comparisons
            MCpv<-array(NA,dim=c(nclasses-1,nclasses-1,length(signif),r))
            for (k in 1:r){
                Exposure <- signexp_obj@Exp[,,'r'=k]
                if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
                for (i in 1:length(signif)){
                    s<-signif[i]
                    pkct<-PMCMRplus::kwAllPairsConoverTest(
                        x=as.vector(Exposure[s,used]), g=used_labels,
                        p.adjust.method=p.adj)
                    MCpv[,,i,k] <- -1*log(pkct$p.value)
                }
            }
            MCpvm <- apply(MCpv,c(1,2,3),quantile,invquant,na.rm=TRUE)
            Mcsig <- MCpvm > lcut
            Allcomp<-Comp_labels(Mcsig)
            names(Allcomp)<-rownames(Lpval)[signif]
            for(k in 1:length(Allcomp)){ names(Allcomp[[k]])<-classes }
            List_sig<-list()
            List_pval<-list()
            for (s in 1:length(signif)){
                Dif<-Mcsig[,,s]
                colnames(Dif)<-classes[-nclasses]
                rownames(Dif)<-classes[-1]
                List_sig<-c(List_sig,list(Dif))
                Pv<-exp(-1*MCpvm[,,s])
                colnames(Pv)<-classes[-nclasses]
                rownames(Pv)<-classes[-1]
                List_pval<-c(List_pval,list(Pv))
            }
            names(List_sig)<-paste("Signature",signif,sep="_")
            names(List_pval)<-paste("Signature",signif,sep="_")
        }
        new_plot <-FALSE
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
            pdf(file,width=7,height=7)
            par(mfrow=c(1,1),mar=c(3.1,4.2,2,2))
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=7)
                new_plot <- TRUE
            }
            par(mfrow=c(1,1),mar=c(4.2,5.2,2,2))
        }
        ####### Plotting
        ####################### ggplot2
        md<-data.frame("Sig"=rep(paste("S",1:n,sep=""),each=r),
                       "Pvalues"=as.vector(t(Lpval)))
        ms = group_by(md, Sig) %>% summarise(q1=min(Pvalues),
                                             q2=quantile(Pvalues,p=0.25),
                                             q3=median(Pvalues),
                                             q4=quantile(Pvalues,p=0.75),
                                             q5=max(Pvalues))
        segments<-data.frame("Sig"=paste("S",1:n,sep=""),"x_bgn"=c(1:n)-0.4,"x_end"=c(1:n)+0.4,
                             "y_bgn"=Lpmed,"y_end"=Lpmed,"q1"=0,"q2"=0,"q3"=0,"q4"=0,"q5"=0)
        g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
            geom_boxplot(stat='identity',show.legend = FALSE) +
            geom_hline(yintercept=lcut,col=col2) +
            geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                         data = segments, show.legend = FALSE)+
            theme_bw()+
            theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
            labs(x="",y=my.ylab)
        #######################
        if(multicompare){
            #if(new_plot) dev.new(width=7, height=7)
            ####################### ggplot2
            Allclass<-rep(as.vector(used_labels),each=n,times=r)
            classdiffs<-rep(NA,length(Allclass))
            fullsigs<-rep(paste("S",1:n,sep=""),times=j*r)
            fullsignif<-rep(c(1:n %in% signif),times=j*r)
            for (c in 1:nclasses){
                cl<-classes[c]
                for (i in 1:length(signif)){
                    s<-signif[i]
                    thissig<-paste("S",s,sep="")
                    classdiffs[Allclass==cl & fullsigs==thissig]<-paste(Allcomp[[i]][[c]],collapse=",")
                }
            }
            md<-data.frame(Sig=fullsigs[fullsignif],
                           class=rep(as.vector(used_labels),each=n,times=r)[fullsignif],
                           classdiffs=classdiffs[fullsignif],
                           #exp=as.vector(Exp[signif,used,]))
                           exp=as.vector(signexp_obj@Exp[signif,used,]))
            ms = group_by(md, Sig, class ) %>% summarize(q1=min(exp),
                                                 q2=quantile(exp,p=0.25),
                                                 q3=median(exp),
                                                 q4=quantile(exp,p=0.75),
                                                 q5=max(exp),
                                                 "fc"=classdiffs[1])
            ymean<-mean(ms$q3)
            g2<-ggplot(ms, aes(x=class,lower=q2, upper=q4, middle=q3, 
                              ymin=q1, ymax=q5)) + 
                geom_boxplot(stat='identity',show.legend = FALSE) +
                facet_wrap(vars(Sig),nrow = n)+
                theme_bw()+
                theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0,face="bold"))+
                labs(x="",y="Exposure")+
                geom_text(aes(x=class,y=ymean,label=fc),data=ms)
            #######################
            figure <- ggarrange(g1, g2, ncol = 1, nrow = 2)
        }else{
            figure<-g1
        }
        plot(figure)
        ####### Plotting end
        if(plot_to_file){
            dev.off()
            cat(paste("Differential exposure analysis results",
                "were plotted to the file",file,
                "on the current directory.",sep=" "),"\n")
        }
        Pmed<-apply(Pval,1,quantile,quant)
        signif_cond <- Pmed<=cutoff
        mainresult<-data.frame(matrix(signif_cond,1,n))
        colnames(mainresult)<-paste("S",1:n,sep="")
        result_list<-list(result=mainresult, Pvquant=Pmed, Pvalues=Pval,
            MostExposed=bigexp)
        if(multicompare){ result_list<-c(result_list,list(Differences=List_sig,
            MCPvalues=List_pval)) }
        return(result_list)
    }
)

Logistic.reg<-function(exp,labels){
    gm <- glm(labels~exp,family=binomial,data=data.frame(exp,labels))
    an <- anova(gm,test="Chisq")
    pval <- an$`Pr(>Chi)`[2]
    return(list(p.value=pval))
}

AUC<-function(exp,labels){
    roc <- roc(labels~exp,data=data.frame(exp,labels))
    auc <- roc$auc
    return(list(p.value=auc))
}


#######################################
#Signature significance test
#######################################
setGeneric("SignLRT",
    def=function(Signatures,Counts,Opp){
        standardGeneric("SignLRT")
    }
)
#Signature significance test for matrices
setMethod("SignLRT",signature(Signatures="matrix",Counts="matrix",Opp="matrix"),
    function(Signatures,Counts,Opp){
        nsig<-ncol(Signatures)
        nfeat<-nrow(Signatures)
        if(ncol(Counts)==nfeat){ Counts<-t(Counts) }
        if(ncol(Opp)==nfeat){ Opp<-t(Opp) }
        nsamp<-ncol(Counts)
        pvalues<-matrix(0,nsig,nsamp)
        LRT<-matrix(0,nsig,nsamp)
        for(k in 1:nsamp){
            found<-Counts[,k]
            opp<-Opp[,k]
            RSS0<-function(vet){
                expect<-Signatures%*%matrix(vet,nsig,1) * opp
                return( sum( expect+lgamma(found+1)-found*log(expect) ) )
                
            }
            Fit0 <- nloptr(x0=rep(1,nsig), eval_f = RSS0, lb=rep(0,nsig), 
                opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel=1e-200, 
                    xtol_abs=1e-200, maxeval = 1e10))
            llh0<-Fit0$objective #-1 x log(lh0)
            llhs<-rep(0,nsig)
            for(s in 1:nsig){
                RSS1<-function(vet){
                    expect<-Signatures[,-s,drop=F]%*%matrix(vet,nsig-1,1) *opp
                    return( sum( expect+lgamma(found+1)-found*log(expect) ) )
                    
                }
                Fit1 <- nloptr(x0=rep(1,nsig-1), eval_f = RSS1, 
                     lb=rep(0,nsig-1), opts = list(algorithm = "NLOPT_LN_SBPLX", 
                     xtol_rel=1e-200, xtol_abs=1e-200, maxeval = 1e10))
                llhs[s]<-Fit1$objective
            }
            llratios<-llh0-llhs # log(lhs/lh0)=log(lhs)-log(lh0)=llh0-llhs
            LRT[,k]<-llratios
            pvalues[,k]<-pchisq(-2*llratios,df=1,lower.tail=F) #test upper tail
        }
        return(list(pvalues,LRT))
    }
)
#Signature significance test for SignExp objects 
setMethod("SignLRT",signature(Signatures="SignExp",
    Counts="matrix",Opp="matrix"),
    function(Signatures,Counts,Opp){
        Ps <- Signatures@Sign
        dm <- dim(Ps)
        nfeat <- as.numeric(dm[1])
        nsig <- as.numeric(dm[2])
        nreal <- as.numeric(dm[3])
        if(ncol(Counts)==nfeat){ Counts<-t(Counts) }
        if(ncol(Opp)==nfeat){ Opp<-t(Opp) }
        nsamp <- ncol(Counts)
        LRT <- array(0,dim=c(n=nsig,j=nsamp,r=nreal))
        Pvalues <- array(0,dim=c(n=nsig,j=nsamp,r=nreal))
        for (r in 1:nreal){
            P <- Ps[,,r,drop=T]
            for(k in 1:nsamp){
                found<-Counts[,k]
                opp<-Opp[,k]
                RSS0<-function(vet){
                    expect<-P%*%matrix(vet,nsig,1) * opp
                    return( sum( expect+lgamma(found+1)-found*log(expect) ))
                    
                }
                Fit0 <- nloptr(x0=rep(1,nsig), eval_f = RSS0, lb=rep(0,nsig), 
                    opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel=1e-200, 
                        xtol_abs=1e-200, maxeval = 1e10))
                llh0<-Fit0$objective #-1 x log(lh0)
                llhs<-rep(0,nsig)
                for(s in 1:nsig){
                    RSS1<-function(vet){
                        expect<-P[,-s,drop=FALSE]%*%matrix(vet,nsig-1,1) * opp
                        return( sum( expect+lgamma(found+1)-found*log(expect) ))
                        
                    }
                    Fit1 <- nloptr(x0=rep(1,nsig-1), eval_f = RSS1, 
                        lb=rep(0,nsig-1),
                        opts = list(algorithm = "NLOPT_LN_SBPLX",  
                            xtol_rel=1e-200, xtol_abs=1e-200, maxeval = 1e10))
                    llhs[s]<-Fit1$objective
                }
                llratios<-llh0-llhs # log(lhs/lh0)=log(lhs)-log(lh0)=llh0-llhs
                LRT[,k,r]<-llratios
                Pvalues[,k,r]<-pchisq(-2*llratios,df=1,lower.tail=F) #upper tail
            }
        }
        Median_pval<-apply(Pvalues,c(1,2),median)
        return(list(Median_pval,Pvalues,LRT))
    }
)
#######################################
#Exposure correlation analysis
#######################################
setGeneric("ExposureCorrelation",
           def=function(Exposures,feature,method="spearman",cutoff_pvalue=0.05,quant=0.5,
                        plot_to_file=FALSE, 
                        file=NA_character_,
                        colors=NA_character_,...){
               standardGeneric("ExposureCorrelation")
           }
)

setMethod("ExposureCorrelation",signature(Exposures="matrix",feature="numeric",
                                          method="ANY",cutoff_pvalue="ANY",
                                          plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,feature,method="spearman",cutoff_pvalue=0.05,
                   plot_to_file=FALSE,file="ExposureCorrelation_plot.pdf",colors=TRUE){
              de <- dim(Exposures) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              Ehat <- Exposures
              if(is.null(rownames(Ehat))){
                  signature_names<-paste("Sig",1:n,sep="")
                  rownames(Ehat)<-signature_names
              }else{
                  signature_names<-rownames(Ehat)
              }
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(ceiling((n+1)/2),2),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(ceiling((n+1)/2),2),mar=c(4.2,5.2,2,2))
              }
              Correlations<-rep(0,n)
              Pvalues<-rep(0,n)
              for(m in 1:n){
                  sig<-Ehat[m,]
                  ct<-cor.test(sig,feature, method=method)
                  Correlations[m]<-ct$estimate
                  Pvalues[m]<-ct$p.value
              }
              ####### Plotting
              #p-values boxplots
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues) & 
                      Pvalues <= cutoff_pvalue  
              }
              Lpval <- -1*log(Pvalues)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              lcut <- -1*log(cutoff_pvalue)
              invquant<-1-quant
              Lpmed<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####################### ggplot2
              md<-data.frame(Sig=paste("S",1:n,sep=""),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              segments<-data.frame(Sig=paste("S",1:n,sep=""),x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpmed,y_end=Lpmed,q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              #correlation plots
              md<-data.frame(Sig=rep(paste("S",1:n,sep=""),times=j),
                             exposure=as.vector(Ehat),
                             Feature=rep(feature,each=n))
              g2<-ggplot(md, aes(x=Feature,y=exposure)) + 
                  geom_point(size=1, shape=19, show.legend = FALSE) +
                  stat_smooth(method="lm", se=FALSE,col="red") +
                  facet_wrap(vars(Sig),nrow = ceiling(n/2)) +
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="Feature",y="Exposure")
              figure <- ggarrange(g1,g2, ncol = 1, nrow = 2)
              figure
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Exposure correlation analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Correlations"=Correlations,
                          "Pvalues"=Pvalues))
          }
)

setMethod("ExposureCorrelation",signature(Exposures="SignExp",feature="numeric",
                                          method="ANY",cutoff_pvalue="ANY",quant="ANY",
                                          plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,feature,method="spearman",cutoff_pvalue=0.05,quant=0.5,
                   plot_to_file=FALSE,file="ExposureCorrelation_plot.pdf",colors=TRUE){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              dp <- dim(Exposures@Sign) #[i,n,r]
              de <- dim(Exposures@Exp) #[n,j,r]
              i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
              Es <- Exposures@Exp
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              Estimates_and_Pvalues<-future_apply(Es,c(1,3),function(thisexposure){
                  ct<-cor.test(thisexposure,feature, method=method)
                  return(c(ct$estimate,ct$p.value))
                  })
              Correlations<-Estimates_and_Pvalues[1,,]
              Pvalues<-Estimates_and_Pvalues[2,,]
              invquant<-1-quant
              Corr_quantiles<-apply(Correlations,1,function(v){quantile(v,probs=c(invquant,0.5,quant), na.rm=T)})
              Corr_quant<-ifelse(Corr_quantiles[2,]>=0,Corr_quantiles[1,],Corr_quantiles[3,])
              Pvalues_quant<-apply(Pvalues,1,function(v){quantile(v,probs=quant, na.rm=T)})
              ####### Plotting
              #p-values boxplots
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues_quant) & 
                      Pvalues_quant <= cutoff_pvalue  
              }
              Lpval <- -1*log(Pvalues)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              lcut <- -1*log(cutoff_pvalue)
              invquant<-1-quant
              Lpmed<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=rep(paste("S",1:n,sep=""),each=r),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              segments<-data.frame(Sig=paste("S",1:n,sep=""),x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpmed,y_end=Lpmed,q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              #Correlation plots
              Em<-Median_exp(Exposures)
              md<-data.frame(Sig=rep(paste("S",1:n,sep=""),times=j),
                             exposure=as.vector(Em),
                             Feature=rep(feature,each=n))
              g2<-ggplot(md, aes(x=Feature,y=exposure)) + 
                  geom_point(size=1, shape=19, show.legend = FALSE) +
                  stat_smooth(method="lm", se=FALSE,col="red") +
                  facet_wrap(vars(Sig),nrow = ceiling(n/2)) +
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="Feature",y="Exposure")
              figure <- ggarrange(g1,g2, ncol = 1, nrow = 2)
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(ceiling((n+1)/2),2),mar=c(3.1,4.2,2,2))
                  plot(figure)
                  dev.off()
                  cat(paste("Exposure correlation analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(ceiling((n+1)/2),2),mar=c(4.2,5.2,2,2))
                  plot(figure)
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Correlation_quantiles"=Corr_quant,
                          "Pvalues_quantiles"=Pvalues_quant,
                          "Correlations"=Correlations,
                          "Pvalues"=Pvalues))
          }
)

#############################################
#Exposure survival analysis
#############################################
setGeneric("ExposureSurvival",
           def=function(Exposures=NA, 
                        surv, 
                        method=NA_character_, #"logrank" or "cox"
                        #byvalue=TRUE, 
                        quant=0.5, #quantile of statistics used to attribute significance. Higher means stricter.  
                        cutoff_pvalue=0.05,
                        cutoff_hr=NA,
                        plot_to_file=FALSE, 
                        file=NA_character_,
                        colors=NA_character_,...){
               standardGeneric("ExposureSurvival")
           }
)

setMethod("ExposureSurvival",signature(Exposures="matrix",surv="ANY",
                                       method="ANY",#byvalue="ANY",
                                       quant="ANY", 
                                       cutoff_pvalue="ANY", cutoff_hr="ANY",
                                       plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,surv,method="logrank",quant=0.5,cutoff_pvalue=0.05, cutoff_hr=NA,
                   plot_to_file=FALSE,file="ExposureSurvival_plot.pdf",colors=TRUE){
              if(is.Surv(surv)){
                  time <- surv[,1]
                  os <- surv[,2]
              }else{
                  if( is.matrix(surv) & all(c("time","status") %in% colnames(surv)) ){
                      time <- surv[,"time"] 
                      os <- surv[,"status"]
                  }else stop("'surv' should be a Surv object or a matrix ",
                             "with 'time' and 'status' columns.\n")
                  surv<-Surv(time,os)
              }
              de <- dim(Exposures) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              Ehat<-Exposures
              if(is.null(rownames(Ehat))){
                  signature_names<-paste("Sig",1:n,sep="")
                  rownames(Ehat)<-signature_names
              }else{
                  signature_names<-rownames(Ehat)
              }
              cutvalues<-rep(NA,n)
              Sigfactors<-matrix(NA,j,n)
              Pvalues_diff<-rep(NA,n)
              Pvalues_prop<-rep(NA,n)
              Pvalues_cox<-rep(NA,n)
              HR_cox<-rep(NA,n)
              univ.tests = data.frame(HR=rep(0,n), Lower_CI=rep(0,n), Upper_CI=rep(0,n),       
                                      Inv_HR=rep(0,n), Inv_Lower_CI=rep(0,n), Inv_Upper_CI=rep(0,n),
                                      P.value=rep(1,n))
              rownames(univ.tests)<-signature_names
              univ.list<-list()
              for(m in 1:n){ # for each signature
                  exposure<-Ehat[m,]
                  thisdata<-data.frame(time,os,exp=exposure)
                  mtHL<-maxstat.test(surv~exposure, data=thisdata, smethod="LogRank",pmethod="HL",
                                     minprop=0.05, maxprop=0.95)
                  cut<-mtHL$estimate
                  group<-ifelse(exposure>cut,1,0)
                  sumup<-sum(group==1)
                  sumdown<-sum(group==0)
                  #plot(mtHL)
                  Sigfactors[,m]<-group
                  Test <- survdiff(surv~group)
                  sf<-survfit(surv~group)
                  Pvalues_diff[m]<-round(1-pchisq(Test$chisq,1),3)
                  cutvalues[m]<-cut
                  const<-min(exposure[exposure>0])*1e-3
                  thisdata$exp<-log2(exposure+const)
                  cph<-coxph(Surv(time,os)~exp, data=thisdata)
                  coxz<-cox.zph(cph)
                  Pvalues_prop[m]<-round(coxz$table[1,3],5)
                  Pvalues_cox[m]<-round(summary(cph)$coefficients[5],5)
                  thisTable<-summary(cph)
                  HR_cox[m]<-round(thisTable$coefficients[2],5)
                  univ.tests[m,]<-as.vector(cox_as_data_frame(thisTable)[1,4:10])
                  univ.list<-c(univ.list,list(cph))
              } # end each signature
              if(method=="logrank") {
                  Pvalues <- Pvalues_diff
              }else{
                  Pvalues <- Pvalues_cox
              }
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues) & 
                      Pvalues <= cutoff_pvalue  
              }
              if(!is.na(cutoff_hr)){
                  signif_signatures <- signif_signatures & 
                      !is.na(HR_cox) & 
                      abs(log(HR_cox)) >= log(cutoff_hr)  
              }
              nsig<-sum(signif_signatures)
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(ceiling((nsig+1)/2),2),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(ceiling((nsig+1)/2),2),mar=c(4.2,5.2,2,2))
              }
              Lpval <- -1*log(Pvalues)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              lcut <- -1*log(cutoff_pvalue)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=paste("S",1:n,sep=""),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              segments<-data.frame(Sig=paste("S",1:n,sep=""),x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpval,y_end=Lpval,q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              plotlist<-list(g1)
              if(method=="logrank"){
                  #KM curve plots
                  survplotlist<-list()
                  for(m in which(signif_signatures)){
                      thisgroup<-Sigfactors[,m]
                      thispvalue<-Pvalues[m]
                      #survdiff plot
                      Test <- survdiff(surv~thisgroup)
                      sf<-survfit(surv~thisgroup)
                      Tbl<-summary(sf)$table
                      surv_bottom<-as.numeric(Tbl[1,5])
                      surv_top<-as.numeric(Tbl[2,5])
                      pval_diff<-round(1-pchisq(Test$chisq,1),3)
                      cutround<-signif(cutvalues[m],digits=4)
                      ###############ggplot
                      maintitle=paste("Data split by exposure to Signature ",m,sep="")
                      legenlabs<-c(paste("exposure <= ",cutround,sep=""),
                                   paste("exposure > ",cutround,sep=""))
                      g2<-ggsurvplot(sf,data=data.frame(as.matrix(surv),thisgroup),
                                     pval=pval_diff,color="thisgroup",legend.labs=legenlabs)+
                          ggtitle(maintitle) +
                          labs(x="Time",y="Survival")
                      survplotlist<-c(survplotlist,list(g2))
                      #####################
                  }
                  plotlist<-c(plotlist,survplotlist)
              }else{
                  forestplotlist<-list()
                  for(m in which(signif_signatures)){
                      fp<-ggforest(univ.list[[m]],main = paste("Signature",m,sep=" "))
                      forestplotlist<-c(forestplotlist,list(fp))
                  }
                  plotlist<-c(plotlist,forestplotlist)
              }
              figure <- ggarrange(plotlist = plotlist,
                                  ncol = 1, nrow = 1+sum(signif_signatures))
              plot(figure)
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Survival vs Exposure analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "pvalues"=Pvalues_diff,
                          "pvalues_linear_model"=Pvalues_prop,
                          "limits"=cutvalues,
                          "Groups"=Sigfactors))
          }
)

#Exposure survival for SignExp object. 
setMethod("ExposureSurvival",signature(Exposures="SignExp",surv="ANY",
                                       method="ANY", #byvalue="ANY",
                                       quant="ANY", 
                                       cutoff_pvalue="ANY", cutoff_hr="ANY",
                                       plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,surv,method="logrank",quant=0.5,cutoff_pvalue=0.05, cutoff_hr=NA,
                   plot_to_file=FALSE,file="ExposureSurvival_plot.pdf",colors=TRUE){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              de <- dim(Exposures@Exp) #[n,j,r]
              n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
              Ehat<-Median_exp(Exposures)
              Es<-Exposures@Exp
              if(is.null(rownames(Ehat))){
                  signature_names<-paste("Sig",1:n,sep="")
                  rownames(Ehat)<-signature_names
              }else{
                  signature_names<-rownames(Ehat)
              }
              if(is.Surv(surv)){
                  time <- surv[,1]
                  os <- surv[,2]
              }else{
                  if( is.matrix(surv) & all(c("time","status") %in% colnames(surv)) ){
                      time <- surv[,"time"] 
                      os <- surv[,"status"]
                      surv<-Surv(time,os)
                  }else{ 
                      stop("'surv' should be a Surv object or a matrix ",
                             "with 'time' and 'status' columns.\n")
                  }
              }
              Res.univ<-future_apply(Es,c(1,3),function(exposure){ #s
                      thisdata<-data.frame(time,os,exp=exposure)
                      const<-min(exposure[exposure>0])*1e-3
                      thisdata$lexp<-log2(exposure+const)
                      mtHL<-maxstat.test(surv~exposure, data=thisdata, smethod="LogRank",pmethod="HL",
                                         minprop=0.1, maxprop=0.9)
                      cut<-mtHL$estimate
                      group<-ifelse(exposure>cut,1,0)
                      sumup<-sum(group==1)
                      sumdown<-sum(group==0)
                      Test <- survdiff(surv~group)
                      sf<-survfit(surv~group)
                      cph<-coxph(surv~lexp, data=thisdata)
                      coxz<-cox.zph(cph)
                      thisTable<-summary(cph)
                      coxdf<-cox_as_data_frame(thisTable)
                      return(c(group,
                                  round(1-pchisq(Test$chisq,1),3),
                                  cut,
                                  round(coxz$table[1,3],5),
                                  as.vector(as.matrix(coxdf[1,4:10])))) 
                                  #HR,Lower_CI,Upper_CI,Inv_HR,Inv_Lower_CI,Inv_Upper_CI,p
              })
              Sigfactors.ar<-Res.univ[c(1:j),,]#group 
              Pvalues_diff<-Res.univ[j+1,,]
              cutvalues<-Res.univ[j+2,,]
              Pvalues_prop<-Res.univ[j+3,,]
              Pvalues_cox<-Res.univ[j+10,,]
              univ.tests.ar<-Res.univ[(j+4):(j+10),,]
              if(method=="logrank") {
                  Pvalues <- Pvalues_diff
              }else{
                  Pvalues <- Pvalues_cox
              }
              invquant<-1-quant
              if(n==1){
                  Pvalues_quant<-quantile(Pvalues,probs=quant)
              }else{
                  Pvalues_quant<-apply(Pvalues,1,function(v){quantile(v,probs=quant)})
              }
              HR_quantiles=apply(univ.tests.ar[1,,],1,function(v){quantile(v,probs=c(invquant,0.5,quant))})
              cond<-HR_quantiles[2,]>=1
              HR_quant<-ifelse(cond,HR_quantiles[1,],HR_quantiles[3,])#if HR>1, takes inferior quantile.
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues_quant) & 
                      Pvalues_quant <= cutoff_pvalue  
              }
              if(!is.na(cutoff_hr)){
                  signif_signatures <- signif_signatures & 
                      !is.na(HR_quant) & 
                      abs(log(HR_quant)) >= log(cutoff_hr)  
              }
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              nsig<-sum(signif_signatures)
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(ceiling((nsig+1)/2),2),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(ceiling((nsig+1)/2),2),mar=c(4.2,5.2,2,2))
              }
              Lpval <- -1*log(Pvalues)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              lcut <- -1*log(cutoff_pvalue)
              Lpmed<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=rep(paste("S",1:n,sep=""),each=r),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              segments<-data.frame(Sig=paste("S",1:n,sep=""),x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpmed,y_end=Lpmed,q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              if(method=="logrank") {
                  #KM curve plots
                  survplotlist<-list()
                  Sigfactors<-apply(Sigfactors.ar,c(1,2),mean,na.rm=T)
                  #Freqs of each sample in top-exposure group 
                  #for each signature, j x n. 
                  Sigfactors[Sigfactors>=quant]<-1
                  Sigfactors[Sigfactors<=invquant]<-0
                  cutvals<-apply(cutvalues,1,median)
                  for(m in which(signif_signatures)){
                      thisgroup<-Sigfactors[,m]
                      cond <- thisgroup == 0 | thisgroup == 1 
                      thisgroup<-thisgroup[cond]
                      this.Surv<-as.matrix(surv)[cond,]
                      #survdiff plot
                      Test <- survdiff(Surv(time,status)~thisgroup,data=data.frame(this.Surv,thisgroup))
                      sf <- survfit(Surv(time,status)~thisgroup,data=data.frame(this.Surv,thisgroup))
                      Tbl<-summary(sf)$table
                      surv_bottom<-as.numeric(Tbl[1,5])
                      surv_top<-as.numeric(Tbl[2,5])
                      pval_diff<-round(1-pchisq(Test$chisq,1),3)
                      cutround<-signif(cutvals[m],digits=4)
                      ###############ggplot
                      maintitle=paste("Cohorts by exposure to Sig ",m,sep="")
                      legenlabs<-c(paste("exposure <= ",cutround,sep=""),
                                   paste("exposure > ",cutround,sep=""))
                      g2<-ggsurvplot(sf,data=data.frame(as.matrix(surv),thisgroup),
                                     pval=pval_diff,color="thisgroup",legend.labs=legenlabs)+
                          ggtitle(maintitle) +
                          labs(x="Time",y="Survival")
                      survplotlist<-c(survplotlist,list(g2$plot))
                      #####################
                  }
                  nrows<-ceiling( (1+sum(signif_signatures))/2 )
                  final_figure <- ggarrange(plotlist=c(list(g1),survplotlist),ncol=2,nrow=nrows)
              }else{
                  HR_quantiles=apply(univ.tests.ar[1,,],1,function(v){quantile(v,probs=c(invquant,0.5,quant))})
                  cond<-HR_quantiles[2,]>=1
                  HR_quant<-ifelse(cond,HR_quantiles[1,],HR_quantiles[3,])#if HR>1, takes inferior quantile.
                  Lower_CI_quant=apply(univ.tests.ar[2,,],1,function(v){quantile(v,probs=invquant)})
                  Upper_CI_quant=apply(univ.tests.ar[3,,],1,function(v){quantile(v,probs=quant)})
                  Inv_HR_quantiles=apply(univ.tests.ar[4,,],1,function(v){quantile(v,probs=c(invquant,quant))})
                  Inv_HR_quant=ifelse(cond,Inv_HR_quantiles[2,],Inv_HR_quantiles[1,])
                  Inv_Lower_CI_quant=apply(univ.tests.ar[5,,],1,function(v){quantile(v,invquant)})
                  Inv_Upper_CI_quant=apply(univ.tests.ar[6,,],1,function(v){quantile(v,probs=quant)})
                  P.value_quant=apply(univ.tests.ar[7,,],1,function(v){quantile(v,probs=quant)})
                  univ.tests = data.frame(HR=HR_quant, 
                                          Lower_CI=Lower_CI_quant, 
                                          Upper_CI=Upper_CI_quant,
                                          Inv_HR=Inv_HR_quant,
                                          Inv_Lower_CI=Inv_Lower_CI_quant,
                                          Inv_Upper_CI=Inv_Upper_CI_quant,
                                          P.value=P.value_quant
                  )
                  rownames(univ.tests)<-signature_names
                  univ.tests$labels<-signature_names
                  univ.tests$colour <- rep(c("white", "gray95"), ceiling(nrow(univ.tests)/2))[1:nrow(univ.tests)]
                  fp <- ggplot(univ.tests, aes(x = HR, y = labels, xmin = Lower_CI, xmax = Upper_CI)) +
                      geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
                      geom_pointrange(shape = 22, fill = "black") +
                      geom_vline(xintercept = 1, linetype = 3) +
                      ylab("") +
                      xlab("Hazard Ratio with 95% CI") +
                      theme_classic() +
                      scale_colour_identity() +
                      scale_y_discrete(limits = rev(univ.tests$labels)) +
                      scale_x_log10(limits = c(0.25, 4), 
                                    breaks = c(0.25, 0.5, 1, 2, 4), 
                                    labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
                      theme(axis.text.y = element_blank(), axis.title.y = element_blank())

                  univ.table<-data.frame(labels=univ.tests[,"labels"],
                                         HR_CI=paste(round(univ.tests$HR,3)," (",
                                                     round(univ.tests$Lower_CI,3),"-",
                                                     round(univ.tests$Upper_CI,3),")",sep=""),
                                         P.value=signif(univ.tests$P.value,3),
                                         colour=univ.tests$colour)
                  
                  ggtable <- ggplot(data = univ.table, aes(y = labels)) +
                      geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
                      geom_text(aes(x = 0, label = labels), hjust = 0) +
                      geom_text(aes(x = 3, label = HR_CI)) +
                      geom_text(aes(x = 3, y=n+0.5, label = "HR(CI)")) +
                      geom_text(aes(x = 7, label = P.value), hjust = 1) +
                      geom_text(aes(x = 7, y=n+0.5, label = "P.value"), hjust = 1) +
                      scale_colour_identity() +
                      theme_void() + 
                      theme(plot.margin = margin(5, 0, 35, 0))
                  
                  figure <- ggarrange(plotlist = list(ggtable,fp),
                                      ncol = 2, nrow = 1)
                  final_figure <- ggarrange(g1,figure,ncol = 1)
                  Sigfactors<-NA
              }
              plot(final_figure)
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Survival vs Exposure analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "pvalues"=Pvalues_diff,
                          "pvalues_linear_model"=Pvalues_prop,
                          "limits"=cutvalues,
                          "Groups"=Sigfactors))
          }
)
