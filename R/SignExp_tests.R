setGeneric("DiffExp",
    def=function(signexp_obj,labels, max_instances=200, method="kruskal.test", contrast="all",
        quant=0.5, cutoff=0.05, p.adj="BH", plot_to_file=FALSE,
        file="Diffexp_boxplot.pdf",colored=TRUE,relative=FALSE,...){
        standardGeneric("DiffExp")
    }
)
setMethod("DiffExp",signature(signexp_obj="SignExp", labels="character", 
                              max_instances="ANY", method="ANY", contrast="ANY", 
                              quant="ANY", cutoff="ANY", p.adj="ANY",
                              plot_to_file="ANY", file="ANY", colored="ANY",
                              relative="ANY"),
    function(signexp_obj, labels, method, contrast, quant, cutoff, plot_to_file,
        file, colored,relative,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        if(r > max_instances){
          select_instances<-sort(sample(1:r,max_instances,replace=FALSE))
          signexp_obj@Exp <-signexp_obj@Exp[,,select_instances]
          signexp_obj@Sign<-signexp_obj@Sign[,,select_instances]
          r<-max_instances
        }
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
            rownames(Lpval)<-signexp_obj@signames
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
            rownames(Lpval)<-signexp_obj@signames
            y.min<-min(Lpval); y.max<-max(Lpval)
            lcut<-cutoff
            Lpmed<-apply(Lpval,1,quantile,quant,na.rm=TRUE)
            cor<-rep("black",n)
            cor[Lpmed>=lcut]<-col1
            signif<-which(Lpmed>=lcut)
            my.ylab <- "AUC"
        }
        bigexp<-rep(NA,n)
        boxnames<-signexp_obj@signames
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
            names(List_sig) <- signexp_obj@signames[signif]
            names(List_pval)<- signexp_obj@signames[signif]
        }
        new_plot <-FALSE
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
            pdf(file,width=6,height=7)
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
        mysignames<-signexp_obj@signames
        if(length(signif)>0){
            mysignames[signif]<-paste(mysignames[signif],"(>",bigexp[signif],")",sep="")
        }                                   
        md<-data.frame("Sig"=rep(mysignames,each=r),
                       "Pvalues"=as.vector(t(Lpval)))
        ms = group_by(md, Sig) %>% summarise(q1=min(Pvalues),
                                             q2=quantile(Pvalues,p=0.25),
                                             q3=median(Pvalues),
                                             q4=quantile(Pvalues,p=0.75),
                                             q5=max(Pvalues))
        sig_order<-order(signexp_obj@signames)
        segments<-data.frame("Sig"=mysignames[sig_order],"x_bgn"=c(1:n)-0.4,"x_end"=c(1:n)+0.4,
                             "y_bgn"=Lpmed[sig_order],
                             "y_end"=Lpmed[sig_order],
                             "q1"=0,"q2"=0,"q3"=0,"q4"=0,"q5"=0)
        g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
            geom_boxplot(stat='identity',show.legend = FALSE) +
            geom_hline(yintercept=lcut,col=col2) +
            geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                         data = segments, show.legend = FALSE)+
            theme_classic(base_size = 15) + #theme_bw()+
            theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=0.5,face="bold")) + 
            theme(strip.background.x=element_blank()) + #7-12-22
            theme(strip.background.y=element_rect(fill='white',color='black')) + #7-12-22
            theme(axis.ticks.x=element_blank())+ #7-12-22
            labs(x="",y=my.ylab)
        #######################
        if(multicompare){
            #if(new_plot) dev.new(width=7, height=7)
            ####################### ggplot2
            Allclass<-rep(as.vector(used_labels),each=n,times=r)
            classdiffs<-rep(NA,length(Allclass))
            fullsigs<-rep(signexp_obj@signames,times=j*r)
            fullsignif<-rep(c(1:n %in% signif),times=j*r)
            for (c in 1:nclasses){
                cl<-classes[c]
                for (i in 1:length(signif)){
                    s<-signif[i]
                    thissig<-signexp_obj@signames[s]
                    classdiffs[Allclass==cl & fullsigs==thissig]<-paste(Allcomp[[i]][[c]],collapse=",")
                }
            }
            v<-as.vector(signexp_obj@Exp[signif,used,])
            const<-min(v[v>0])/100
            Logexpr<-log(signexp_obj@Exp[signif,used,] + const)
            md<-data.frame(Sig=fullsigs[fullsignif],
                           class=rep(as.vector(used_labels),each=n,times=r)[fullsignif],
                           classdiffs=classdiffs[fullsignif],
                           #exp=as.vector(Exp[signif,used,]))
                           exp=as.vector(Logexpr))
            ms = group_by(md, Sig, class ) %>% summarize(q1=min(exp),
                                                 q2=quantile(exp,p=0.25),
                                                 q3=median(exp),
                                                 q4=quantile(exp,p=0.75),
                                                 q5=max(exp),
                                                 "fc"=classdiffs[1])
            ms_ord<-c()
            for(s in signexp_obj@signames){
              for (c in 1:nclasses){
                cl<-classes[c]
                ms_ord<-c(ms_ord,which(ms$Sig==s & ms$class==cl))
              }
            }
            ms<-ms[ms_ord,]
            ymean<-min(ms$q1)
            g2<-ggplot(ms, aes(x=class,lower=q2, upper=q4, middle=q3, 
                              ymin=q1, ymax=q5)) + 
                geom_boxplot(stat='identity',show.legend = FALSE) +
                facet_wrap(vars(Sig),ncol = 3)+
                theme_classic(base_size = 8) + #theme_bw()+
                theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0,face="bold"))+
                labs(x="",y="log(Exposure)")+
                geom_text(aes(x=class,y=ymean,label=fc),data=ms)
            #######################
            final_figure <- list(g1,g2) #ggarrange(g1, g2, ncol = 1, nrow = 2)
        }else{
            final_figure<-list(g1)
        }
        #plot(figure)
        for(kk in 1:length(final_figure)){
            plot(final_figure[[kk]])
        }            
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
        colnames(mainresult)<-signexp_obj@signames
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
                opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel=1e-100, 
                    xtol_abs=1e-100, maxeval = 1e7))
            llh0<-Fit0$objective #-1 x log(lh0)
            llhs<-rep(0,nsig)
            for(s in 1:nsig){
                RSS1<-function(vet){
                    expect<-Signatures[,-s,drop=F]%*%matrix(vet,nsig-1,1) *opp
                    return( sum( expect+lgamma(found+1)-found*log(expect) ) )
                    
                }
                Fit1 <- nloptr(x0=rep(1,nsig-1), eval_f = RSS1, 
                     lb=rep(0,nsig-1), opts = list(algorithm = "NLOPT_LN_SBPLX", 
                     xtol_rel=1e-100, xtol_abs=1e-100, maxeval = 1e7))
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
                    opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel=1e-100, 
                        xtol_abs=1e-100, maxeval = 1e7))
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
                            xtol_rel=1e-100, xtol_abs=1e-100, maxeval = 1e7))
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
           def=function(Exposures,feature,method="spearman",  
                        max_instances=200,
                        cutoff_pvalue=0.05,
                        quant=0.5,
                        plot_to_file=FALSE, 
                        file="ExposureCorrelation_plot.pdf",
                        colors=TRUE,...){
               standardGeneric("ExposureCorrelation")
           }
)

setMethod("ExposureCorrelation",signature(Exposures="matrix",feature="numeric",
                                          method="ANY", max_instances="ANY",
                                          cutoff_pvalue="ANY",
                                          plot_to_file="ANY", file="ANY",
                                          colors="ANY"),
          function(Exposures, feature, method, max_instances, cutoff_pvalue,
                   plot_to_file, file, colors,...){
              de <- dim(Exposures) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              Ehat <- Exposures
              signature_names<-rownames(Ehat)
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
              boxnames<-signature_names
              boxlines<-rep(0.5,n)
              ####################### ggplot2
              md<-data.frame(Sig=signature_names,
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              sig_order<-order(signature_names)
              segments<-data.frame(Sig=signature_names[sig_order],x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpmed[sig_order],
                                   y_end=Lpmed[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_classic(base_size = 15) + #theme_bw()+
                  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              #correlation plots
              md<-data.frame(Sig=rep(signature_names,times=j),
                             exposure=as.vector(Ehat),
                             Feature=rep(feature,each=n))
              g2<-ggplot(md, aes(x=Feature,y=exposure)) + 
                  geom_point(size=1, shape=19, show.legend = FALSE) +
                  stat_smooth(method="lm", se=FALSE,col="red") +
                  facet_wrap(vars(Sig),nrow = ceiling(n/2)) +
                  theme_classic(base_size = 15) + #theme_bw()+
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
                                          method="ANY",  max_instances="ANY",
                                          cutoff_pvalue="ANY",quant="ANY",
                                          plot_to_file="ANY", file="ANY",
                                          colors="ANY"),
          function(Exposures,feature, method, max_instances,cutoff_pvalue,quant,
                   plot_to_file,file,colors,...){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              dp <- dim(Exposures@Sign) #[i,n,r]
              de <- dim(Exposures@Exp) #[n,j,r]
              i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
              if(r > max_instances){
                select_instances<-sort(sample(1:r,max_instances,replace=FALSE))
                Exposures@Exp <-Exposures@Exp[,,select_instances]
                Exposures@Sign<-Exposures@Sign[,,select_instances]
                r<-max_instances
              }
              Es <- Exposures@Exp
              Em <- Median_exp(Exposures)
              if(is.null(rownames(Em))){
                signature_names<-Exposures@signames
                rownames(Em)<-signature_names
              }else{
                signature_names<-rownames(Em)
              }
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
              boxnames<-Exposures@signames
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=rep(Exposures@signames,each=r),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              sig_order<-order(Exposures@signames)
              segments<-data.frame(Sig=Exposures@signames[sig_order],x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpmed[sig_order],
                                   y_end=Lpmed[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_classic(base_size = 15) + #theme_bw()+
                  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
              #######################
              #Correlation plots
              if(any(signif_signatures)){
                md<-data.frame(Sig=rep(Exposures@signames[signif_signatures],times=j),
                               exposure=as.vector(Em[signif_signatures,]),
                               Feature=rep(feature,each=sum(signif_signatures)))
                g2<-ggplot(md, aes(x=Feature,y=exposure)) + 
                  geom_point(size=1, shape=19, show.legend = FALSE) +
                  stat_smooth(method="lm", se=FALSE,col="red") +
                  facet_wrap(vars(Sig),nrow = ceiling(sum(signif_signatures)/2),scales="free") +
                  theme_classic(base_size = 15) + #theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="Feature",y="Exposure")
                figure <- ggarrange(g1,g2, ncol = 1, nrow = 2)
              }else{
                figure <- g1
              }
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
                        max_instances=200,
                        method="logrank", #"logrank" or "cox", #byvalue=TRUE, 
                        quant=0.5, #quantile of statistics used to attribute significance. Higher means stricter.  
                        cutoff_pvalue=0.05,
                        cutoff_hr=NA,
                        plot_to_file=FALSE, 
                        file="ExposureSurvival_plot.pdf",
                        colors=TRUE,...){
               standardGeneric("ExposureSurvival")
           }
)

setMethod("ExposureSurvival",signature(Exposures="matrix",surv="ANY",
                                       max_instances="ANY", method="ANY",
                                       quant="ANY", cutoff_pvalue="ANY", 
                                       cutoff_hr="ANY", plot_to_file="ANY", 
                                       file="ANY",colors="ANY"),
          function(Exposures, surv, max_instances, method, quant, cutoff_pvalue,
                   cutoff_hr, plot_to_file, file, colors){
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
              signature_names<-rownames(Ehat)
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
                  pdf(file,width=6,height=7)
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
              boxnames<-signature_names
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=signature_names,
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              sig_order<-order(signature_names)
              segments<-data.frame(Sig=signature_names[sig_order],x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpval[sig_order],
                                   y_end=Lpval[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_classic(base_size = 15) + #theme_bw()+
                  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=0,face="bold")) + 
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
                      pval_diff<-signif(Test$pvalue,digits=3)
                      cutround<-signif(cutvalues[m],digits=3)
                      ###############ggplot
                      maintitle=paste("Data split by exposure to Signature ",m,sep="")
                      legenlabs<-c(paste("exposure <= ",cutround,sep=""),
                                   paste("exposure > ",cutround,sep=""))
                      g2<-ggsurvplot(sf,data=data.frame(as.matrix(surv),thisgroup),
                                     ggtheme = theme_classic(base_size = 8),#theme_bw(),
                                     #font.main = c(12, "bold", "darkblue"),
                                     #font.x = c(14, "bold.italic", "red"),
                                     #font.y = c(14, "bold.italic", "darkred"),
                                     pval=pval_diff, pval.size=2,
                                     color="thisgroup",legend.labs=legenlabs)+
                          ggtitle(maintitle) +
                          labs(x="Time",y="Survival")
                      survplotlist<-c(survplotlist,list(g2$plot))
                      #####################
                  }
                  plotlist<-c(plotlist,survplotlist)
                  nextplotlist<-list()
                  final_figure <- list()
                  ppage1<-4
                  ppage2<-8
                  if(length(plotlist)>1){
                      if(length(plotlist)>(ppage1+1)){
                          thisfig<-ggarrange(plotlist=plotlist[2:(ppage1+1)], ncol=2, nrow= ppage1/2 )
                          final_figure[[1]] <- ggarrange(plotlist[[1]],thisfig,ncol=1,nrow=2)
                          nextplotlist<-plotlist[-c(1:(ppage1+1))]
                      }else{
                          thisfig<-ggarrange(plotlist=plotlist[-1], ncol=2, nrow= ppage1/2 )
                          final_figure[[1]] <- ggarrange(plotlist[[1]],thisfig,ncol=1,nrow=2)
                      }                
                  }else{
                      final_figure[[1]] <-plotlist[[1]]
                  }
                  while(length(nextplotlist)>0){
                      if(length(nextplotlist)>(ppage2)){
                          thisfig<-ggarrange(plotlist=nextplotlist[1:(ppage2)], ncol=2, nrow= ppage2/2 )
                          nextplotlist<-nextplotlist[-c(1:ppage2)]
                      }else{
                          thisfig<-ggarrange(plotlist=nextplotlist, ncol=2, nrow= ppage2/2 )
                          nextplotlist<-list()
                      }
                      final_figure<-c(final_figure,list(thisfig))     
                  }
              }else{
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
                  final_figure <- list(ggarrange(g1,figure,ncol = 1))
                  Sigfactors<-NA
              }
              for(kk in 1:length(final_figure)){
                plot(final_figure[[kk]])
              }
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Survival vs Exposure analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              #Output
              return(list("Significance"=signif,
                          "pvalues"=Pvalues_diff,
                          "pvalues_linear_model"=Pvalues_prop,
                          "limits"=cutvalues,
                          "Groups"=Sigfactors))
          }
)

#Exposure survival for SignExp object. 
setMethod("ExposureSurvival",signature(Exposures="SignExp",surv="ANY",
                                       max_instances="ANY", method="ANY",
                                       quant="ANY", cutoff_pvalue="ANY", 
                                       cutoff_hr="ANY", plot_to_file="ANY", 
                                       file="ANY",colors="ANY"),
          function(Exposures, surv, max_instances, method, quant=0.5,
                   cutoff_pvalue=0.05, cutoff_hr=NA, plot_to_file=FALSE,
                   file, colors=TRUE){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              de <- dim(Exposures@Exp) #[n,j,r]
              n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
              if(r > max_instances){
                select_instances<-sort(sample(1:r,max_instances,replace=FALSE))
                Exposures@Exp <-Exposures@Exp[,,select_instances]
                Exposures@Sign<-Exposures@Sign[,,select_instances]
                r<-max_instances
              }
              Ehat<-Median_exp(Exposures)
              Es<-Exposures@Exp
              if(is.null(rownames(Ehat))){
                  signature_names<-Exposures@signames
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
              ### Tests 
              invquant<-1-quant
              if(method=="logrank") {
                  Res.univ<-future_apply(Es,c(1,3),function(exposure){ #s
                          thisdata<-data.frame(time,os,exp=exposure)
                          const<-min(exposure[exposure>0])*1e-3
                          thisdata$lexp<-log2(exposure+const)
                          mtHL<-maxstat.test(surv~exposure, data=thisdata, smethod="LogRank",pmethod="HL", #Lr
                                             minprop=0.1, maxprop=0.9) #Lr
                          cut<-mtHL$estimate #Lr
                          group<-ifelse(exposure>cut,1,0) #Lr
                          sumup<-sum(group==1) #Lr
                          sumdown<-sum(group==0) #Lr
                          Test <- survdiff(surv~group) #Lr
                          return(c(group,#Lr
                                      round(Test$pvalue,3),#Lr
                                      cut)) #Lr
                                      #round(coxz$table[1,3],5),#Nada
                                      #as.vector(as.matrix(coxdf[1,4:10])))) #Cx
                                      #HR,Lower_CI,Upper_CI,Inv_HR,Inv_Lower_CI,Inv_Upper_CI,p
                  })
                  Sigfactors.ar<-Res.univ[c(1:j),,]#group 
                  Pvalues<-Res.univ[j+1,,] #Pvalues_diff
                  cutvalues<-Res.univ[j+2,,]
                  if(n==1){
                      Pvalues_quant<-quantile(Pvalues,probs=quant)
                  }else{
                      Pvalues_quant<-apply(Pvalues,1,function(v){quantile(v,probs=quant)})
                  }
                  #find signif signatures
                  signif_signatures <- rep(TRUE,n)
                  if(!is.na(cutoff_pvalue)){
                      signif_signatures <- signif_signatures & 
                          !is.na(Pvalues_quant) & 
                          Pvalues_quant <= cutoff_pvalue  
                  }
              }else{ #Cox model
                  Res.univ<-future_apply(Es,c(1,3),function(exposure){ #s
                          thisdata<-data.frame(time,os,exp=exposure)
                          const<-min(exposure[exposure>0])*1e-3
                          thisdata$lexp<-log2(exposure+const)
                          cph<-coxph(surv~lexp, data=thisdata) #Cx
                          thisTable<-summary(cph) #Cx
                          coxdf<-cox_as_data_frame(thisTable) #Cx
                          return(as.vector(as.matrix(coxdf[1,4:10])))
                          #HR,Lower_CI,Upper_CI,Inv_HR,Inv_Lower_CI,Inv_Upper_CI,p
                  })
                  univ.tests.ar<-Res.univ[1:7,, ]#[(j+4):(j+10),,]
                  Pvalues <-Res.univ [7,, ]#[j+10,,] #Pvalues_cox
                  HR_quantiles=apply(univ.tests.ar[1,,],1,function(v){quantile(v,probs=c(invquant,0.5,quant))})
                  cond<-HR_quantiles[2,]>=1
                  HR_quant<-ifelse(cond,HR_quantiles[1,],HR_quantiles[3,])#if HR>1, takes inferior quantile.
                  if(n==1){
                      Pvalues_quant<-quantile(Pvalues,probs=quant)
                  }else{
                      Pvalues_quant<-apply(Pvalues,1,function(v){quantile(v,probs=quant)})
                  }
                  #find signif signatures
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
                  cutvalues<-NA
              }
              #Plot prepare 
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              nsig<-sum(signif_signatures)
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=6,height=7)
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
              boxnames<-signature_names
              boxlines<-rep(0.5,n)
              ####### Plotting
              ####################### ggplot2
              md<-data.frame(Sig=rep(signature_names,each=r),
                             Pvalues=as.vector(t(Lpval)))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              sig_order<-order(signature_names)
              segments<-data.frame(Sig=signature_names[sig_order],x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lpval[sig_order],
                                   y_end=Lpval[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_classic(base_size = 15) + #theme_bw()+
                  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=0,face="bold")) + 
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
                      pval_diff<-signif(Pvalues_quant[m],digits=3)
                      this_pval<-signif(Test$pvalue,digits=3)
                      cutround<-signif(cutvals[m],digits=3)
                      ###############ggplot
                      maintitle <- signature_names[m]
                      legenlabs <- c(paste("exposure <= ",cutround,sep=""),
                                   paste("exposure > ",cutround,sep=""))
                      g2<-ggsurvplot(sf,data=data.frame(as.matrix(surv),thisgroup),
                                     ggtheme = theme_classic(base_size = 8),#theme_bw(),
                                     #font.main = c(12, "bold", "darkblue"),
                                     #font.x = c(14, "bold.italic", "red"),
                                     #font.y = c(14, "bold.italic", "darkred"),
                                     pval=pval_diff, pval.size=2,
                                     color="thisgroup",legend.labs=legenlabs)+
                          ggtitle(maintitle) +
                          labs(x="Time",y="Survival")
                      survplotlist<-c(survplotlist,list(g2$plot))
                      #####################
                  }
                  plotlist<-c(list(g1),survplotlist)
                  nextplotlist<-list()
                  final_figure <- list()
                  ppage1<-4
                  ppage2<-8
                  if(length(plotlist)>1){
                      if(length(plotlist)>(ppage1+1)){
                          thisfig<-ggarrange(plotlist=plotlist[2:(ppage1+1)], ncol=2, nrow= ppage1/2 )
                          final_figure[[1]] <- ggarrange(plotlist[[1]],thisfig,ncol=1,nrow=2)
                          nextplotlist<-plotlist[-c(1:(ppage1+1))]
                      }else{
                          thisfig<-ggarrange(plotlist=plotlist[-1], ncol=2, nrow= ppage1/2 )
                          final_figure[[1]] <- ggarrange(plotlist[[1]],thisfig,ncol=1,nrow=2)
                      }                
                  }else{
                      final_figure[[1]] <-plotlist[[1]]
                  }
                  while(length(nextplotlist)>0){
                      if(length(nextplotlist)>(ppage2)){
                          thisfig<-ggarrange(plotlist=nextplotlist[1:(ppage2)], ncol=2, nrow= ppage2/2 )
                          nextplotlist<-nextplotlist[-c(1:ppage2)]
                      }else{
                          thisfig<-ggarrange(plotlist=nextplotlist, ncol=2, nrow= ppage2/2 )
                          nextplotlist<-list()
                      }                
                      final_figure<-c(final_figure,list(thisfig))     
                  }
              }else{ #Cox model
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
                  final_figure <- list(ggarrange(g1,figure,ncol = 1))
                  Sigfactors<-NA
              }
              for(kk in 1:length(final_figure)){
                  plot(final_figure[[kk]])
              }            
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Survival vs Exposure analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              #Output
              return(list("Significance"=signif,
                          "pvalues"=Pvalues,
                          #"pvalues_linear_model"=Pvalues_prop,
                          "limits"=cutvalues,
                          "Groups"=Sigfactors))
          }
)

