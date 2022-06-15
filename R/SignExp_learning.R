################################################################################
# Classification:
################################################################################
setGeneric("ExposureClassify",
           def=function(signexp_obj=NA, labels, addata=NA, method="knn", k=3, weights=NA,
                        plot_to_file=FALSE, file="Classification_barplot.pdf",
                        colors=NA_character_, min_agree=0.75,...){
               standardGeneric("ExposureClassify")
           }
)
setMethod("ExposureClassify",signature(signexp_obj="ANY",labels="character",
                                       addata="ANY",
                                       method="ANY", k="ANY", weights="ANY", plot_to_file="ANY",
                                       file="ANY", colors="ANY", min_agree="ANY"),
          function(signexp_obj, labels, addata, method, k, plot_to_file, file, colors,
                   min_agree,...){
              SEinput<-FALSE
              Addinput<-FALSE
              if(class(signexp_obj)=="SignExp"){
                  if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
                  dp <- dim(signexp_obj@Sign) #[i,n,r]
                  de <- dim(signexp_obj@Exp) #[n,j,r]
                  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
                  SEinput<-TRUE
                  Exposures<-signexp_obj@Exp
              }else if(class(signexp_obj) %in% c("matrix","data.frame")){
                      n<-nrow(signexp_obj)
                      j<-ncol(signexp_obj)
                      r<-1
                      Exposures<-array(as.vector(signexp_obj),dim=c('n'=n,'j'=j,'r'=1))
              }else  if(class(signexp_obj)=="logical" & is.na(signexp_obj)){
                      r<-1
              }else{
                      stop("Provided signexp_obj should be of class SignExp or matrix.\n")
              }
              if(class(addata) %in% c("matrix","data.frame")){
                  Addinput<-TRUE
              }
              knearest<-k
              if(!length(labels)==j){
                  stop(paste("Labels should be provided for each sample.",
                             "Samples that will be classified should be",
                             "labeled with 'NA'.",sep=" "))
              }
              totest<-is.na(labels)
              totrain<-!totest
              ntest<-sum(totest)
              ntrain<-sum(totrain)
              if(ntest==0){
                  stop("All samples are already classified. There is nothing to be done.\n")
              }else{
                  if(identical(method,"knn") & ntrain<k){
                      stop(paste("The number of labeled samples is less then the number ",
                                 "of nearest neighbors used for classification, k=",k,
                                 ". Please run classification with k<=",ntrain,".",sep=""))
                  }
                  testsamples<-signexp_obj@samples[totest]
                  cl<-as.factor(labels[totrain])
                  classes<-as.vector(levels(as.factor(cl)))
                  nclass<-length(classes)
                  is.binary <- nclass<=2
                  Allprobs<-array(data = as.vector(future_apply(Exposures,3,function(D){
                      n0<-n
                      if(Addinput){ 
                          D<-rbind(D,t(addata))
                          n<-nrow(D)
                      }
                      if(all(is.na(weights))){
                          weights<-rep(1,n)
                      }
                      if (length(weights)<n){
                          weights<-rep(weights,n)[1:n]
                      }
                      if (any(weights == "standardize")){
                          cond <- which(weights == "standardize")
                          weights[cond]<-1/apply(D[cond,],1,sd)
                          weights<-as.numeric(weights)
                      }
                      D<-D*weights
                      rownames(D)<-c(paste("Sig",1:n0,sep="_"),colnames(addata))
                      D<-D/median(D[D>0]) #avoid too small values
                      Train<-t(D[,!is.na(labels),drop=FALSE])
                      Test<-t(D[,is.na(labels),drop=FALSE])
                      method<-tolower(method)
                      switch(method,
                             knn = {
                                 Classific<-kknn(cl ~ ., train=data.frame(Train,cl),
                                                 test=data.frame(Test),k=knearest,kernel="rectangular")
                                 signexp_objclass<-as.vector(Classific$fitted.values)
                                 Class_probs<-Classific$prob
                             }, #K-NearestNeighbors - pkg kknn
                             lvq = {
                                 Codebook<-lvqinit(Train,cl)
                                 LVQmod<-olvq1(Train,cl,Codebook)
                                 Classific<-kknn(cl ~ ., train=data.frame(LVQmod$x,cl=LVQmod$cl),
                                                 test=data.frame(Test),k=1,kernel="rectangular")
                                 #Classific<-lvqtest(LVQmod, Test)
                                 signexp_objclass<-as.vector(Classific$fitted.values)
                                 Class_probs<-Classific$prob
                             }, #LearningVectorQuantization - pkg class
                             logreg = {
                                 if(is.binary){
                                     lnMod <- glm(cl ~ ., data=data.frame(Train,cl), family=binomial(link="logit"))
                                     Classific<-predict(lnMod, data.frame(Test),type="response")
                                     signexp_objclass<-rep(classes[1],ntest)
                                     signexp_objclass[round(Classific)==1]<-classes[2]
                                     Class_probs<-data.frame(matrix(NA,ntest,2))
                                     colnames(Class_probs)<-classes
                                     rownames(Class_probs)<-testsamples
                                     Class_probs[,2]<-Classific
                                     Class_probs[,1]<-1-Classific
                                 }else{ 
                                     lnMod <- vglm(cl ~ ., data=data.frame(Train,cl), family=multinomial(refLevel = 1))
                                     Classific<-predict(lnMod, data.frame(Test),type="response")
                                     signexp_objclass<-apply(Classific,1,function(v){colnames(Classific)[which.max(v)]})
                                     Class_probs<-Classific
                                 }
                             }, #LogisticRegression - pkg VGAM
                             lda = {
                                 ldaMod <- lda(cl ~ ., data=data.frame(Train,cl))
                                 Classific<-predict(ldaMod, data.frame(Test))
                                 signexp_objclass<-as.vector(Classific$class)
                                 Class_probs<-Classific$posterior
                             }, #LinearDiscriminantAnalysis - pkg MASS
                             lasso = {
                                 if(is.binary){
                                     cvglm<-cv.glmnet(x=as.matrix(Train), y=cl, family="binomial")
                                     gn<-try(glmnet(x=as.matrix(Train), y=cl,
                                                    lambda=cvglm$lambda.min,family="binomial"),silent=TRUE)
                                     if(!class(gn)[1]=="try-error"){
                                         Classific<-predict(gn,as.matrix(Test),type="response")
                                         signexp_objclass<-rep(classes[1],ntest)
                                         signexp_objclass[round(Classific)==1]<-classes[2]
                                         Class_probs<-data.frame(matrix(NA,ntest,2))
                                         colnames(Class_probs)<-classes
                                         rownames(Class_probs)<-testsamples
                                         Class_probs[,2]<-Classific
                                         Class_probs[,1]<-1-Classific
                                     }else{
                                         signexp_objclass<-rep(classes[1],ntest)
                                         Class_probs<-data.frame(matrix(NA,ntest,2))
                                         colnames(Class_probs)<-classes
                                         rownames(Class_probs)<-testsamples
                                         Class_probs[,1]<-1
                                         Class_probs[,2]<-0
                                     }
                                 }else{
                                     cvglm<-cv.glmnet(x=as.matrix(Train), y=cl, family="multinomial")
                                     gn<-try(glmnet(x=as.matrix(Train), y=cl,
                                                    lambda=cvglm$lambda.min,family="multinomial"),silent=TRUE)
                                     if(!class(gn)[1]=="try-error"){
                                         Classific<-predict(gn,as.matrix(Test),type="class")
                                         signexp_objclass<-as.vector(Classific)
                                         Class_probs<-predict(gn,as.matrix(Test),type="response")
                                     }else{
                                         signexp_objclass<-rep(classes[1],ntest)
                                         Class_probs<-data.frame(matrix(NA,ntest,2))
                                         colnames(Class_probs)<-classes
                                         rownames(Class_probs)<-testsamples
                                         Class_probs[,1]<-1
                                         Class_probs[,2]<-0
                                     }
                                 }
                             }, #Lasso - pkg glmnet
                             nb = {
                                 NBmod<-naiveBayes(cl ~ ., data=data.frame(Train,cl))
                                 Classific<-predict(NBmod,data.frame(Test),type="class")
                                 signexp_objclass<-as.vector(Classific)
                                 Class_probs<-predict(NBmod,data.frame(Test),type="raw")
                             }, #NaiveBayes - pkg e1071
                             svm = {
                                 SVMmod<-svm(cl ~ ., data=data.frame(Train,cl),probability=T)
                                 Classific<-predict(SVMmod,data.frame(Test),probability=T)
                                 signexp_objclass<-as.vector(Classific)
                                 Class_probs<-attr(Classific,"probabilities")
                             }, #SupportVectorMachine - pkg e1071
                             rf = {
                                 RF<-randomForest(cl ~ ., data=data.frame(Train,cl))
                                 Classific<-predict(RF,data.frame(Test),type="response")
                                 signexp_objclass<-as.vector(Classific)
                                 Class_probs<-predict(RF,data.frame(Test),type="prob")
                             }, #RandomForest - pkg randomForest
                             ab = {
                                 if(is.binary){
                                     AdaMod <- ada(cl ~ ., data=data.frame(Train,cl), loss="logistic")
                                 }else{ 
                                     stop("Adaboost algorithm is implemented only for binary classificaton.\n")
                                 }
                                 Classific<-predict(AdaMod, data.frame(Test),type="both")
                                 signexp_objclass<-as.vector(Classific$class)
                                 Class_probs<-data.frame(Classific$probs)
                             }, #Adaboost - pkg ada
                             {
                                 Classific<-method(Train,Test,cl,...)
                                 signexp_objclass<-as.vector(Classific)
                                 Class_probs<-data.frame(matrix(0,ntest,nclass))
                                 colnames(Class_probs)<-classes
                                 rownames(Class_probs)<-testsamples
                                 Class_probs[signexp_objclass==classes[1],1]<-1
                                 Class_probs[signexp_objclass==classes[2],2]<-1
                             }#Other
                      )
                      return(as.vector(Class_probs))
                  })), dim=c('sample'=ntest,'class'=nclass,'rep'=r),
                  dimnames = list(testsamples,classes,paste("r",1:r,sep="_")))
                  FoundClasses<-future_apply(Allprobs,c(1,3),function(v){
                      resp<-rep(0,length(v))
                      resp[which.max(v)]<-1
                      return(resp)
                  })
                  Freqs<-future_apply(FoundClasses,c(1,2),sum)
                  rownames(Freqs)<-classes
                  #Summaries
                  result<-rep("",ntest)#Final classification
                  winner_freq<-rep(0,ntest)#Frequencies of classifications equal to final class
                  FinalProbs<-apply(Allprobs,c(1,2),mean)
                  for (k in 1:ntest){ #Frequencies and final classifications
                      counts<-Freqs[,k]
                      result[k] <- classes[which.max(counts)][1]
                      winner_freq[k]   <- max(counts)/sum(counts)
                  }
                  result_final<-result
                  result_final[winner_freq < min_agree] <- "undefined"
                  colnames(Freqs)<-paste(testsamples,result_final,sep="\n")
                  names(result_final)<-testsamples
                  names(winner_freq)<-testsamples
                  #Plots:
                  if(is.na(colors[1])){
                      #cols<-rainbow(nclass,start=0.5,end=0.8)
                      cols<-terrain.colors(nclass+1,0.6)[1:nclass]
                  }else{
                      if(length(colors)==nclass){
                          cols<-colors
                      }else{
                          cols<-rep(colors,nclass)[1:nclass]
                      }
                  }
                  if(plot_to_file){
                      if(length(grep("\\.pdf$",file))==0){
                          file<-paste(file,"pdf",sep=".")
                      }
                      pdf(file,width=max(5,2*ntest),height=7)
                  }else{
                      if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                                names(dev.cur()),perl=TRUE)){
                          dev.new(width=max(5,2*ntest), height=7)
                      }
                  }
                  ####### Plotting - ggplot2
                  g<-ggplot(Mymelt(Freqs), aes(x=Col,y=value,fill=Row)) + 
                      geom_col(position='stack') +  
                      theme(axis.text.x=element_text(angle=0,vjust=.5)) + 
                      coord_cartesian(expand=0) + 
                      scale_y_continuous(labels=percent_format()) + 
                      labs(x="",fill="Class")
                  plot(g)
                  if(plot_to_file){
                      dev.off()
                      outmessage<-paste("Classification results were",
                                        "plotted to the file", file,
                                        "on the current directory.",sep=" ")
                      cat(outmessage,"\n")
                  }
                  colnames(Freqs)<-paste(testsamples,result_final,sep="-")
                  return(list(class=result_final,freq=winner_freq,allfreqs=Freqs,
                              probs=FinalProbs))
              }
          }
)

Mymelt<-function(M){
    D<-data.frame(M)
    Di<-NROW(D)
    Dj<-NCOL(D)
    N<-data.frame(matrix(0.1,Di*Dj,3))
    colnames(N)<-c("Row","Col","value")
    for(i in 1:Di){
        for(j in 1:Dj){
            N[(j-1)*Di+i,1:2]<-c(rownames(D)[i],colnames(D)[j])
            N[(j-1)*Di+i,3]<-D[i,j]
        }
    }
    return(N)
}


################################################################################
# Exposure GLM multiv.:
################################################################################
setGeneric("ExposureGLM",
           def=function(Exposures,feature,cutoff_pvalue=0.05,quant=0.5,
                        plot_to_file=FALSE, 
                        file=NA_character_,
                        colors=NA_character_,...){
               standardGeneric("ExposureGLM")
           }
)

setMethod("ExposureGLM",signature(Exposures="matrix",feature="numeric",
                                  cutoff_pvalue="ANY",plot_to_file="ANY", 
                                  file="ANY",colors="ANY"),
          function(Exposures,feature,cutoff_pvalue=0.05, plot_to_file=FALSE,
                   file="ExposureGLM_plot.pdf",colors=TRUE,...){
              de <- dim(Exposures) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              Es <- data.frame(t(Exposures),target=feature)
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(1,1),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(1,1),mar=c(4.2,5.2,2,2))
              }
              gl<-glm("target~.", data=Es)
              thisTable<-summary(gl)$coefficients
              Pvalues<-thisTable[-1,4]
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
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####### Plotting
              plot(1:n,rep(lcut,n),type="n",main="",xlab="",ylab="-log(pvalue)",
                   xlim=c(0.5,n+0.5),ylim=c(y.min,y.max),xaxt="n",cex.lab=1.2)
              boxplot(data.frame(t(Lpval)),at=1:n,add=TRUE,border=cor,
                      names=rep("",n),pch=45)
              mtext(boxnames,side=1,line=boxlines,at=1:n,cex=1,las=1)
              lines(x=c(0,n+1),y=rep(lcut,2),col=col2)
              for (k in 1:n){segments(k-0.39,Lpval[k],k+0.39,Lpval[k],col=col3,lwd=2)}
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Exposure GLM analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Stats"=thisTable,
                          "Pvalues"=Pvalues))
          }
)

setMethod("ExposureGLM",signature(Exposures="SignExp",feature="numeric",
                                  cutoff_pvalue="ANY",quant="ANY",
                                  plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,feature,cutoff_pvalue=0.05,quant=0.5, plot_to_file=FALSE,
                   file="ExposureGLM_plot.pdf",colors=TRUE){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              dp <- dim(Exposures@Sign) #[i,n,r]
              de <- dim(Exposures@Exp) #[n,j,r]
              i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
              Es <- Exposures@Exp
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(1,1),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(1,1),mar=c(4.2,5.2,2,2))
              }
              Stats<-array(as.vector(
                  future_apply(Es,3,function(D){
                      thisexposures<-data.frame(t(D),target=feature)
                      gl<-glm("target~.", data=thisexposures)
                      thisTable<-summary(gl)$coefficients
                      return(thisTable)
                  })
              ),dim=c((n+1),4,r))
              Pvalues.df<-Stats[-1,4,,drop=T] #2:(n+1)
              
              Pvalues_quant<-apply(Pvalues.df,1,function(v){quantile(v,probs=quant, na.rm=T)})
              #p-values boxplots
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues_quant) & 
                      Pvalues_quant <= cutoff_pvalue  
              }
              Lpval <- -1*log(Pvalues.df)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              invquant<-1-quant
              Lp_quant<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
              lcut <- -1*log(cutoff_pvalue)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              ####### Plotting
              plot(1:n,rep(lcut,n),type="n",main="",xlab="",ylab="-log(pvalue)",
                   xlim=c(0.5,n+0.5),ylim=c(y.min,y.max),xaxt="n",cex.lab=1.2)
              boxplot(data.frame(t(Lpval)),at=1:n,add=TRUE,border=cor,
                      names=rep("",n),pch=45)
              mtext(boxnames,side=1,line=boxlines,at=1:n,cex=1,las=1)
              lines(x=c(0,n+1),y=rep(lcut,2),col=col2)
              for (k in 1:n){segments(k-0.39,Lp_quant[k],k+0.39,Lp_quant[k],col=col3,lwd=2)}
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Exposure GLM analysis results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Stats"=Stats,
                          "Pvalues"=Pvalues.df))
          }
)

################################################################################
# Coxmodel multiv.:
################################################################################
setGeneric("ExposureSurvModel",
           def=function(Exposures=NA, 
                        surv, 
                        addata=NA, 
                        quant=0.5, #quantile of statistics used to attribute significance. Higher means stricter.  
                        cutoff_pvalue=0.05,
                        cutoff_hr=NA,
                        plot_to_file=FALSE, 
                        file=NA_character_,
                        colors=NA_character_,...){
               standardGeneric("ExposureSurvModel")
           }
)

setMethod("ExposureSurvModel",signature(Exposures="matrix",surv="ANY",
                                        addata="ANY", 
                                        quant="ANY", 
                                        cutoff_pvalue="ANY", cutoff_hr="ANY",
                                        plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,surv,addata,quant=0.5,cutoff_pvalue=0.05,cutoff_hr=NA,
                   plot_to_file=FALSE,file="ExposureSurvModel_plot.pdf",colors=TRUE){
              if(is.Surv(surv)){
                  dtime <- surv[,1]
                  os <- surv[,2]
              }else{
                  if( is.matrix(surv) & all(c("time","status") %in% colnames(surv)) ){
                      dtime <- surv[,"time"] 
                      os <- surv[,"status"]
                  }else stop("'surv' should be a Surv object or a matrix ",
                             "with 'time' and 'status' columns.\n")
                  surv<-Surv(dtime,os)
              }
              de <- dim(Exposures) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              Ehat<-t(Exposures)
              if(is.null(colnames(Ehat))){
                  signature_names<-paste("Sig",1:n,sep="")
                  colnames(Ehat)<-signature_names
              }else{
                  signature_names<-colnames(Ehat)
              }
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(2,1),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(2,1),mar=c(4.2,5.2,2,2))
              }
              const<-min(Ehat[Ehat>0])*1e-3
              model0<-!is.na(addata[[1]][1])
              if(model0){ 
                  cph0<-coxph(Surv(dtime,os)~., data=data.frame(dtime,os,addata))
                  #addata<-matrix(1,j,1)
                  #colnames(addata)<-"Intercept0"
                  thisdata<-data.frame(dtime,os,log2(Ehat+const),addata)
                  colnames(thisdata)<-c("time","os",signature_names,
                                        colnames(as.data.frame(addata)))
              }else{
                  thisdata<-data.frame(dtime,os,log2(Ehat+const))
                  colnames(thisdata)<-c("time","os",signature_names)
              }
              cph<-coxph(Surv(dtime,os)~., data=thisdata)
              coxz<-cox.zph(cph)
              Pvalues_prop<-round(coxz$table[,3],5)
              thisTable<-summary(cph)
              Pvalues<-round(thisTable$coefficients[,5],5)
              HR<-round(thisTable$coefficients[,2],5)
              multiv.tests<-cox_as_data_frame(thisTable)[,4:10]
              colnames(multiv.tests)[7]<-"P.value"
              if(model0){
                  thisAnova<-anova(cph0,cph)
                  pvalAnova<-thisAnova[[4]][2]
              }
              signif_signatures <- rep(TRUE,n)
              if(!is.na(cutoff_pvalue)){
                  signif_signatures <- signif_signatures & 
                      !is.na(Pvalues) & 
                      Pvalues <= cutoff_pvalue  
              }
              if(!is.na(cutoff_hr)){
                  signif_signatures <- signif_signatures & 
                      !is.na(HR) & 
                      abs(log(HR)) >= log(cutoff_hr)  
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
              pos=c(1:n)
              box_width=0.8
              data<-multiv.tests
              np <- ifelse((data$P.value < 0.05), ifelse((data$P.value < 0.01), paste0(round(data$P.value, 3)," **"), 
                                                         paste0(round(data$P.value, 3)," *")), round(data$P.value,3))
              hr.ic <- ifelse((!is.na(data$HR)), paste0(round(data$HR, 4),' [',round(data$Lower_CI,3),',',round(data$Upper_CI,3),']'), NA)
              tabletext <- cbind(c("Variable",rownames(data)), 
                                 c("Hazard Ratio [95% CI]",hr.ic),
                                 c("P-Value",np))
              #p-values boxplots &
              #forestplots univariate
              ####### Plotting
              grid.newpage()
              pushViewport(viewport(layout = grid.layout(5,3,
                                    heights=unit(rep(1,5),c("lines","null","lines","null","lines")),
                                    widths=unit(c(4,1,3),c("lines","null","lines"))),
                                    just=c("center","center")))
              pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2, 
                                    xscale = c(0, n+1), yscale = c(0,1.5*y.max),
                                    width = n+1, height = 1.5*y.max))
              grid.rect(x = unit((n+1)/2, "native"), y = unit(0.75*y.max, "native"),
                        width = unit(n+0.5, "native"), height = unit(1.5*y.max, "native"),
                        just = "centre", hjust = NULL, vjust = NULL,
                        default.units = "native", name = NULL,
                        gp=gpar(), draw = TRUE, vp = NULL)
              grid.xaxis(at=c(1:n),label=paste("S",1:n,sep=""))
              grid.yaxis()
              grid.text("-log(pvalue)", x = unit(-2.75, "lines"),
                        gp = gpar(fontsize = 14), rot = 90)
              grid.segments(0.25, lcut, 
                            n+0.75, lcut, 
                            default.units = "native", gp = gpar(col = "red"))
              for(k in 1:n){
                  grid.segments(pos[k] - 0.5 * box_width, Lpval[k], 
                                pos[k] + 0.5 * box_width, Lpval[k], 
                                default.units = "native", gp = gpar(col = col3))
              }
              popViewport()
              pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
              forestplot(labeltext=tabletext[-2,], graph.pos=3, 
                         mean=c(NA,data$HR), 
                         lower=c(NA,data$Lower_CI), upper=c(NA,data$Upper_CI),
                         xlab="Hazard ratio [95% CI]",
                         hrzl_lines=list("2" = gpar(lwd=1, col="black")),
                         txt_gp=fpTxtGp(label=gpar(cex=1, fontface="bold"),
                                        ticks=gpar(cex=1.1),
                                        xlab=gpar(cex = 1.1),
                                        title=gpar(cex = 1.1)),
                         col=fpColors(box="black", lines="black", zero = "gray50"),
                         zero=1, cex=0.9, boxsize=0.5, colgap=unit(6,"mm"),
                         lwd.ci=2, ci.vertices=TRUE, title="",new_page = FALSE)
              
              popViewport(2)
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Exposure Survival model results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Stats"=multiv.tests))
          })

setMethod("ExposureSurvModel",signature(Exposures="SignExp",surv="ANY",
                                        addata="ANY", 
                                        quant="ANY", 
                                        cutoff_pvalue="ANY", cutoff_hr="ANY",
                                        plot_to_file="ANY", file="ANY",colors="ANY"),
          function(Exposures,surv,addata,quant=0.5,cutoff_pvalue=0.05, cutoff_hr=NA,
                   plot_to_file=FALSE,file="ExposureSurvModel_plot.pdf",colors=TRUE){
              if(!Exposures@normalized) Exposures<-Normalize(Exposures)
              if(is.Surv(surv)){
                  dtime <- surv[,1]
                  os <- surv[,2]
              }else{
                  if( is.matrix(surv) & all(c("time","status") %in% colnames(surv)) ){
                      dtime <- surv[,"time"] 
                      os <- surv[,"status"]
                  }else stop("'surv' should be a Surv object or a matrix ",
                             "with 'time' and 'status' columns.\n")
                  surv<-Surv(dtime,os)
              }
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
              model0<-!is.na(addata[[1]][1])
              if(model0){ 
                  cph0<-coxph(Surv(dtime,os)~., data=data.frame(dtime,os,addata))
              }else{
                  addata<-matrix(0,j,0)
                  # empty matrix
              }
              Allstats<-future_apply(Es,3,function(D){
                  E<-t(D)
                  const<-min(E[E>0])*1e-3
                  thisdata<-data.frame(dtime,os,log2(E+const),addata)
                  cph<-coxph(Surv(dtime,os)~., data=thisdata)
                  coxz<-cox.zph(cph)
                  pval_prop.df<-round(coxz$table[,3],5)[1:(n+1)] #<-output n+1 real
                  thisTable<-summary(cph)
                  pval.df<-round(thisTable$coefficients[,5],5)[1:n] #<-output n real
                  hr.df<-round(thisTable$coefficients[,2],5)[1:n] #<-output n real
                  mult.tests.ar<-as.vector(as.matrix(cox_as_data_frame(thisTable)[1:n,4:10])) #<-output nx7 real
                  if(model0){
                      thisAnova<-anova(cph0,cph)
                      pval.anova<-thisAnova[[4]][2] #<-output 1 real
                  }else{
                      pval.anova<-NA
                  }
                  return(c(pval.anova,pval_prop.df,pval.df,hr.df,mult.tests.ar))
              })
              pvalAnova<-as.vector(Allstats[1,])
              Pvalues_prop.df<-data.frame(Allstats[2:(n+2),])
              Pvalues.df<-data.frame(Allstats[(n+3):(2*n+2),])
              HR.df<-data.frame(Allstats[(2*n+3):(3*n+2),])
              multiv.tests.ar<-array(as.vector(Allstats[(3*n+3):(10*n+2),]),dim=c(n,7,r))
              if(model0){
                  LpvalAnova <- -1*log(pvalAnova)
                  dev.new(width=7, height=7)
                  par(mfrow=c(1,1))
                  boxplot(LpvalAnova,ylab="-log(p-value)",main="Anova.coxph comparison of surv models")
              }
              if(colors){ col1<-"darkgreen"; col2<-"red"; col3<-"blue"
              }else{ col1<-"black"; col2<-"black"; col3<-"black"  }
              if(plot_to_file){
                  if(length(grep("\\.pdf$",file))==0){file<-paste(file,"pdf",sep=".")}
                  pdf(file,width=7,height=7)
                  par(mfrow=c(n,2),mar=c(3.1,4.2,2,2))
              }else{
                  if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                            names(dev.cur()),perl=TRUE)){
                      dev.new(width=7, height=7)
                  }
                  par(mfrow=c(n,2),mar=c(4.2,5.2,2,2))
              }
              invquant<-1-quant
              Pvalues_quant<-apply(Pvalues.df,1,function(v){quantile(v,probs=quant)})
              HR_quantiles=apply(multiv.tests.ar[,1,],1,function(v){quantile(v,probs=c(invquant,0.5,quant))})
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
              Lpval <- -1*log(Pvalues.df)
              Lpval[Lpval>750]<-750 ### avoid Inf
              y.min<-min(Lpval)
              y.max<-max(Lpval)
              lcut <- -1*log(cutoff_pvalue)
              Lp_quant<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
              cor<-rep("black",n)
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-paste("S",1:n,sep="")
              boxlines<-rep(0.5,n)
              pos=c(1:n)
              box_width=0.8
              bplt<-boxplot(data.frame(t(Lpval)),at=1:n,plot=FALSE,
                            names=rep("",n),pch=45)
              boxplot_stats<-bplt$stats
              Lower_CI_quant=apply(multiv.tests.ar[,2,],1,function(v){quantile(v,probs=invquant)})
              Upper_CI_quant=apply(multiv.tests.ar[,3,],1,function(v){quantile(v,probs=quant)})
              Inv_HR_quantiles=apply(multiv.tests.ar[,4,],1,function(v){quantile(v,probs=c(invquant,quant))})
              Inv_HR_quant=ifelse(cond,Inv_HR_quantiles[2,],Inv_HR_quantiles[1,])
              Inv_Lower_CI_quant=apply(multiv.tests.ar[,5,],1,function(v){quantile(v,invquant)})
              Inv_Upper_CI_quant=apply(multiv.tests.ar[,6,],1,function(v){quantile(v,probs=quant)})
              P.value_quant=apply(multiv.tests.ar[,7,],1,function(v){quantile(v,probs=quant)})
              multiv.tests = data.frame(HR=HR_quant, 
                                        Lower_CI=Lower_CI_quant, 
                                        Upper_CI=Upper_CI_quant,
                                        Inv_HR=Inv_HR_quant,
                                        Inv_Lower_CI=Inv_Lower_CI_quant,
                                        Inv_Upper_CI=Inv_Upper_CI_quant,
                                        P.value=P.value_quant
              )
              rownames(multiv.tests)<-signature_names
              data<-multiv.tests
              np <- ifelse((data$P.value < 0.05), ifelse((data$P.value < 0.01), paste0(round(data$P.value, 3)," **"), 
                                                         paste0(round(data$P.value, 3)," *")), round(data$P.value,3))
              hr.ic <- ifelse((!is.na(data$HR)), paste0(round(data$HR, 4),' [',round(data$Lower_CI,3),',',round(data$Upper_CI,3),']'), NA)
              tabletext <- cbind(c("Variable",rownames(data)), 
                                 c("Hazard Ratio [95% CI]",hr.ic),
                                 c("P-Value",np))
              #p-values boxplots &
              #forestplots multivariate
              ####### Plotting
              grid.newpage()
              pushViewport(viewport(layout = grid.layout(5,3,
                                                         heights=unit(rep(1,5),c("lines","null","lines","null","lines")),
                                                         widths=unit(c(4,1,3),c("lines","null","lines"))),
                                    just=c("center","center")))
              pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2, 
                                    xscale = c(0, n+1), yscale = c(0,1.5*y.max),
                                    width = n+1, height = 1.5*y.max))
              grid.rect(x = unit((n+1)/2, "native"), y = unit(0.75*y.max, "native"),
                        width = unit(n+0.5, "native"), height = unit(1.5*y.max, "native"),
                        just = "centre", hjust = NULL, vjust = NULL,
                        default.units = "native", name = NULL,
                        gp=gpar(), draw = TRUE, vp = NULL)
              grid.xaxis(at=c(1:n),label=paste("S",1:n,sep=""))
              grid.yaxis()
              grid.text("-log(pvalue)", x = unit(-2.75, "lines"),
                        gp = gpar(fontsize = 14), rot = 90)
              grid.segments(0.25, lcut, 
                            n+0.75, lcut, 
                            default.units = "native", gp = gpar(col = "red"))
              for(k in 1:n){
                  grid.rect(x = pos[k], y = boxplot_stats[2, k], 
                            height = boxplot_stats[4,k] - boxplot_stats[2, k], 
                            width = 1 * box_width,
                            just = "bottom", default.units = "native", gp = gpar(col = cor[k]))
                  grid.segments(pos[k] - 0.25 * box_width, boxplot_stats[5,k], 
                                pos[k] + 0.25 * box_width, boxplot_stats[5,k], 
                                default.units = "native", gp = gpar(fill = "#CCCCCC"))
                  grid.segments(pos[k], boxplot_stats[5, k], pos[k], boxplot_stats[4,k], 
                                default.units = "native", gp = gpar(fill = "#CCCCCC"))
                  grid.segments(pos[k], boxplot_stats[1, k], pos[k], boxplot_stats[2,k], 
                                default.units = "native", gp = gpar(fill = "#CCCCCC"))
                  grid.segments(pos[k] - 0.25 * box_width, boxplot_stats[1,k], 
                                pos[k] + 0.25 * box_width, boxplot_stats[1,k], 
                                default.units = "native", gp = gpar(fill = "#CCCCCC"))
                  grid.segments(pos[k] - 0.5 * box_width, boxplot_stats[3,k], 
                                pos[k] + 0.5 * box_width, boxplot_stats[3,k], 
                                default.units = "native", gp = gpar(fill = "#CCCCCC"))
                  grid.segments(pos[k] - 0.5 * box_width, Lp_quant[k], 
                                pos[k] + 0.5 * box_width, Lp_quant[k], 
                                default.units = "native", gp = gpar(col = col3))
                  outliers<-bplt$out[bplt$group==k]
                  if(length(outliers)>0){
                      grid.points(x = rep(pos[k], length(outliers)), y = outliers, 
                                  default.units = "native", pch=45)
                  }
              }
              popViewport()
              pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
              forestplot(labeltext=tabletext, graph.pos=3, 
                         mean=c(NA,data$HR), 
                         lower=c(NA,data$Lower_CI), upper=c(NA,data$Upper_CI),
                         xlab="Hazard ratio [95% CI]",
                         hrzl_lines=list("2" = gpar(lwd=1, col="black")),
                         txt_gp=fpTxtGp(label=gpar(cex=1, fontface="bold"),
                                        ticks=gpar(cex=1.1),
                                        xlab=gpar(cex = 1.1),
                                        title=gpar(cex = 1.1)),
                         col=fpColors(box="black", lines="black", zero = "gray50"),
                         zero=1, cex=0.9, boxsize=0.5, colgap=unit(6,"mm"),
                         lwd.ci=2, ci.vertices=TRUE, title="",new_page = FALSE)
              
              popViewport(2)
              ####### Plotting end
              if(plot_to_file){
                  dev.off()
                  cat(paste("Exposure Survival model results",
                            "were plotted to the file",file,
                            "on the current directory.",sep=" "),"\n")
              }
              signif<-as.integer(signif_signatures)
              return(list("Significance"=signif,
                          "Stats"=multiv.tests))
          })



################################################################################
# Fuzzy Clustering:
################################################################################
setGeneric("FuzzyClustExp",
    def=function(signexp_obj, Med_exp=NA, Clim=NA ,method.dist="euclidean", 
        method.clust="fcm", relative=FALSE, m=2, plot_to_file=FALSE,
        file="FuzzyClustExp.pdf",colored=TRUE, iseed=NA_integer_, try_all=FALSE, 
        fast=TRUE, ...){
        standardGeneric("FuzzyClustExp")
    }
)
setMethod("FuzzyClustExp",signature(signexp_obj="SignExp", Med_exp="ANY", 
    Clim="ANY", method.dist="ANY", method.clust="ANY", 
    relative="ANY", m="ANY", plot_to_file="ANY",
    file="ANY",colored="ANY",iseed="ANY",try_all="ANY",fast="ANY"),
    function(signexp_obj, Med_exp, method.dist, method.clust,
        relative=FALSE, plot_to_file, file, colored=TRUE, iseed=NA, 
        try_all=FALSE,fast=TRUE,...){
        if(!is.na(iseed)) set.seed(seed = iseed)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        Es<-signexp_obj@Exp
        if (is.na(Clim[1])){
            Clim <- c(2,j-1)
        }else{
            if(length(Clim)==1) Clim<-rep(Clim,2)    
        } 
        Cmin<-Clim[1]
        Cmax<-Clim[2]
        if(try_all){
            step0 <- 1
        }else{
            step0 <- 2^max((floor(log2(Cmax-Cmin+1))-2),0)
        }
        if(is.na(Med_exp[1])){
            Med_exp<-as.matrix(Median_exp(signexp_obj))
        }
        if(fast){
            my_obj<-Med_exp
        }else{
            my_obj<-signexp_obj    
        }
        aseed<-iseed
        if(Cmin<Cmax){
            cat(paste("Evaluating models with the number of groups ranging from ",
                      Cmin," to ",Cmax,", please be patient.\n",sep=""))
            Ops<-Optimal_sigs(testfun=function(n){
                cat(paste("Testing ",n," groups.\n",sep=""))
                Cm<-CmeansExp(my_obj, Med_exp, C=n, method.dist, method.clust,
                              relative, iseed=aseed,...)
                U<-Cm[[1]]
                PBMF<-as.vector(future_apply(Es,3,function(E){PBMFindex(U,Data=E,m)}))
                return(list(median(PBMF),PBMF))
            },
            liminf=Cmin,limsup=Cmax,step=step0,significance=FALSE
            )
            bestn<-as.numeric(Ops[[1]])
            rm(Ops)
            cat(paste("Optimum number of groups is ",bestn,
                      ". Performing final clustering.\n",sep=""))
        }else{
            bestn<-Cmin
            cat(paste("Performing clustering in ",bestn," groups.\n",sep=""))
        }
        # cat(class(signexp_obj));cat("\n")
        # cat(class(Med_exp));cat("\n")
        # cat(class(bestn));cat("\n")
        # cat(bestn);cat("\n")
        # cat(method.dist);cat("\n")
        # cat(method.clust);cat("\n")
        # cat(relative);cat("\n")
        # cat(iseed);cat("\n")
        Cm<-CmeansExp(signexp_obj, Med_exp, C=bestn, method.dist, method.clust,
            relative, iseed=aseed,...)
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=max(5,2*j),height=7)
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                      names(dev.cur()),perl=TRUE)){
                dev.new(width=max(5,2*j), height=7)
            }
        }
        if(colored){
            clr = rev(colorRampPalette(brewer.pal(name="Spectral",n=11))(100))
        }else{
            clr=rev(colorRampPalette(brewer.pal(name="Greys",n=9))(100)) 
        }
        pheatmap(Cm$Meanfuzzy,border_color=NA, color=clr, 
                 clustering_method='ward.D2',
                 clustering_distance_rows='canberra',
                 cluster_cols = FALSE)
        if(plot_to_file){
            dev.off()
        }
        return(Cm)
    }
)

#CmeansExp
setGeneric("CmeansExp",
    def=function(signexp_obj, Med_exp=NA, C, method.dist="euclidean", 
        method.clust="fcm", relative=FALSE, iseed=NA_integer_, ...){
        standardGeneric("CmeansExp")
    }
)
#For matrix:
setMethod("CmeansExp",signature(signexp_obj="matrix", Med_exp="ANY", C="ANY",
                                method.dist="ANY", method.clust="ANY", 
                                relative="ANY", iseed="ANY"),
          function(signexp_obj, Med_exp, C, method.dist, method.clust,
                   relative, plot_to_file, file, colored, iseed=NA,...){
              # initialize random seed
              if(!is.na(iseed)) set.seed(seed = iseed)
              de <- dim(signexp_obj) #[n,j]
              n<-de[[1]]; j<-de[[2]]
              if(is.na(Med_exp[1])){
                  Med_exp<-signexp_obj
              }
              if(relative){ Med_exp<-t(t(Med_exp)/colSums(Med_exp)) }
              colnames(Med_exp)<-colnames(signexp_obj)
              if(method.clust=="km"){
                  baseclust<-kmeans(t(Med_exp),centers=C)
                  basefuzzy<-t(sapply(baseclust$cluster,function(n){
                      as.numeric(c(1:C)==n)
                  }))
              }else{ 
                  if(method.clust=="fcm"){
                      baseclust<-ppclust::fcm(t(Med_exp),centers=C)
                  }else if (method.clust=="pcm"){
                      baseclust<-ppclust::pcm(t(Med_exp),centers=C)
                  }else if (method.clust=="fpcm"){
                      baseclust<-ppclust::fpcm(t(Med_exp),centers=C)
                  }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
                  basefuzzy<-baseclust$u
              }
              return(list(Meanfuzzy=basefuzzy,AllFuzzy=basefuzzy,Fuzzy=basefuzzy))
          })
#For SignExp
setMethod("CmeansExp",signature(signexp_obj="SignExp", Med_exp="ANY", C="ANY",
    method.dist="ANY", method.clust="ANY", relative="ANY",iseed="ANY"),
    function(signexp_obj, Med_exp, C, method.dist, method.clust,
        relative, iseed=NA,...){
        # initialize random seed
        if(!is.na(iseed)) set.seed(seed = iseed)
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        if(is.na(Med_exp[1])){
            Med_exp<-Median_exp(signexp_obj)
        }
        if(relative){ Med_exp<-t(t(Med_exp)/colSums(Med_exp)) }
        colnames(Med_exp)<-signexp_obj@samples
        if(method.clust=="km"){
            baseclust<-kmeans(t(Med_exp),centers=C)
            basefuzzy<-t(sapply(baseclust$cluster,function(n){
                as.numeric(c(1:C)==n)
            }))
        }else{ 
            if(method.clust=="fcm"){
                baseclust<-ppclust::fcm(t(Med_exp),centers=C)
            }else if (method.clust=="pcm"){
                baseclust<-ppclust::pcm(t(Med_exp),centers=C)
            }else if (method.clust=="fpcm"){
                baseclust<-ppclust::fpcm(t(Med_exp),centers=C)
            }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
            basefuzzy<-baseclust$u
        }
        Es<-signexp_obj@Exp
        Fuzzy2<-future_apply(Es,3,function(Exposure){
            if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
            if(relative){ Exposure<-t(t(Exposure)/colSums(Exposure)) }
            colnames(Exposure)<-signexp_obj@samples
            if(method.clust=="km"){
                thisclust<-kmeans(t(Exposure),centers=C)
                thisfuzzy<-t(sapply(thisclust$cluster,function(n){
                    as.numeric(c(1:C)==n)
                }))
            }else{ 
                if(method.clust=="fcm"){
                    thisclust<-ppclust::fcm(t(Exposure),centers=C)
                }else if (method.clust=="pcm"){
                    thisclust<-ppclust::pcm(t(Exposure),centers=C)
                }else if (method.clust=="fpcm"){
                    thisclust<-ppclust::fpcm(t(Exposure),centers=C)
                }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
                thisfuzzy<-thisclust$u
            }
            #hungarian algorithm to assign clusters
            D<-sapply(1:C,function(j){
                sapply(1:C,function(i){
                    sum(abs(basefuzzy[,i]-thisfuzzy[,j]))
                })
            })#rows correspond to basefuzzy clusters, cols to thisfuzzy clusters 
            assignment <- clue::solve_LSAP(D)
            thisfuzzy <- thisfuzzy[,as.vector(assignment)]
            colnames(thisfuzzy)<-colnames(basefuzzy)
            return(thisfuzzy)
        })
        Fuzzy<-array(as.vector(Fuzzy2),dim=c('j'=j, 'c'=C, 'r'=r))
        rownames(Fuzzy)<-rownames(basefuzzy)
        colnames(Fuzzy)<-colnames(basefuzzy)
        Meanfuzzy<-apply(Fuzzy,c(1,2),mean)
        return(list(Meanfuzzy=Meanfuzzy,
                    AllFuzzy=Fuzzy,
                    Fuzzy=basefuzzy))
    })
####
PBMFindex<-function(U,Data,m=2){
    U<-as.matrix(U)
    Data<-as.matrix(Data)
    k<-NCOL(U)
    n<-NROW(U)
    Um<-U^m
    Centroids<-t(t(Data%*%Um)/colSums(Um))
    v0<-rowMeans(Data)
    D0<-Data-v0
    e<-sum(sqrt(colSums(D0^2)))
    dk<-max(dist(t(Centroids)))
    Dc<-(as.matrix(dist(t(cbind(Data,Centroids))))[1:n,(n+1):(n+k)])^2
    #Dc<-(proxy::dist(t(Data),t(Centroids)))^2
    J<-sum(Um*Dc)
    PBMF<-((e*dk)/(k*J))^2
    return(PBMF)
}

################################################################################
# Hierarquical Clustering:
################################################################################
setGeneric("HClustExp",
    def=function(signexp_obj, Med_exp=NA, method.dist="euclidean", 
        method.hclust="average", use.cor=FALSE, relative=FALSE, 
        plot_to_file=FALSE, file="HClustExp_dendrogram.pdf", colored=TRUE,...){
        standardGeneric("HClustExp")
    }
)
setMethod("HClustExp",signature(signexp_obj="SignExp", Med_exp="ANY",
    method.dist="ANY", method.hclust="ANY",
    use.cor="ANY", relative="ANY", plot_to_file="ANY",
    file="ANY",colored="ANY"),
    function(signexp_obj, Med_exp, method.dist, method.hclust,
        use.cor, relative, plot_to_file,
        file,colored,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        if(is.na(Med_exp[1])){
            Med_exp<-Median_exp(signexp_obj)
        }
        if(relative){ Med_exp<-t(t(Med_exp)/colSums(Med_exp)) }
        # Med_exp: (n,p) matrix, n-samples, p-variables
        # Use custom distance function
        if(is.function(method.dist)) {
            distance <- method.dist(Med_exp)
        } else {
            distance<-dist.pvclust(Med_exp,method=method.dist,use.cor=use.cor)
        }
        
        # ward -> ward.D
        # only if R >= 3.1.0
        if(method.hclust == "ward" && getRversion() >= '3.1.0') {
            method.hclust <- "ward.D"
        }
        Med_exp.hclust <- hclust(distance, method=method.hclust)
        # multiscale bootstrap
        r <- list(1.0)
        mboot <-list(Exp.hclust(signexp_obj=signexp_obj, 
            object.hclust=Med_exp.hclust, method.dist=method.dist, 
            use.cor=use.cor, method.hclust=method.hclust, relative=relative)) 
        
        result<-Expclust.merge(Med_exp,object.hclust=Med_exp.hclust,mboot=mboot)
        if(plot_to_file){
            pdf(file)
        }
        ####### Plotting
        plot(result)  # plot.pvclust
        ####### Plotting end
        if(plot_to_file){
            dev.off()
        }
        return(result)
    }
    
)
Exp.hclust <- function(signexp_obj, object.hclust, method.dist, use.cor, 
    method.hclust, relative=FALSE)
{ 
    if(!signexp_obj@normalized){
        signexp_obj<-Normalize(signexp_obj)
    }
    Es<-signexp_obj@Exp
    n<-dim(Es)[[1]]
    j<-dim(Es)[[2]]
    r<-dim(Es)[[3]]
    
    pattern   <- hc2split(object.hclust)$pattern
    edges.cnt <- table(factor(pattern)) - table(factor(pattern))
    
    for(k in 1:r){
        Exposure <- Es[,,'r'=k]
        if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
        if(relative){ Exposure<-t(t(Exposure)/colSums(Exposure)) }
        colnames(Exposure)<-signexp_obj@samples
        if(is.function(method.dist)) {
            suppressWarnings(distance  <- method.dist(Exposure))
        } else {
            suppressWarnings(distance  <- dist.pvclust(Exposure,
                method=method.dist,use.cor=use.cor))
        }
        if(all(is.finite(distance))) { # check if distance is valid
            x.hclust  <- hclust(distance,method=method.hclust)
            pattern.i <- hc2split(x.hclust)$pattern # split
            edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
        } else {
            x.hclust <- NULL
            warning(paste("inappropriate distance matrices",
                " are omitted in computation: r = ", k), call.=FALSE)
        }
    }
    boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
        method.hclust=method.hclust, nboot=r, size=n, r=1, store=NA)
    class(boot) <- "boot.hclust"
    return(boot)
}

Expclust.merge <- function(data, object.hclust, mboot){
    pattern <- hc2split(object.hclust)$pattern
    r     <- unlist(lapply(mboot,"[[","r"))
    nboot <- unlist(lapply(mboot,"[[","nboot"))
    store <- lapply(mboot,"[[", "store")
    rl <- length(mboot)
    ne <- length(pattern)
    edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
    row.names(edges.bp) <- pattern
    names(edges.cnt) <- paste("r", 1:rl, sep="")
    for(j in 1:rl) {
        edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
        edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
    }
    ms.fitted <- lapply(as.list(1:ne),
        function(x, edges.bp, r, nboot){
            msfit(as.vector(t(edges.bp[x,])), r, nboot)},
        edges.bp, r, nboot)
    class(ms.fitted) <- "mslist"
    p    <- lapply(ms.fitted,"[[","p")
    se   <- lapply(ms.fitted,"[[","se")
    coef <- lapply(ms.fitted,"[[","coef")
    au    <- unlist(lapply(p,"[[","au"))
    bp    <- unlist(lapply(p,"[[","bp"))
    se.au <- unlist(lapply(se,"[[","au"))
    se.bp <- unlist(lapply(se,"[[","bp"))
    v     <- unlist(lapply(coef,"[[","v"))
    cc    <- unlist(lapply(coef,"[[","c"))
    pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
    edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
        v=v, c=cc, pchi=pchi)
    row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne
    result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
        msfit=ms.fitted, nboot=nboot, r=r, store=store)
    class(result) <- "pvclust"
    return(result)
}

#Internal functions from package pvclust:

hc2split <- function(x)
{
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()
    
    for(i in 1:(n-1)){
        ai <- A[i,1]
        
        if(ai < 0)
            B[[i]] <- -ai
        else
            B[[i]] <- B[[ai]]        
        
        ai <- A[i,2]
        
        if(ai < 0)
            B[[i]] <- sort(c(B[[i]],-ai))
        else
            B[[i]] <- sort(c(B[[i]],B[[ai]]))
    }
    
    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
    
    for(i in 1:(n-1)){
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
            CC[i,bi[j]] <- 1
    }
    
    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
    
    return(split)
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
    if(!is.na(pmatch(method,"correlation"))){
        res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
        attr(res,"method") <- "correlation"
        return(res)
    }
    else if(!is.na(pmatch(method,"abscor"))){
        res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
        attr(res,"method") <- "abscor"
        return(res)
    }
    else if(!is.na(pmatch(method,"uncentered"))){
        if(sum(is.na(x)) > 0){
            x <- na.omit(x)
            warning("Rows including NAs were omitted")
        }
        x  <- as.matrix(x)
        P  <- crossprod(x)
        qq <- matrix(diag(P),ncol=ncol(P))
        Q  <- sqrt(crossprod(qq))
        res <- as.dist(1 - P/Q)
        attr(res,"method") <- "uncentered"
        return(res)
    }
    else
        dist(t(x),method)
}

