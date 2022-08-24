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
                      rownames(D)<-c(signexp_obj@signames,colnames(addata))
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
              if(!is.null(rownames(Exposures)[1])){
                  boxnames<-rownames(Exposures)
              }else{
                  boxnames<-paste("S",1:n,sep="")
              }
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
              boxnames<-Exposures@signames
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
              invquant<-1-quant
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
                  LpvalAnova <- -1*log(pvalAnova)
                  y.max_anova<-max(LpvalAnova)
                  y.min_anova<-min(LpvalAnova)
                  Lp_quant_Anova<-quantile(LpvalAnova,invquant,na.rm=TRUE)
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
              lcut <- -1*log(cutoff_pvalue)
              Lp_quant<-Lpval
              cor[signif_signatures]<-col1
              bigexp<-rep(NA,n)
              boxnames<-signature_names
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
              ####### Plotting p-values boxplots & forestplots univariate
              ####################### ggplot2
              md<-data.frame(Sig=signature_names,
                             Pvalues=as.vector(Lpval))
              ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                   q2=quantile(Pvalues,p=0.25),
                                                   q3=median(Pvalues),
                                                   q4=quantile(Pvalues,p=0.75),
                                                   q5=max(Pvalues))
              sig_order<-order(signature_names)
              segments<-data.frame(Sig=signature_names[sig_order],x_bgn=c(1:n)-0.4,x_end=c(1:n)+0.4,
                                   y_bgn=Lp_quant[sig_order],
                                   y_end=Lp_quant[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                geom_hline(yintercept=lcut,col=col2) +
                geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                             data = segments, show.legend = FALSE)+
                theme_bw()+
                theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                labs(x="",y="-log(pvalue)")
              
              fp <- ggplot(multiv.tests, aes(x = HR, y = labels, xmin = Lower_CI, xmax = Upper_CI)) +
                geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
                geom_pointrange(shape = 22, fill = "black") +
                geom_vline(xintercept = 1, linetype = 3) +
                ylab("") +
                xlab("Hazard Ratio with 95% CI") +
                theme_classic() +
                scale_colour_identity() +
                scale_y_discrete(limits = rev(multiv.tests$labels)) +
                scale_x_log10(limits = c(0.25, 4), 
                              breaks = c(0.25, 0.5, 1, 2, 4), 
                              labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
                theme(axis.text.y = element_blank(), axis.title.y = element_blank())
              
              multiv.table<-data.frame(labels=multiv.tests[,"labels"],
                                       HR_CI=paste(round(multiv.tests$HR,3)," (",
                                                   round(multiv.tests$Lower_CI,3),"-",
                                                   round(multiv.tests$Upper_CI,3),")",sep=""),
                                       P.value=signif(multiv.tests$P.value,3),
                                       colour=multiv.tests$colour)
              
              ggtable <- ggplot(data = multiv.table, aes(y = labels)) +
                geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
                geom_text(aes(x = 0, label = labels), hjust = 0) +
                geom_text(aes(x = 3, label = HR_CI)) +
                geom_text(aes(x = 3, y=n+0.5, label = "HR(CI)")) +
                geom_text(aes(x = 7, label = P.value), hjust = 1) +
                geom_text(aes(x = 7, y=n+0.5, label = "P.value"), hjust = 1) +
                scale_colour_identity() +
                theme_void() + 
                theme(plot.margin = margin(5, 0, 35, 0))
              
              forestplot <- ggarrange(plotlist = list(ggtable,fp),
                                      ncol = 2, nrow = 1)
              if(model0){
                md<-data.frame(Sig="Anova",
                               Pvalues=as.vector(LpvalAnova))
                ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                     q2=quantile(Pvalues,p=0.25),
                                                     q3=median(Pvalues),
                                                     q4=quantile(Pvalues,p=0.75),
                                                     q5=max(Pvalues))

                segments<-data.frame(Sig="Anova",x_bgn=0.6,x_end=1.4,
                                     y_bgn=Lp_quant_Anova,y_end=Lp_quant_Anova,
                                     q1=Lp_quant_Anova,q2=Lp_quant_Anova,
                                     q3=Lp_quant_Anova,q4=Lp_quant_Anova,
                                     q5=Lp_quant_Anova)
                g2<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
                final_figure <- ggarrange(g1,forestplot,g2,ncol = 1)
              }else{
                final_figure <- ggarrange(g1,forestplot,ncol = 1)
              } 
              plot(final_figure)
              
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
              invquant<-1-quant
              Ehat<-Median_exp(Exposures)
              Es<-Exposures@Exp
              signature_names<-Exposures@signames
              rownames(Ehat)<-signature_names
              model0<-!is.na(addata[[1]][1])
              if(model0){ 
                  cph0<-coxph(Surv(dtime,os)~., data=data.frame(dtime,os,addata))
              }else{
                  addata<-matrix(0,j,0) # empty matrix
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
                 y.max_anova<-max(LpvalAnova)
                 y.min_anova<-min(LpvalAnova)
                 Lp_quant_Anova<-quantile(LpvalAnova,invquant,na.rm=TRUE)
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
              boxnames<-signature_names
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
              multiv.tests$labels<-signature_names
              multiv.tests$colour <- rep(c("white", "gray95"), ceiling(nrow(multiv.tests)/2))[1:nrow(multiv.tests)]
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
                                   y_bgn=Lp_quant[sig_order],
                                   y_end=Lp_quant[sig_order],
                                   q1=0,q2=0,q3=0,q4=0,q5=0)
              g1<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                geom_boxplot(stat='identity',show.legend = FALSE,col=cor) +
                geom_hline(yintercept=lcut,col=col2) +
                geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                             data = segments, show.legend = FALSE)+
                theme_bw()+
                theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                labs(x="",y="-log(pvalue)")
              
              fp <- ggplot(multiv.tests, aes(x = HR, y = labels, xmin = Lower_CI, xmax = Upper_CI)) +
                geom_hline(aes(yintercept = labels, colour = colour), size = 7) + 
                geom_pointrange(shape = 22, fill = "black") +
                geom_vline(xintercept = 1, linetype = 3) +
                ylab("") +
                xlab("Hazard Ratio with 95% CI") +
                theme_classic() +
                scale_colour_identity() +
                scale_y_discrete(limits = rev(multiv.tests$labels)) +
                scale_x_log10(limits = c(0.25, 4), 
                              breaks = c(0.25, 0.5, 1, 2, 4), 
                              labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
                theme(axis.text.y = element_blank(), axis.title.y = element_blank())
              
              multiv.table<-data.frame(labels=multiv.tests[,"labels"],
                                       HR_CI=paste(round(multiv.tests$HR,3)," (",
                                                   round(multiv.tests$Lower_CI,3),"-",
                                                   round(multiv.tests$Upper_CI,3),")",sep=""),
                                       P.value=signif(multiv.tests$P.value,3),
                                       colour=multiv.tests$colour)
              
              ggtable <- ggplot(data = multiv.table, aes(y = labels)) +
                geom_hline(aes(yintercept = labels, colour = colour), size = 7) +
                geom_text(aes(x = 0, label = labels), hjust = 0) +
                geom_text(aes(x = 3, label = HR_CI)) +
                geom_text(aes(x = 3, y=n+0.5, label = "HR(CI)")) +
                geom_text(aes(x = 7, label = P.value), hjust = 1) +
                geom_text(aes(x = 7, y=n+0.5, label = "P.value"), hjust = 1) +
                scale_colour_identity() +
                theme_void() + 
                theme(plot.margin = margin(5, 0, 35, 0))
              
              forestplot <- ggarrange(plotlist = list(ggtable,fp),
                                  ncol = 2, nrow = 1)
              if(model0){
                md<-data.frame(Sig=rep("Anova",r),
                               Pvalues=as.vector(LpvalAnova))
                ms = group_by(md, Sig) %>% summarize(q1=min(Pvalues),
                                                     q2=quantile(Pvalues,p=0.25),
                                                     q3=median(Pvalues),
                                                     q4=quantile(Pvalues,p=0.75),
                                                     q5=max(Pvalues))
                sig_order<-order(signature_names)
                segments<-data.frame(Sig=signature_names[sig_order],x_bgn=0.6,x_end=1.4,
                                     y_bgn=Lp_quant_Anova[sig_order],
                                     y_end=Lp_quant_Anova[sig_order],
                                     q1=0,q2=0,q3=0,q4=0,q5=0)
                g2<-ggplot(ms, aes(x=Sig,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                  geom_boxplot(stat='identity',show.legend = FALSE,col='black') +
                  geom_hline(yintercept=lcut,col=col2) +
                  geom_segment(aes(x = x_bgn, y = y_bgn, xend = x_end, yend = y_end), col = col3,
                               data = segments, show.legend = FALSE)+
                  theme_bw()+
                  theme(axis.text.x=element_text(angle=0,vjust=.5,hjust=0,face="bold")) + 
                  labs(x="",y="-log(pvalue)")
                final_figure <- ggarrange(g1,forestplot,g2,ncol = 1)
              }else{
                final_figure <- ggarrange(g1,forestplot,ncol = 1)
              } 
              plot(final_figure)
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
    def=function(signexp_obj, Med_exp=NA, 
                 Clim=NA_integer_,
                 method.dist="euclidean", method.clust="fcm", relative=FALSE, 
                 m=2, plot_to_file=FALSE, file="FuzzyClustExp.pdf",colored=TRUE,
                 iseed=NA_integer_, try_all=FALSE,fast=TRUE,
                 parplan="multisession",...){
        standardGeneric("FuzzyClustExp")
    }
)
setMethod("FuzzyClustExp",signature(signexp_obj="SignExp", Med_exp="ANY", 
    Clim="ANY", method.dist="ANY", method.clust="ANY", 
    relative="ANY", m="ANY", plot_to_file="ANY",
    file="ANY",colored="ANY",iseed="ANY",try_all="ANY",
    fast="ANY",parplan="ANY"),
    function(signexp_obj, Med_exp, Clim, method.dist, method.clust,
        relative=FALSE, m, plot_to_file, file, colored=TRUE, iseed=NA, 
        try_all=FALSE,fast=TRUE,parplan){
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
                              relative, iseed=aseed)
                U<-Cm[[1]]
                PBMF<-as.vector(apply(Cm[[4]],3,function(E){PBMFindex(U,Data=E,m)}))
                return(list(median(PBMF),PBMF))
            },
            liminf=Cmin,limsup=Cmax,step=step0,significance=FALSE,parplan=parplan
            )
            bestn<-as.numeric(Ops[[1]])
            rm(Ops)
            cat(paste("Optimum number of groups is ",bestn,
                      ". Performing final clustering.\n",sep=""))
        }else{
            bestn<-Cmin
            cat(paste("Performing clustering in ",bestn," groups.\n",sep=""))
        }
        Cm<-CmeansExp(signexp_obj, Med_exp, C=bestn, method.dist, method.clust,
            relative, iseed=aseed,parplan=parplan)
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
        # Heatmap
        MF<-Cm$Meanfuzzy
        colnames(MF) = paste("Cl",1:NCOL(MF),sep="")
        rownames(MF)<-signexp_obj@samples
        show_samples<- j<30
        ph<-pheatmap(t(MF),border_color=NA, color=clr, 
                 clustering_method='ward.D2',
                 clustering_distance_cols='canberra',
                 cluster_rows = FALSE,
                 show_colnames=show_samples,
                 show_rownames=TRUE,
                 angle_col=90,
                 silent=TRUE)
        gt<-ph[[4]]
        ####################### ggplot2 boxplot
        AllF = Cm$AllFuzzy[ph[[2]]$order,,]
        rownames(AllF) = signexp_obj@samples[ph[[2]]$order]
        colnames(AllF) = paste("Cl",1:NCOL(MF),sep="")
        m = reshape2::melt(AllF)
        colnames(m)<-c("Sample","Cluster","r","value")
        m$Cluster<-as.factor(m$Cluster)
        ms = group_by(m, Sample, Cluster) %>% summarize(q1=min(value),
                                                           q2=quantile(value,p=0.25),
                                                           q3=median(value),
                                                           q4=quantile(value,p=0.75),
                                                           q5=max(value))
        ms_ord<-c()
        for(sp in rownames(AllF)){
          for(cl in colnames(AllF)){
              ms_ord<-c(ms_ord,which(ms$Sample==sp & ms$Cluster==cl))
          }
        }
        ms<-ms[ms_ord,]
        # rownames(m)<-paste(m$Sample,m$Cluster,sep="_")
        # comb_rownames<-paste(rep(rownames(AllF),each=NCOL(AllF)),rep(colnames(AllF),NROW(AllF)),sep="_")
        # m<-m[comb_rownames,]
        if(show_samples){
          g<-ggplot(ms, aes(x=Sample,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
            geom_boxplot(stat='identity') + 
            facet_grid(Cluster~.,scales="free") + 
            theme_cowplot() +
            theme(axis.text.x=element_text(size=6,vjust=.5,hjust=0,face="bold")) +
            scale_y_log10()+
            labs(x="")
        }else{
          g<-ggplot(ms, aes(x=Sample,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
            geom_boxplot(stat='identity') + 
            facet_grid(Cluster~.,scales="free") + 
            theme_cowplot() +
            theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
            scale_y_log10()+
            labs(x="")
        }
        #plotting 
        final_figure <- ggarrange(gt,g,ncol = 1)
        plot(final_figure)
        #######################
        if(plot_to_file){
            dev.off()
        }
        return(Cm)
    }
)

#CmeansExp
setGeneric("CmeansExp",
    def=function(signexp_obj, Med_exp=NA, C, method.dist="euclidean", 
        method.clust="fcm", relative=FALSE, iseed=NA_integer_,
        parplan="multisession",...){
        standardGeneric("CmeansExp")
    }
)
#For matrix:
setMethod("CmeansExp",signature(signexp_obj="matrix", Med_exp="ANY", C="ANY",
                                method.dist="ANY", method.clust="ANY", 
                                relative="ANY", iseed="ANY",parplan="ANY"),
          function(signexp_obj, Med_exp, C, method.dist, method.clust,
                   relative, plot_to_file, file, colored, iseed=NA, parplan){
            # initialize random seed
            if(!is.na(iseed)) set.seed(seed = iseed)
            de <- dim(signexp_obj) #[n,j]
            n<-de[[1]]; j<-de[[2]]
            if(is.na(Med_exp[1])){
              Med_exp<-signexp_obj
            }
            if(relative){ 
              signexp_obj<-t(t(signexp_obj)/colSums(signexp_obj))
              Med_exp<-t(t(Med_exp)/colSums(Med_exp))
              }
            colnames(Med_exp)<-colnames(signexp_obj)
            # if(method.clust=="km"){
            #     baseclust<-kmeans(t(Med_exp),centers=C)
            #     basefuzzy<-t(sapply(baseclust$cluster,function(n){
            #         as.numeric(c(1:C)==n)
            #     }))
            # }else{ 
            #     if(method.clust=="fcm"){
            #         baseclust<-ppclust::fcm(t(Med_exp),centers=C)
            #     }else if (method.clust=="pcm"){
            #         baseclust<-ppclust::pcm(t(Med_exp),centers=C)
            #     }else if (method.clust=="fpcm"){
            #         baseclust<-ppclust::fpcm(t(Med_exp),centers=C)
            #     }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
            #     basefuzzy<-baseclust$u
            # }
            if(method.clust=="km"){
              m <- 1        
            }else{ 
              if(method.clust=="fcm"){
                m <- 2
              }else stop("method.clust should be 'km' or 'fcm'!\n")
            }
            Es<-array(0,dim=c("n"=n,"j"=j,"r"=1))
            Es[,,1]<-signexp_obj
            Fuzzy<-FuzzyClusterCpp(Es,Med_exp,C,m,0.01)
            #return(list(Meanfuzzy=basefuzzy,AllFuzzy=basefuzzy))
            Meanfuzzy<-apply(Fuzzy[[1]],c(1,2),mean)
            return(list(Meanfuzzy=Meanfuzzy,
                        AllFuzzy=Fuzzy[[1]],
                        Centroids=Fuzzy[[2]],
                        Es=Es))
          })
#For SignExp
setMethod("CmeansExp",signature(signexp_obj="SignExp", Med_exp="ANY", C="ANY",
                                method.dist="ANY", method.clust="ANY", relative="ANY",iseed="ANY",parplan="ANY"),
          function(signexp_obj, Med_exp, C, method.dist, method.clust,
                   relative, iseed=NA,parplan){
            # initialize random seed
            if(!is.na(iseed)) set.seed(seed = iseed)
            if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
            dp <- dim(signexp_obj@Sign) #[i,n,r]
            de <- dim(signexp_obj@Exp) #[n,j,r]
            i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
            if(is.na(Med_exp[1])){
              Med_exp<-Median_exp(signexp_obj)
            }
            Es<-signexp_obj@Exp
            if(relative){ 
              Med_exp<-t(t(Med_exp)/colSums(Med_exp)) 
              colnames(Med_exp)<-signexp_obj@samples
              for(s in 1:r){
                thisE<-Es[,,s]
                Es[,,s]<-t(t(thisE)/colSums(thisE))
              }
            }
            if(method.clust=="km"){
              m <- 1        
            }else{ 
              if(method.clust=="fcm"){
                m <- 2
              }else stop("method.clust should be 'km' or 'fcm'!\n")
            }
            Fuzzy<-FuzzyClusterCpp(Es,Med_exp,C,m,0.01)
            
            
            #        if(method.clust=="km"){
            #            baseclust<-kmeans(t(Med_exp),centers=C)
            #            basefuzzy<-t(sapply(baseclust$cluster,function(n){
            #                as.numeric(c(1:C)==n)
            #            }))
            #        }else{ 
            #            if(method.clust=="fcm"){
            #                baseclust<-ppclust::fcm(t(Med_exp),centers=C)
            #            }else if (method.clust=="pcm"){
            #                baseclust<-ppclust::pcm(t(Med_exp),centers=C)
            #            }else if (method.clust=="fpcm"){
            #                baseclust<-ppclust::fpcm(t(Med_exp),centers=C)
            #            }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
            #            basefuzzy<-baseclust$u
            #        }
            #        avail.cores<-as.numeric(availableCores()-1)
            #        future::plan(parplan,workers=avail.cores)
            #        Fuzzy2<-apply(Es,3,function(Exposure){
            #            if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
            #            if(relative){ Exposure<-t(t(Exposure)/colSums(Exposure)) }
            #            colnames(Exposure)<-signexp_obj@samples
            #            if(method.clust=="km"){
            #                thisclust<-kmeans(t(Exposure),centers=C)
            #                thisfuzzy<-t(sapply(thisclust$cluster,function(n){
            #                    as.numeric(c(1:C)==n)
            #                rm(thisclust)
            #                }))
            #            }else{ 
            #                if(method.clust=="fcm"){
            #                    thisclust<-ppclust::fcm(t(Exposure),centers=C)
            #                }else if (method.clust=="pcm"){
            #                    thisclust<-ppclust::pcm(t(Exposure),centers=C)
            #                }else if (method.clust=="fpcm"){
            #                    thisclust<-ppclust::fpcm(t(Exposure),centers=C)
            #                }else stop("method.clust should be 'fcm', 'pcm' or 'fpcm'!\n")
            #                thisfuzzy<-thisclust$u
            #                rm(thisclust)
            #            }
            #            rm(Exposure)
            #            #hungarian algorithm to assign clusters
            #            D<-as.matrix(dist(rbind(t(basefuzzy),t(thisfuzzy)),"manhattan")) #L1 distance among clusters
            #            D<-D[1:C,(1:C)+C]#rows contain basefuzzy clusters, cols thisfuzzy.
            #            assignment <- clue::solve_LSAP(D)
            #            rm(D)
            #            thisfuzzy <- thisfuzzy[,as.vector(assignment)]
            #            rm(assignment)
            #            colnames(thisfuzzy)<-colnames(basefuzzy)
            #            return(thisfuzzy)
            #        })
            #        Fuzzy<-array(as.vector(Fuzzy2),dim=c('j'=j, 'c'=C, 'r'=r))
            #        rownames(Fuzzy)<-rownames(basefuzzy)
            #        colnames(Fuzzy)<-colnames(basefuzzy)
            
            
            Meanfuzzy<-apply(Fuzzy[[1]],c(1,2),mean)
            return(list(Meanfuzzy=Meanfuzzy,
                        AllFuzzy=Fuzzy[[1]],
                        Centroids=Fuzzy[[2]],
                        Es=Es))
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
        if(j<=30){
          colnames(Exposure)<-signexp_obj@samples
        }else{
          colnames(Exposure)<-NULL
        }
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

