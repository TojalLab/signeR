BICboxplot<-function(signeRout, plot_to_file=FALSE, 
    file="Model_selection_BICs.pdf"){
    nopt <- signeRout$Nsign
    tested_n <- signeRout$tested_n
    if(is.na(tested_n[1])){
        Test_BICs <- signeRout$BICs
        tested_n <- nopt
    }else{
        Test_BICs <- signeRout$Test_BICs
    }
    ymin <- min(unlist(Test_BICs))
    ymax <- max(unlist(Test_BICs))
    xmin <- min(tested_n)-0.6
    xmax <- max(tested_n)+0.6
    ticks <- seq(ymin,ymax,length=6)
    labels <- format(ticks, scientific=TRUE,digits=3)
    labels <- gsub("+0","+",labels,fixed=TRUE)
    labels <- gsub("-0","-",labels,fixed=TRUE)
    labels <- gsub("+","",labels,fixed=TRUE)
    labels <- gsub("e","%*%10^",labels)
    if(length(grep("\\.pdf$",file))==0){
        file<-paste(file,"pdf",sep=".")
    }
    if(plot_to_file){
        pdf(file,width=7,height=7)
    }
    par(mar=c(5,6,3,2))
    plot(tested_n, seq(ymin,ymax,length=length(tested_n)), xlim=c(xmin,xmax),
        ylim=c(ymin,ymax),type='n', main="Bayesian Information Criterion", 
        xlab="Number of processes",ylab="",xaxt="n",yaxt="n")
    axis(side=2,at=ticks,labels=parse(text=labels),tick=TRUE,las=1,cex=0.5)
    boxplot(Test_BICs, add=TRUE, at=tested_n, names=tested_n, 
        xaxt="s", yaxt="n",cex=0.7)
    if(length(tested_n)==1){
        axis(side=1,at=tested_n,labels=tested_n,las=1,cex=0.5)
    }
    abline(v=nopt,lty=3)
    if(plot_to_file){
        dev.off()
        outmessage<-paste("BIC boxplots were exported to the file",file,
            "on the current directory.",sep=" ")
        cat(outmessage,"\n")
    }
}
