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
    #ggplot2
    reps<-unlist(lapply(Test_BICs,length))
    tested_n_rep<-c()
    for(k in 1:length(tested_n)){
        tested_n_rep<-c(tested_n_rep,rep(tested_n[k],reps[k]))
    }
    Boxdata<-data.frame(BIC=unlist(Test_BICs),n=tested_n_rep)
    p <- ggplot(Boxdata, aes(x=n, y=BIC, group=n)) + 
        scale_x_continuous(breaks=tested_n, labels=as.character(tested_n))+
        scale_y_continuous(breaks=ticks, labels=parse(text=labels))+
        geom_boxplot(outlier.shape = 19,notch=F)+
        labs(title="Bayesian Information Criterion",
            subtitle = "", 
            x="Number of processes", 
            y = "")+
        geom_vline(xintercept=nopt,lty=3)+
        theme_classic(base_size = 15)
    plot(p)
    if(plot_to_file){
        dev.off()
        outmessage<-paste("BIC boxplots were exported to the file",file,
            "on the current directory.",sep=" ")
        cat(outmessage,"\n")
    }
}
