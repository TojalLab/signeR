setGeneric("Paths",
    def=function(signexp_obj,plot_to_file=FALSE,
        file_suffix="plot.pdf",plots_per_page=4,...){
        standardGeneric("Paths")
    }
)
setMethod("Paths",signature(signexp_obj="SignExp",plot_to_file="ANY",
    file_suffix="ANY",plots_per_page="ANY"),
    function(signexp_obj,plot_to_file,file_suffix,
        plots_per_page,...){
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        cols_p<-terrain.colors(i)
        cols_e<-terrain.colors(j)
        if(plot_to_file){
            if(length(grep("\\.pdf$",file_suffix))==0){
                file_suffix<-paste(file_suffix,"pdf",sep=".")
            }
            P_file <- paste("Signature_paths",file_suffix,sep="_")
            E_file <- paste("Exposure_paths",file_suffix,sep="_")
            #P matrix plot
            pdf(file=P_file,width=7,height=2*plots_per_page)
            par(mfcol=c(plots_per_page,1),mar=c(3.8, 3, 1.9, 2) )
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=8, height=2*n)
            }
            par(mfcol=c(n,2),mar=c(3.8, 3, 1.9, 2) )
        }
        for (k in 1:n){
            if( k==1 | (k %% plots_per_page == 1 & plot_to_file)){
                maintitle <-"P matrix, signatures"
            }else{
                maintitle <-""
            }
            y.max <- max(signexp_obj@Sign[,k,])
            y.min <- min(signexp_obj@Sign[,k,])
            plot(c(0), c(0), type="l", xlim=c(0,r), ylim=c(y.min,y.max),
                xlab="Gibbs sampler iterations", ylab="", main=maintitle)
            for (s in 1:i){
                lines(signexp_obj@Sign[s, k, ], col=cols_p[s], ...)
            }
        }
        if(plot_to_file){
            dev.off()
            outmessage<-paste("Signature plots were exported to the file",
                P_file,"on the current directory.",sep=" ")
            cat(outmessage,"\n")
            #E matrix plot
            pdf(file=E_file,width=7,height=1.75*plots_per_page)
            par(mfcol=c(plots_per_page,1),mar=c(3.8, 3, 1.9, 2) )
        }
        for (k in 1:n){
            if( k==1 | (k %% plots_per_page == 1 & plot_to_file)){
                maintitle <-"E matrix, exposures"
            }else{
                maintitle <-""
            }
            y.max <- max(signexp_obj@Exp[k,,])
            y.min <- min(signexp_obj@Exp[k,,])
            plot(c(0), c(0), type="l", xlim=c(0,r), ylim=c(y.min,y.max),
                xlab="Gibbs sampler iterations", ylab="",
                main=maintitle)
            for (g in 1:j){
                lines(signexp_obj@Exp[k, g, ], col=cols_e[g], ...)
            }
        }
        if(plot_to_file){
            dev.off()
            outmessage2<-paste("Exposure plots were exported to the file",
                E_file,"on the current directory.",sep=" ")
            cat(outmessage2,"\n")
        }
    }
)

setGeneric("SignPlot",
    def=function(signexp_obj, plot_to_file=FALSE,
        file="Signature_plot.pdf",pal='bcr1',threshold=0,
        plots_per_page=4,gap=1,reord=NA_real_,...){
        standardGeneric("SignPlot")
    }
)
setMethod("SignPlot",signature(signexp_obj="SignExp",plot_to_file="ANY",
    file="ANY", pal="ANY", threshold="ANY",
    plots_per_page="ANY", gap="ANY",reord="ANY"),
    function(signexp_obj, plot_to_file, file, pal, threshold,
        plots_per_page, gap, reord,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        x <- c(1,16,32,48,64,80,96)
        rmut<-read.snv.context(signexp_obj@mutations)
        mutord<-order(rmut[[2]],rmut[[1]])
        x.names <- rmut[[1]][mutord] #triplets
        muts<-rmut[[2]][mutord]
        muttypes <- unique(rmut[[2]][mutord])
        if(is.na(reord[1])) reord<-1:n
        # define color palette
        if (pal == 'brew'){ xcolores <- c("#a6cee3","#1f78b4","#b2df8a",
            "#33a02c","#fb9a99","#e31a1c")
        }else if (pal =='lba'){ xcolores <- c("#29b4f4ff","#000000ff",
            "#f41d09ff","#bfc0bfff","#76d248ff","#f9b6b7ff")
        }else if (pal == 'bcr1'){ xcolores <- c("#006040FF","#BF7C40FF",
            "#BF0040FF","#406000FF","#4000BFFF","#FF7C00FF")
        }else if (pal == 'bcr2'){ xcolores <- c("#006040FF","#BF7C40FF",
            "#FF0000FF","#406000FF","#0000FFFF","#FF7C00FF")
        }else if (pal=="rdh"){
            xcolores<-colorRampPalette(c("#fee0d2","#a50f15"))(6)
        }else if (pal=="roh"){
            xcolores<-c("#fee0d2","#fcbba1","#fc9272","#fb6a4a",
                "#ef3b2c","#cb181d")
        }else if (pal=="blh"){
            xcolores<-colorRampPalette(c("#ece7f2","#d0d1e6","#a6bddb","#74a9cf",
                "#3690c0","#0570b0","#045a8d","#023858"))(6)
        }else if (pal=="bph"){
            xcolores<-rainbow(6,start=0.5,end=0.8)
        }else if (pal == 'bw'){ xcolores <- c("#949494FF", "#4D4D4DFF",
            "#FFFFFFFF","#B7B7B7FF","#DBDBDBFF", "#707070FF")
        }else stop("Unknown pallete.")
        if(i==96){
            bar.col <- rep(xcolores,each=16)
        }else{
            bar.col <- rep(xcolores,i,each=1)[1:i]
        }
        if (pal == 'bw'){
            border <- "black"
            top_border<-"black"
        }else{
            border <- bar.col
            top_border<-"white"
        }
        #Plot
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=7,height=1.75*plots_per_page)
        }else{
            plots_per_page<-n
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=1.75*n)
            }
        }
        hg<-rep(1,plots_per_page)
        hg[c(1,plots_per_page)]<-1.15
        for(k in 1:n){
            #get data
            Pdata<-signexp_obj@Psummary[mutord,reord[k],1:6,drop=TRUE]
            Pdata[Pdata<threshold]<-0
            Pdata<-data.frame(Pdata)
            Pdata$Sig<-signexp_obj@signames[reord[k]]
            Pdata$y.max <- max(Pdata[,6])
            Pdata$y.width <- Pdata$y.max*1.05
            Pdata<-Pdata[,-6]
            if (k==1){ 
                Pdata_all <- Pdata 
            }else{
                Pdata_all <- rbind(Pdata_all,Pdata)
            }
        }
        colnames(Pdata_all)<-c("q05","q25","medians","q75","q95","Sig","y.max","y.width")
        Pdata_all$triplets<-rep(x.names,n)
        Pdata_all$muts<-rep(muts,n)
        g<-ggplot(Pdata_all,aes(x=triplets,y=medians,fill=muts))+ 
          scale_fill_manual(values=xcolores) +
          geom_col(width=0.75) +
          facet_grid(Sig~muts,scales="free") +
          theme_classic(base_size = 15) +
          theme(axis.text.x=element_text(angle=90,size=6,vjust=.5,hjust=1,face="bold")) + 
          labs(x="SNV type",y="") +
          theme(legend.position="none") +
          coord_cartesian(expand=0) +
          geom_segment( y = Inf, yend = Inf, aes(color = muts), x = -Inf, xend = Inf, size = 4) +
          scale_color_manual(values=xcolores) +
          geom_errorbar(aes(ymin=q25, ymax=q75), width=.2, position=position_dodge(.9))+
          theme(strip.background.x=element_blank()) +
          theme(strip.background.y=element_rect(fill='white',color='black')) +
          theme(axis.ticks.x=element_blank())
        plot(g)
        if(plot_to_file){
            dev.off()
            outmess<-paste("Signature barplots were exported to the file",
                file,"on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)

setGeneric("ExposureBoxplot",
    def=function(signexp_obj, plot_to_file=FALSE,
        file="Exposure_boxplot.pdf", col='tan2', threshold=0,
        show_samples=NA_integer_, plots_per_page=4, reord=NA_real_,...){
        standardGeneric("ExposureBoxplot")
    }
)
setMethod("ExposureBoxplot",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", col="ANY", threshold="ANY", show_samples="ANY",
    plots_per_page="ANY",reord="ANY"),
    function(signexp_obj, plot_to_file, file, col, threshold,
        plots_per_page,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        bar.col <- col
        #Plot
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=7,height=1.75*plots_per_page)
        }else{
            plots_per_page<-n
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=1.75*n)
            }
        }
        if(is.na(reord[1])) reord<-1:n
        if(is.na(show_samples)){show_samples <- (j<=50)}
        ####################### ggplot2
        m = signexp_obj@Exp
        colnames(m) = signexp_obj@samples
        rownames(m) = signexp_obj@signames
        m = reshape2::melt(m)
        colnames(m)<-c("Signatures","Samples","r","value")
        m$Samples<-as.factor(m$Samples)
        m = group_by(m, Signatures, Samples) %>% summarize(q1=min(value),
                                                           q2=quantile(value,p=0.25),
                                                           q3=median(value),
                                                           q4=quantile(value,p=0.75),
                                                           q5=max(value))
        if(show_samples){
            g<-ggplot(m, aes(x=Samples,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
            geom_boxplot(stat='identity') + 
            facet_grid(Signatures~.,scales="free") + 
            theme_cowplot() +
            theme(axis.text.x=element_text(angle=90,size=8,vjust=.5,hjust=0,face="bold")) +
            scale_y_log10()+
            labs(x="")
        }else{
            g<-ggplot(m, aes(x=Samples,ymin=q1,lower=q2,middle=q3,upper=q4,ymax=q5)) + 
                geom_boxplot(stat='identity') + 
                facet_grid(Signatures~.,scales="free") + 
                theme_cowplot() +
                theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
                scale_y_log10()+
                labs(x="")
        }
        plot(g)
        #######################
        if(plot_to_file){
            dev.off()
            outmess<-paste("Exposure boxplots were exported to the file",
                file, "on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)

setGeneric("ExposureBarplot",
    def=function(signexp_obj, plot_to_file=FALSE,
        file="Exposure_barplot.pdf",col='tan2',threshold=0,relative=FALSE,
        title="", show_samples=NA_integer_, ...){
        standardGeneric("ExposureBarplot")
    }
)
setMethod("ExposureBarplot",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", col="ANY", threshold="ANY",relative="ANY",title="ANY",
    show_samples="ANY"),
    function(signexp_obj, plot_to_file, file, col, threshold,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
        bar.col <- col
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        Ehat<-Median_exp(signexp_obj)
        colnames(Ehat)<-signexp_obj@samples
        Ehat[Ehat<threshold]<-0
        if(relative){
            Ehat<-t(t(Ehat)/colSums(Ehat))
            ylabel<-"Relative Signature contribution to mutations"
        }else{
            ylabel<-"Signature contribution to mutation counts"
        }
        if(is.na(show_samples)){show_samples <- (j<=50)}
        mycolors<-terrain.colors(n)[n:1]
        #Plot
        if(plot_to_file){
            if(length(grep("\\.pdf$",file,perl=TRUE))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=7,height=7)
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=7)
            }
        }
        ################ggplot2
        m = Median_exp(signexp_obj)
        colnames(m) = signexp_obj@samples
        rownames(m) = signexp_obj@signames
        m = reshape2::melt(m)
        colnames(m)<-c("Signatures","Samples","value")
        if(relative){
          if(show_samples){
            g<-ggplot(m, aes(x=Samples,y=value,fill=Signatures)) + 
                geom_bar(stat='identity',position='fill') + 
                theme_cowplot()+
                theme(axis.text.x=element_text(angle=90,vjust=.5)) + 
                coord_cartesian(expand=0) + 
                scale_y_continuous(labels=percent_format()) +
                labs(x="",fill="Signatures")
            }else{
              g<-ggplot(m, aes(x=Samples,y=value,fill=Signatures)) + 
                geom_bar(stat='identity',position='fill') + 
                theme_cowplot()+
                theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) + 
                coord_cartesian(expand=0) + 
                scale_y_continuous(labels=percent_format()) +
                labs(x="",fill="Signatures")
            }
        }else{
          if(show_samples){
            g<-ggplot(m, aes(x=Samples,y=value,fill=Signatures)) + 
                geom_col(position='stack') + 
                theme_cowplot()+
                theme(axis.text.x=element_text(angle=90,vjust=.5)) + 
                coord_cartesian(expand=0) + 
                #scale_y_continuous(labels=percent_format()) + 
                labs(x="",fill="Signatures")
          }else{
            g<-ggplot(m, aes(x=Samples,y=value,fill=Signatures)) + 
              geom_col(position='stack') + 
              theme_cowplot()+
              theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) + 
              coord_cartesian(expand=0) + 
              #scale_y_continuous(labels=percent_format()) + 
              labs(x="",fill="Signatures")
          }
        }
        plot(g)
        #######################
        if(plot_to_file){
            dev.off()
            outmess<-paste("Exposure barplots were exported to the file",
                file, "on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)

setGeneric("SignHeat",
    def=function(signexp_obj, plot_to_file=FALSE,
        file="Signature_heatmap.pdf", nbins=50, pal="roh",...){
        standardGeneric("SignHeat")
    }
)
setMethod("SignHeat",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", nbins="ANY", pal="ANY"),
    function(signexp_obj, plot_to_file, file, nbins, pal,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        i<-dp[[1]]; n<-dp[[2]]; r<-dp[[3]]
        signat<-paste("S",1:n,sep="")
        rmut<-read.snv.context(signexp_obj@mutations)
        mutord<-order(rmut[[2]],rmut[[1]])
        x.names <- rmut[[1]][mutord]
        muttypes <- rmut[[2]][mutord]
        x.names.del<-unlist(lapply(strsplit(x.names,""),
            function(v){
                paste(v[c(1,3)],collapse=".")
            })
        )
        mut.names <- paste(muttypes,x.names.del,sep=":")
        if (pal == 'brew'){ colors <- colorRampPalette(c("#a6cee3","#1f78b4",
            "#b2df8a","#33a02c","#fb9a99","#e31a1c"))(nbins)
        }else if (pal =='lba'){ colors <- colorRampPalette(c("#29b4f4ff","#000000ff",
            "#f41d09ff","#bfc0bfff","#76d248ff","#f9b6b7ff"))(nbins)
        }else if (pal == 'bcr1'){ colors <- colorRampPalette(c("#006040FF","#BF7C40FF",
            "#BF0040FF","#406000FF","#4000BFFF","#FF7C00FF"))(nbins)
        }else if (pal == 'bcr2'){ colors <- colorRampPalette(c("#006040FF","#BF7C40FF",
            "#FF0000FF","#406000FF","#0000FFFF","#FF7C00FF"))(nbins)
        }else if (pal =="rdh"){ colors<-colorRampPalette(c("#fee0d2","#a50f15"))(nbins)
        }else if (pal =="roh"){ colors<-colorRampPalette(c("#fee0d2","#fcbba1",
                "#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15"))(nbins)
        }else if (pal == "blh"){ colors<-colorRampPalette(c("#ece7f2","#d0d1e6",
                "#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858"))(nbins)
        }else if (pal == "bph"){ colors<-rainbow(nbins,start=0.5,end=0.8)
        }else if (pal == "bw"){ colors<-colorRampPalette(c("#f0f0f0","#000000"))(nbins)
        }else stop("Unknown pallete.")
        Phat<-Median_sign(signexp_obj,normalize=TRUE)
        minval<-min(Phat)
        maxval<-max(Phat)
        lims<-seq(minval,maxval,length=nbins+1)
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=7,height=7)
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=7)
            }
        }
        ###################ggplot2
        m <- Median_sign(signexp_obj)[]
        rownames(m) <- signexp_obj@mutations
        colnames(m) <- signexp_obj@signames
        rmut<-read.snv.context(signexp_obj@mutations)
        mutord<-order(rmut[[2]],rmut[[1]])
        m = t(m[mutord,])
        triplets <- rmut[[1]][mutord]
        muts<-rmut[[2]][mutord]
        muttypes <- unique(rmut[[2]][mutord])
        mutannot<-data.frame(mutation=muttypes)
        rownames(mutannot)<-muttypes
        colnames(m) <- muts
        #clr=rev(colorRampPalette(brewer.pal(name="Spectral",n=11))(100))
        pheatmap(m,border_color=NA, color=colors, 
                 clustering_method='ward.D2',
                 clustering_distance_rows='canberra',
                 cluster_cols=F,annotation_col = mutannot,
                 labels_col=triplets,
                 fontsize_col=6)
        ##########################
        if(plot_to_file){
            dev.off()
            outmess<-paste("Signatures heatmap was exported to the file",
                file,"on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)

setGeneric("ExposureHeat",
    def=function(signexp_obj, plot_to_file=FALSE,
        file="Exposure_heatmap.pdf",nbins=50,pal="roh",distmethod="euclidean",
        clustermethod="complete",show_samples=NA_integer_,...){
        standardGeneric("ExposureHeat")
    }
)
setMethod("ExposureHeat",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", nbins="ANY", pal="ANY",distmethod="ANY",clustermethod="ANY",
    show_samples="ANY"),
    function(signexp_obj, plot_to_file, file, nbins, pal, distmethod,
        clustermethod,show_samples,...){
        de <- dim(signexp_obj@Exp) #[n,j,r]
        n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
        signat<-paste("S",1:n,sep="")
        if (pal == 'brew'){ colors <- colorRampPalette(c("#a6cee3","#1f78b4",
            "#b2df8a","#33a02c","#fb9a99","#e31a1c"))(nbins)
        }else if (pal =='lba'){ colors <- colorRampPalette(c("#29b4f4ff","#000000ff",
            "#f41d09ff","#bfc0bfff","#76d248ff","#f9b6b7ff"))(nbins)
        }else if (pal == 'bcr1'){ colors <- colorRampPalette(c("#006040FF","#BF7C40FF",
            "#BF0040FF","#406000FF","#4000BFFF","#FF7C00FF"))(nbins)
        }else if (pal == 'bcr2'){ colors <- colorRampPalette(c("#006040FF","#BF7C40FF",
            "#FF0000FF","#406000FF","#0000FFFF","#FF7C00FF"))(nbins)
        }else if (pal =="rdh"){ colors<-colorRampPalette(c("#fee0d2","#a50f15"))(nbins)
        }else if (pal =="roh"){ colors<-colorRampPalette(c("#fee0d2","#fcbba1",
                "#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15"))(nbins)
        }else if (pal == "blh"){ colors<-colorRampPalette(c("#ece7f2","#d0d1e6",
                "#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858"))(nbins)
        }else if (pal == "bph"){ colors<-rainbow(nbins,start=0.5,end=0.8)
        }else if (pal == "bw"){ colors<-colorRampPalette(c("#f0f0f0","#000000"))(nbins)
        }else stop("Unknown pallete.")
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        if(is.na(show_samples)){show_samples <- (j<=50)}
        Ehat<-Median_exp(signexp_obj)
        colnames(Ehat)<-signexp_obj@samples
        Dis<-dist(t(Ehat),method=distmethod)
        Hc <- hclust(Dis,method=clustermethod)
        maxh<-max(Hc$height)
        relh<-maxh/n
        plotM<-t(Ehat[,Hc$order,drop=FALSE])
        ddr <- as.dendrogram(Hc)
        minval<-log(min(plotM))
        maxval<-log(max(plotM))
        lims<-seq(minval,maxval,length=nbins+1)
        if(plot_to_file){
            if(length(grep("\\.pdf$",file))==0){
                file<-paste(file,"pdf",sep=".")
            }
            pdf(file,width=7,height=7)
        }else{
            if(!grepl("pdf|postscript|cairo_|png|tiff|jpeg|bmp",
                names(dev.cur()),perl=TRUE)){
                dev.new(width=7, height=7)
            }
        }
        ####################ggplot2
        m = Median_exp(signexp_obj)
        colnames(m) = signexp_obj@samples
        rownames(m) = signexp_obj@signames
        #clr=rev(colorRampPalette(brewer.pal(name="Spectral",n=11))(100))
        pheatmap(log(m), border_color=NA, color=colors, 
                 clustering_method='ward.D2',
                 clustering_distance_rows='canberra',
                 clustering_distance_cols='canberra',
                 show_colnames=show_samples)
        ###########################
        # layout(matrix(1:2,1,2),widths=c(4,1),heights=1)
        # ##Dendrogram##
        # par(mar=c(5, 1, 4, 4),las=2,cex=1,cex.axis=1)
        # plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",
        #     xlim=c(maxh,-1*maxh))
        # ##Heatmap##
        # for(s in 1:j){
        #     for(g in 1:n){
        #         color<-colors[max(sum(lims<log(plotM[s,g])),1)]
        #         rect(-1*g*relh,s-0.5,(1-g)*relh,s+0.5, col=color, border=color)
        #     }
        #     mtext(rownames(plotM)[s],side=4,at=s,cex=0.7)
        # }
        # for(g in 1:n){
        #     mtext(signat[g],side=1,at=(0.5-g)*relh,las=1,cex=0.8)
        # }
        # mtext("Exposure Heatmap",side=3,at=-1*maxh/2,cex=1,las=1,line=2)
        # mtext("Signatures",side=1,at=-1*maxh/2,cex=0.8,las=1,line=2)
        # ##Legend##
        # par(mar=c(3, 0, 3, 2),cex=0.7)#mai=c(1.6, 0.2, 1, 0.2)
        # plot(rep(1,nbins+1), lims, xlim=c(0,2), ylim=c(0,nbins+1),
        #     type="n", axes=FALSE, xlab="", ylab="")
        # for(k in 1:(nbins+1)){
        #     color<-colors[k]
        #     rect(0,k-1,0.8,k,col=color,border=color)
        #     text(1.1,k-1,format(lims[k],digits=4,scientific=FALSE),
        #         cex=0.9,pos=4)
        # }
        # mtext("log scale",side=1,at=1,las=1,cex=0.8,padj=0,line=-1)
        if(plot_to_file){
            dev.off()
            outmess<-paste("Exposures heatmap was exported to the file",
                file,"on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)
