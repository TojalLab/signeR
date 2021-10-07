#### S4 Class SignExp: pair of tensors with signatures and exposures ####
setClass("SignExp",
    slots = c(Sign="array",
        Exp="array",
        samples="character",
        mutations="character",
        sigSums="matrix",
        normalized="logical",
        Psummary="array",
        Esummary="array",
        Eoutliers="list"),
    prototype = list(Sign=array(NA,dim=c('i'=1,'n'=1,'k'=1)),
        Exp=array(NA,dim=c('n'=1,'j'=1,'k'=1)),
        samples=NA_character_,
        mutations=NA_character_,
        sigSums=matrix(NA_real_,1,1),
        normalized=FALSE,
        Psummary=array(NA,dim=c('i'=1,'n'=1,'q'=6)),
        Esummary=array(NA,dim=c('n'=1,'j'=1,'q'=6)),
        Eoutliers=list())
)

SignExpConstructor<-function(Ps=NA,Es=NA,samplenames=NA,mutnames=NA){
    if(!(is.na(Ps[1]) | is.na(Es[1]))){
        if(!(is.array(Ps) & is.array(Es))){
            stop("Signatures and exposures must be arrays.")
        }
        dp<-dim(Ps) #[i,n,r]
        de<-dim(Es) #[n,j,r]
        if(!(dp[[3]]==de[[3]] & dp[[2]]==de[[1]])){
            stop("Signatures and exposures are not compatible")
        }
        n<-dp[[2]]
        sSums<-apply(Ps,c(2,3),sum)
        if(all(is.na(samplenames))){
            samplenames<-paste("sample",1:de[[2]],sep="_")
        }
        if(all(is.na(mutnames))){
            bases=c("A","C","G","T")
            mutnames<-paste(rep(bases,each=4,6),rep(c("C","T"),each=48),
                rep(bases,24),">",
                c(rep(bases[-2],each=16),rep(bases[-4],each=16)),
                sep="")
        }
        SE<-new("SignExp",Sign=Ps,Exp=Es,samples=samplenames,mutations=mutnames,
            sigSums=sSums,normalized=FALSE,
            Psummary=array(NA,dim=c('i'=1,'n'=1,'q'=6)),
            Esummary=array(NA,dim=c('n'=1,'j'=1,'q'=6)),
            Eoutliers=list())
        SE<-Normalize(SE)
    }else{
        SE<-new("SignExp")
    }
    return(SE)
}

setGeneric("setSamples",
    def=function(signexp_obj,names){
        standardGeneric("setSamples")
    }
)
setMethod("setSamples",signature(signexp_obj="SignExp",names="ANY"),
    function(signexp_obj,names){
        signexp_obj@samples<-names
        return(signexp_obj)
    }
)

setGeneric("setMutations",
    def=function(signexp_obj,mutations){
        standardGeneric("setMutations")
    }
)
setMethod("setMutations",signature(signexp_obj="SignExp",mutations="ANY"),
    function(signexp_obj,mutations){
        signexp_obj@mutations<-mutations
        return(signexp_obj)
    }
)

setGeneric("Normalize",
    def=function(signexp_obj){
        standardGeneric("Normalize")
    }
)
setMethod("Normalize",signature("SignExp"), function(signexp_obj){
    #Normalize signatures.
    if(!signexp_obj@normalized){
        Ps<-signexp_obj@Sign
        Es<-signexp_obj@Exp
        dp<-dim(Ps)
        de<-dim(Es)
        i<-dp[[1]]
        n<-dp[[2]]
        j<-de[[2]]
        r<-de[[3]]
        for(k in 1:r){
            P<-Ps[,,k]
            E<-Es[,,k]
            vm<-signexp_obj@sigSums[,k,drop=TRUE]
            signexp_obj@Sign[,,k] <- t(t(matrix(as.vector(P),i,n))/vm)
            signexp_obj@Exp[,,k] <- matrix(as.vector(E),n,j)*vm
        }
        signexp_obj@normalized=TRUE
        Psum<-array(NA,dim=c('i'=i,'n'=n,'q'=6))
        Esum<-array(NA,dim=c('n'=n,'j'=j,'q'=5))
        Eout<-list()
        for(k in 1:n){
            P<-signexp_obj@Sign[,k,,drop=TRUE]
            Psum[,k,]<-t(apply(P,1,quantile,c(0.05,0.25,0.50,0.75,0.95,1)))
            E<- signexp_obj@Exp[k,,,drop=TRUE]
            bpstats<-apply(E,1,function(v){
                bp<-boxplot(v,plot=FALSE)
                return(list(bp$stats,bp$out))
            })
            Esum[k,,]<-t(matrix(unlist(lapply(bpstats,function(l){
                l[[1]]})),5,j))
            Eout[[k]]<-lapply(bpstats,function(l){l[[2]]})
        }
        signexp_obj@Psummary <- Psum
        signexp_obj@Esummary <- Esum
        signexp_obj@Eoutliers <- Eout
    }
    return(signexp_obj)
})

setGeneric("Reorder_signatures",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_signatures")
    }
)
setMethod("Reorder_signatures",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change signatures order.
        if(length(ord)==dim(signexp_obj@Sign)[[2]]){
            signexp_obj@Sign<-signexp_obj@Sign[,ord,]
            signexp_obj@Exp<-signexp_obj@Exp[ord,,]
            if(!all(is.na(signexp_obj@sigSums))){
                signexp_obj@sigSums<-signexp_obj@sigSums[ord,]
            }
            if(signexp_obj@normalized){
                signexp_obj@Psummary <- signexp_obj@Psummary[,ord,]
                signexp_obj@Esummary <- signexp_obj@Esummary[ord,,]
                signexp_obj@Eoutliers <- signexp_obj@Eoutliers[ord]
            }
            return(signexp_obj)
        }else{
            stop(paste("'ord' needs to be a vector of length",
                "equal to the number of signatures.",sep=" "))
        }
    }
)

setGeneric("Reorder_samples",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_samples")
    }
)
setMethod("Reorder_samples",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change sample order or take subsets.
        if(!length(ord)==dim(signexp_obj@Sign)[[1]]){
            warning("Reorder_samples will generate a new SignExp object",
                " with the sample subset enumerated in 'ord'.\n")
        }
        signexp_obj@Exp<-signexp_obj@Exp[,ord,,drop=FALSE]
        if(!all(is.na(signexp_obj@samples))){
            signexp_obj@samples<-signexp_obj@samples[ord]
        }
        if(signexp_obj@normalized){
            signexp_obj@Esummary <- signexp_obj@Esummary[,ord,,drop=FALSE]
            for(k in 1:length(signexp_obj@Eoutliers)){
                signexp_obj@Eoutliers[[k]] <- signexp_obj@Eoutliers[[k]][ord]
            }
        }
        return(signexp_obj)
    }
)

setGeneric("Reorder_mutations",
    def=function(signexp_obj,ord){
        standardGeneric("Reorder_mutations")
    }
)
setMethod("Reorder_mutations",signature(signexp_obj="SignExp",ord="numeric"),
    function(signexp_obj,ord){
        #Change mutation order.
        if(length(ord)==dim(signexp_obj@Sign)[[1]]){
            signexp_obj@Sign<-signexp_obj@Sign[ord,,]
            if(!all(is.na(signexp_obj@mutations))){
                signexp_obj@mutations<-signexp_obj@mutations[ord]
            }
            if(signexp_obj@normalized){
                signexp_obj@Psummary <- signexp_obj@Psummary[ord,,]
            }
            return(signexp_obj)
        }else{
            stop(paste("'ord' needs to be a vector of length",
                "equal to the number of mutations.",sep=" "))
        }
    }
)

setGeneric("Average_sign",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Average_sign")
    }
)
setMethod("Average_sign",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        Ps<-signexp_obj@Sign #[i,n,r]
        Phat<-apply(Ps,c(1,2),mean)
        return(Phat)
    }
)

setGeneric("Median_sign",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Median_sign")
    }
)
setMethod("Median_sign",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        i<-dp[[1]]; n<-dp[[2]]; r<-dp[[3]]
        Phat<-signexp_obj@Psummary[,,3,drop=TRUE]
        if(n==1) Phat<-matrix(as.vector(Phat),i,n)
        return(Phat)
    }
)

setGeneric("Average_exp",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Average_exp")
    }
)
setMethod("Average_exp",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        Es<-signexp_obj@Exp #[n,j,r]
        Ehat<-apply(Es,c(1,2),mean)
        return(Ehat)
    }
)

setGeneric("Median_exp",
    def=function(signexp_obj,normalize=TRUE){
        standardGeneric("Median_exp")
    }
)
setMethod("Median_exp",signature(signexp_obj="SignExp",normalize="ANY"),
    function(signexp_obj,normalize){
        if(normalize & !signexp_obj@normalized){
            signexp_obj<-Normalize(signexp_obj)
        }
        de <- dim(signexp_obj@Exp) #[n,j,r]
        n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
        Ehat<-signexp_obj@Esummary[,,3,drop=TRUE]
        if(n==1) Ehat<-matrix(as.vector(Ehat),n,j)
        return(Ehat)
    }
)

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
        x.names <- rmut[[1]][mutord]
        muttypes <- unique(rmut[[2]][mutord])
        if(is.na(reord[1])) reord<-1:n
        # define color palette
        if (pal == 'brew'){ xcolores <- c("#a6cee3","#1f78b4","#b2df8a",
            "#33a02c","#fb9a99","#e31a1c")
        }else if (pal =='lba'){ xcolores <- c("#29b4f4ff","#000000ff",
            "#f41d09ff","#bfc0bfff","#76d248ff","#f9b6b7ff")
        }else if (pal == 'bw'){ xcolores <- c("#949494FF", "#4D4D4DFF",
            "#FFFFFFFF","#B7B7B7FF","#DBDBDBFF", "#707070FF")
        }else if (pal == 'bcr1'){ xcolores <- c("#006040FF","#BF7C40FF",
            "#BF0040FF","#406000FF","#4000BFFF","#FF7C00FF")
        }else if (pal == 'bcr2'){ xcolores <- c("#006040FF","#BF7C40FF",
            "#FF0000FF","#406000FF","#0000FFFF","#FF7C00FF")
        }else stop("Unknown pallete.")
        bar.col <- rep(xcolores,each=16)
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
        par(cex= 1, cex.axis=1, mfrow=c(plots_per_page,1),
            mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA)
        hg<-rep(1,plots_per_page)
        hg[c(1,plots_per_page)]<-1.15
        layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
        for(k in 1:n){
            Pdata<-signexp_obj@Psummary[mutord,reord[k],1:6,drop=TRUE]
            Pdata[Pdata<threshold]<-0
            q05<-Pdata[,1]
            q25<-Pdata[,2]
            medians<-Pdata[,3]
            q75<-Pdata[,4]
            q95<-Pdata[,5]
            y.max <- max(Pdata[,6])
            y.width <- y.max*1.05
            bottom <- k==n | (k %% plots_per_page==0 & plot_to_file)
            top <- k==1 | (k %% plots_per_page==1 & plot_to_file)
            topmar<-1
            bottommar<-1
            if(bottom) bottommar<-3
            if(top){
                topmar<-3
                lastk <- n-plots_per_page*floor(k/plots_per_page)
                if (lastk<plots_per_page){
                    hg[plots_per_page]<-1
                    hg[lastk]<-1.15
                    layout(mat=matrix(1:plots_per_page,plots_per_page,1),
                        heights=hg)
                }
            }
            par(mar=c(bottommar,2,topmar,1),cex=0.6,cex.axis=0.6)
            mp <- barplot(medians, names=NULL, axes=FALSE, cex.names=1,
                cex.axis=0.8,las=2, col=bar.col, border=border,
                ylim=c(0,y.max), xlim=c(0,100*(1+gap)),
                space=c(0,rep(gap,95)),tcl=NA,family="mono",
                font=2)
            MP <- mp + (mp[2,] - mp[1,])/2 #positions
            if(top){ #big text mutations
                adj<-16*(1+gap)
                text(c(MP[8]+adj*(0:5)), c(rep(y.max+y.max*0.15, 6)),
                    labels=muttypes,
                    cex=1, col=c(rep("black",6)),family="mono",font=2)
            }
            axis(side=2,pos=-2,cex.axis=1,font=2)
            if(bottom){ #trinucleotides
                nameletters<-strsplit(x.names,"")
                for (s in 1:length(x.names)){
                    letters<-nameletters[[s]]
                    if(pal=="bw"){
                        middlecolor<-"black"
                    }else{
                        middlecolor<-xcolores[floor((s-1)/16)+1]
                    }
                    mtext(letters[1],side=1,line=0,at=mp[s],cex=0.55,
                        las=2,col="grey30",adj=3.8,family="mono")
                    mtext(letters[2],side=1,line=0,at=mp[s],cex=0.55,
                        las=2,col=middlecolor,font=2,adj=2.7,
                        family="mono")
                    mtext(letters[3],side=1,line=0,at=mp[s],cex=0.55,
                        las=2,col="grey30",adj=1.6,family="mono")
                }
            }
            #vertical lines
            segments(mp[1] - (mp[2,] - mp[1,])/2, 0,
                mp[1] - (mp[2,] - mp[1,])/2, y.max, lwd=0.4)
            for (j in 2:7) segments(MP[x[j]], 0, MP[x[j]], y.max, lwd=0.4)
            #top bars
            rect(mp[1] - (mp[2,] - mp[1,])/2, y.max, MP[x[2]], y.width,
                col=xcolores[1], border=top_border)
            for (j in 2:6) rect(MP[x[j]], y.max, MP[x[j+1]], y.width,
                col=xcolores[j], border=top_border)
            #Error bars
            segments(x0=mp,y0=q05,x1=mp,y1=q95,col="grey50",lwd=0.2)
            cap <- (mp[2,] - mp[1,])/6
            segments(x0=mp-(cap/3),y0=q05,x1=mp+(cap/3),y1=q05,
                col="grey40",lwd=0.2)
            segments(x0=mp-cap,y0=q25,x1=mp+cap,y1=q25,
                col="grey40",lwd=0.2)
            segments(x0=mp-cap,y0=q75,x1=mp+cap,y1=q75,
                col="grey40",lwd=0.2)
            segments(x0=mp-(cap/3),y0=q95,x1=mp+(cap/3),y1=q95,
                col="grey40",lwd=0.2)
            #side bars
            rect(97*(1+gap),0.01*y.max,100*(1+gap),0.99*y.max,
                col="grey90",border="grey10")
            text(98.5*(1+gap),y.max/2,labels=paste("S",k,sep=""),cex=0.9)
        }
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
        plots_per_page=4,reord=NA_real_,...){
        standardGeneric("ExposureBoxplot")
    }
)
setMethod("ExposureBoxplot",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", col="ANY", threshold="ANY",
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
        par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1),
            mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2)
        hg<-rep(1,plots_per_page)
        hg[plots_per_page]<-1.4
        layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
        for(k in 1:n){
            bp<-t(signexp_obj@Esummary[reord[k],,1:5,drop=TRUE])
            outliers<-signexp_obj@Eoutliers[[reord[k]]]
            y.max <- max(c(as.vector(bp),unlist(outliers)))
            y.width <- y.max+y.max*.05
            bp[bp<threshold]<-0
            q05<-bp[1,]
            q25<-bp[2,]
            medians<-bp[3,]
            q75<-bp[4,]
            q95<-bp[5,]
            bottom <- k==n | (k %% plots_per_page==0 & plot_to_file)
            top <- k==1 | (k %% plots_per_page==1 & plot_to_file)
            topmar<-1
            bottommar<-1
            if(bottom) bottommar<-5
            if(top){
                lastk <- n-plots_per_page*floor(k/plots_per_page)
                if (lastk<plots_per_page){
                    hg[plots_per_page]<-1
                    hg[lastk]<-1.4
                    layout(mat=matrix(1:plots_per_page,plots_per_page,1),
                        heights=hg)
                }
            }
            par(mar=c(bottommar,3.5,topmar,1),cex=0.7,cex.axis=1)
            plot(1:j, medians, ylim=c(0,y.max),xlim=c(0.6,j+1.4),
                type="n", main="", xlab="", ylab="", xaxt="n",font=2)
            boxplot(bp, col=bar.col, border="black", add=TRUE, cex=0.5,
                at=1:j, xaxt="n", lwd=0.2, pch=45)
            for (g in 1:j){
                points(rep(g,length(outliers[[g]])),outliers[[g]],
                    col="grey40", lwd=0.2, pch=45)
            }
            if (bottom){
                axis(side=1, at=1:j, labels=signexp_obj@samples, font=2)
            }
            #Error bars
            mp<-1:j
            segments(x0=mp,y0=q05,x1=mp,y1=q95,col="grey50",lwd=0.2)
            cap <- 1/6
            segments(x0=mp-(cap/3),y0=q05,
                x1=mp+(cap/3),y1=q05,col="grey40",lwd=0.2)
            segments(x0=mp-cap,y0=q25,
                x1=mp+cap,y1=q25,col="grey40",lwd=0.2)
            segments(x0=mp-cap,y0=q75,
                x1=mp+cap,y1=q75,col="grey40",lwd=0.2)
            segments(x0=mp-(cap/3),y0=q95,
                x1=mp+(cap/3),y1=q95,col="grey40",lwd=0.2)
            #side bars
            rect(j+0.6,0.01*y.max, j+1.4, 0.99*y.max, col="grey90",
                border="grey10")
            text(j+1, y.max/2, labels=paste("S",k,sep=""), cex=0.9)
        }
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
        title="", samplenames=TRUE, ...){
        standardGeneric("ExposureBarplot")
    }
)
setMethod("ExposureBarplot",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", col="ANY", threshold="ANY",relative="ANY",title="ANY",
    samplenames="ANY"),
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
        if(j>50){ samplenames<-FALSE }
        if(!samplenames){ colnames(Ehat)<-NULL }
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
        layout(mat=matrix(1:2,1,2),widths=c(5,1))
        par(cex= 1, cex.axis=0.9, mar=c(4,5,2,1), mgp=c(4,0.35,0), xpd=NA,las=2)
        barplot(Ehat[n:1,],col=mycolors,las=3,xlab="Samples",ylab="",main=title)
        mtext(ylabel,side=2,at=0.5*max(colSums(Ehat)),
            cex=1,las=3,padj=-3)
        par(mar=c(15,1,3,0))
        plot(rep(1,n),1:n,type="n",xlim=c(-2,4),ylim=c(-5,n),main="Signatures",
            xlab="",ylab="",cex.main=0.8,xaxt="n",yaxt="n",bty="n")
        for(k in 1:n){
            rect(0,k-1,1,k,col=mycolors[k],border="black")
            text(2,k-0.5,paste("S",n-k+1,sep=""),cex=0.6)

        }
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
        file="Signature_heatmap.pdf", nbins=20, pal="roh",...){
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
        if (pal=="rdh"){
            colors<-colorRampPalette(c("#fee0d2","#a50f15"))(nbins)
        }else if (pal=="roh"){
            colors<-colorRampPalette(c("#fee0d2","#fcbba1","#fc9272","#fb6a4a",
                "#ef3b2c","#cb181d","#a50f15"))(nbins)
        }else if (pal=="blh"){
            colors<-colorRampPalette(c("#ece7f2","#d0d1e6","#a6bddb","#74a9cf",
                "#3690c0","#0570b0","#045a8d","#023858"))(nbins)
        }else if (pal=="bph"){
            colors<-rainbow(nbins,start=0.5,end=0.8)
        }else if (pal=="bw"){
            colors<-colorRampPalette(c("#f0f0f0","#000000"))(nbins)
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
        layout(matrix(1:2,1,2),widths=c(4,1),heights=2)
        par(mar=c(2,1.2,2,1.2))
        plot(1:n, xlim=c(-0.3,n), ylim=c(-1,i), type="n", axes=FALSE,
            main="Signatures Heatmap", xlab="Signatures", ylab="")
        for(s in 1:i){
            for(r in 1:n){
                color<-colors[max(sum(lims<Phat[i+1-s,r]),1)]
                rect(r-1, s-1, r, s, col=color, border=color)
                if(s==1){
                    text(r-0.5,-1,signat[r],cex=0.6)
                }
            }
            text(-0.3,s+0.1,mut.names[i+1-s],cex=0.4)
        }
        par(mar=c(5,1,5,1))
        plot(rep(1,nbins+1), lims, xlim=c(0,2), ylim=c(0,nbins+1),
            type="n", axes=FALSE, xlab="", ylab="")
        for(i in 1:(nbins+1)){
            color<-colors[i]
            rect(1,i-1,2,i,col=color,border=color)
            text(0.5,i-1,format(lims[i],digits=3,scientific=TRUE),cex=0.5)
        }
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
        file="Exposure_heatmap.pdf",nbins=20,pal="roh",distmethod="euclidean",
        clustermethod="complete",...){
        standardGeneric("ExposureHeat")
    }
)
setMethod("ExposureHeat",signature(signexp_obj="SignExp", plot_to_file="ANY",
    file="ANY", nbins="ANY", pal="ANY",distmethod="ANY",clustermethod="ANY"),
    function(signexp_obj, plot_to_file, file, nbins, pal, distmethod,
        clustermethod,...){
        de <- dim(signexp_obj@Exp) #[n,j,r]
        n<-de[[1]]; j<-de[[2]]; r<-de[[3]]
        signat<-paste("S",1:n,sep="")
        if (pal=="rdh"){
            colors<-colorRampPalette(c("#fee0d2","#a50f15"))(nbins)
        }else if (pal=="roh"){
            colors<-colorRampPalette(c("#fee0d2","#fcbba1","#fc9272","#fb6a4a",
                "#ef3b2c","#cb181d","#a50f15"))(nbins)
        }else if (pal=="blh"){
            colors<-colorRampPalette(c("#ece7f2","#d0d1e6","#a6bddb","#74a9cf",
                "#3690c0","#0570b0","#045a8d","#023858"))(nbins)
        }else if (pal=="bph"){
            colors<-rainbow(nbins,start=0.5,end=0.8)
        }else if (pal=="bw"){
            colors<-colorRampPalette(c("#f0f0f0","#000000"))(nbins)
        }else stop("Unknown pallete.")
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
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
        layout(matrix(1:2,1,2),widths=c(4,1),heights=1)
        ##Dendrogram##
        par(mar=c(5, 1, 4, 4),las=2,cex=1,cex.axis=1)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",
            xlim=c(maxh,-1*maxh))
        ##Heatmap##
        for(s in 1:j){
            for(g in 1:n){
                color<-colors[max(sum(lims<log(plotM[s,g])),1)]
                rect(-1*g*relh,s-0.5,(1-g)*relh,s+0.5, col=color, border=color)
            }
            mtext(rownames(plotM)[s],side=4,at=s,cex=0.7)
        }
        for(g in 1:n){
            mtext(signat[g],side=1,at=(0.5-g)*relh,las=1,cex=0.8)
        }
        mtext("Exposure Heatmap",side=3,at=-1*maxh/2,cex=1,las=1,line=2)
        mtext("Signatures",side=1,at=-1*maxh/2,cex=0.8,las=1,line=2)
        ##Legend##
        par(mar=c(3, 0, 3, 2),cex=0.7)#mai=c(1.6, 0.2, 1, 0.2)
        plot(rep(1,nbins+1), lims, xlim=c(0,2), ylim=c(0,nbins+1),
            type="n", axes=FALSE, xlab="", ylab="")
        for(k in 1:(nbins+1)){
            color<-colors[k]
            rect(0,k-1,0.8,k,col=color,border=color)
            text(1.1,k-1,format(lims[k],digits=4,scientific=FALSE),
                cex=0.9,pos=4)
        }
        mtext("log scale",side=1,at=1,las=1,cex=0.8,padj=0,line=-1)
        if(plot_to_file){
            dev.off()
            outmess<-paste("Exposures heatmap was exported to the file",
                file,"on the current directory.",sep=" ")
            cat(outmess,"\n")
        }
    }
)

setGeneric("DiffExp",
    def=function(signexp_obj,labels, method=kruskal.test, contrast="all",
        quant=0.5, cutoff=0.05, p.adj="BH", plot_to_file=FALSE,
        file="Diffexp_boxplot.pdf",colored=TRUE,relative=FALSE,...){
        standardGeneric("DiffExp")
    }
)
setMethod("DiffExp",signature(signexp_obj="SignExp", labels="character",
    method="ANY", contrast="ANY", quant="ANY", cutoff="ANY", p.adj="ANY",
    plot_to_file="ANY", file="ANY", colored="ANY",relative="ANY"),
    function(signexp_obj, labels, method, contrast, quant, cutoff, plot_to_file,
        file, colored,...){
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
        Pval<-matrix(NA,n,r)
        MoreExp<-data.frame(matrix(NA,r,n))
        for (k in 1:r){
            Exposure <- signexp_obj@Exp[,,'r'=k]
            if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
            if(relative){ Exposure<-t(t(Exposure)/colSums(Exposure)) }
            Pval[,k]<-sapply(1:n,function(s){
                Test<-method(as.vector(Exposure[s,used]),used_labels,...)
                return(Test$p.value)
            })
            MoreExp[k,]<-sapply(1:n,function(m){
                vet<-as.vector(Exposure[m,used])
                mean_exp<-sapply(classes,function(cl){
                    mean(vet[used_labels==cl])
                })
                maxexp <- max(mean_exp)
                return(classes[which(mean_exp==maxexp)[1]])
            })
        }
        Lpval <- -1*log(Pval) #n x r
        rownames(Lpval)<-paste("S",1:n,sep="")
        y.min<-min(Lpval); y.max<-max(Lpval[Lpval<Inf])
        lcut <- -1*log(cutoff)
        invquant<-1-quant
        Lpmed<-apply(Lpval,1,quantile,invquant,na.rm=TRUE)
        cor<-rep("black",n)
        cor[Lpmed>=lcut]<-col1
        bigexp<-rep(NA,n)
        boxnames<-paste("S",1:n,sep="")
        boxlines<-rep(0.5,n)
        signif<-which(Lpmed>=lcut)
        for(k in signif){
            morevet<-MoreExp[,k]
            freqcl<-sapply(classes, function(cl){ sum(morevet==cl) })
            maxfreq <- max(freqcl)
            bigexp[k]<-classes[which(freqcl==maxfreq)[1]]
            boxnames[k]<-paste(boxnames[k],bigexp[k],sep="\n")
            boxlines[k]<-1.5
        }
        multicompare <- length(signif)>0 & nclasses>2
        if(multicompare){ #Multiple comparisons
            MCpv<-array(NA,dim=c(nclasses-1,nclasses-1,length(signif),r))
            for (k in 1:r){
                Exposure <- signexp_obj@Exp[,,'r'=k]
                if(n==1) Exposure <- matrix(as.vector(Exposure),n,j)
                for (i in 1:length(signif)){
                    s<-signif[i]
                    pkct<-kwAllPairsConoverTest(
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
        plot(1:n,rep(lcut,n),type="n",main="",xlab="",ylab="-log(pvalue)",
            xlim=c(0.5,n+0.5),ylim=c(y.min,y.max),xaxt="n",cex.lab=1.2)
        Lpval[Lpval==Inf]<-NA
        boxplot(data.frame(t(Lpval)),at=1:n,add=TRUE,border=cor,
            names=rep("",n),pch=45)
        mtext(boxnames,side=1,line=boxlines,at=1:n,cex=1,las=1)
        lines(x=c(0,n+1),y=rep(lcut,2),col=col2)
        for (k in 1:n){segments(k-0.39,Lpmed[k],k+0.39,Lpmed[k],col=col3,lwd=2)}
        if(multicompare){
            if(new_plot) dev.new(width=7, height=7)
            par(mfrow=c(length(signif),1),mar=c(3.1,4.2,3,2))
            for (i in 1:length(signif)){
                multicomp<-Allcomp[[i]]
                s<-signif[i]
                boxdata<-list()
                for(cl in classes){  boxdata<-c(boxdata,
                    list(as.vector(signexp_obj@Exp[s,labels==cl,]))) }
                names(boxdata)=classes
                boxplot(boxdata,main=paste("Signature",s),xaxt="n",pch=45)
                for(c in 1:nclasses){
                    mtext(classes[c],side=3,at=c,cex=0.7)
                    mtext(paste(multicomp[[c]],collapse=","),side=1,at=c,padj=1,
                        cex=0.8)
                }
            }
        }
        if(plot_to_file){
            dev.off()
            cat(paste("Differential exposure analysis results",
                "were plotted to the file",file,
                "on the current directory.",sep=" "),"\n")
        }
        Pmed<-apply(Pval,1,quantile,quant)
        signif <- Pmed<=cutoff
        mainresult<-data.frame(matrix(signif,1,n))
        colnames(mainresult)<-paste("S",1:n,sep="")
        result_list<-list(result=mainresult, Pvquant=Pmed, Pvalues=Pval,
            MostExposed=bigexp)
        if(multicompare){ result_list<-c(result_list,list(Differences=List_sig,
            MCPvalues=List_pval)) }
        return(result_list)
    }
)

setGeneric("Classify",
    def=function(signexp_obj, labels, method=knn, k=3, weights=NA,
        plot_to_file=FALSE, file="Classification_barplot.pdf",
        colors=NA_character_, min_agree=0.75,...){
        standardGeneric("Classify")
    }
)
setMethod("Classify",signature(signexp_obj="SignExp",labels="character",
    method="ANY", k="ANY", weights="ANY", plot_to_file="ANY",
    file="ANY", colors="ANY", min_agree="ANY"),
    function(signexp_obj, labels, method, k, plot_to_file, file, colors,
        min_agree,...){
        if(!signexp_obj@normalized) signexp_obj<-Normalize(signexp_obj)
        dp <- dim(signexp_obj@Sign) #[i,n,r]
        de <- dim(signexp_obj@Exp) #[n,j,r]
        i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
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
        if(identical(method,knn) & ntrain<k){
            stop(paste("The number of labeled samples is less then the number ",
                "of nearest neighbors used for classification, k=",k,
                ". Please run classification with k<=",ntrain,".",sep=""))
        }
        testsamples<-signexp_obj@samples[totest]
        cl<-labels[totrain]
        classes<-as.vector(levels(as.factor(cl)))
        nclass<-length(classes)
        if(all(is.na(weights))){
            weights<-rep(1,n)
        }
        if (length(weights)<n){
            weights<-rep(weights,n)[1:n]
        }
        Freqs<-matrix(0,nclass,ntest)
        rownames(Freqs)<-classes
        for (smp in 1:r){ # generates classifications
            D<-matrix(as.vector(signexp_obj@Exp[,,'r'=smp,drop=FALSE]),n,j)
            if (any(weights == "standardize")){
                cond <- which(weights == "standardize")
                weights[cond]<-1/apply(D[cond,],1,sd)
                weights<-as.numeric(weights)
            }
            D<-D*weights
            Train<-t(D[,!is.na(labels),drop=FALSE])
            Test<-t(D[,is.na(labels),drop=FALSE])
            if(identical(method,knn)){ # kNN classification
                Classific<-knn(Train,Test,cl,k=knearest)
                signexp_objclass<-as.vector(Classific)
            }else{
                Classific<-method(Train,Test,cl,...)
                signexp_objclass<-as.vector(Classific)
            }
            for (k in 1:ntest){
                currentcount <- Freqs[signexp_objclass[k],k]
                Freqs[signexp_objclass[k],k] <- currentcount+1
            }
        }
        result<-rep("",ntest)
        prob<-rep(0,ntest)
        for (k in 1:ntest){
            counts<-Freqs[,k]
            result[k] <- classes[which(counts==max(counts))][1]
            prob[k]   <- max(counts)/sum(counts)
        }
        result_final<-result
        result_final[prob < min_agree] <- "undefined"
        colnames(Freqs)<-paste(testsamples,result_final,sep="\n")
        names(result)<-testsamples
        names(prob)<-testsamples
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
        layout(matrix(1:2,1,2),widths=c(4,1),heights=1)
        bpl<-barplot(Freqs,main="Sample Classification",
            ylab="Frequencies",col=cols, las=1)
        for(k in 1:ntest){
            j<-which(classes==result[k])
            if(j==1){
                hei <- Freqs[1,k]/2
            }else{
                hei <- sum(Freqs[1:(j-1),k])+ Freqs[j,k]/2
            }
            text(bpl[k],hei,paste(round(100*prob[k],1),"%",sep=""))
        }
        par(mar=c(5,1,5,1))
        plot(1:4, 1:4, xlim=c(0,4), ylim=c(0,2*nclass), axes=FALSE,
            type="n", main="", xlab="", ylab="")
        for (k in 1:nclass){
            rect(0,nclass+k-2,1,nclass+k-1,col=cols[k])
            text(2,nclass+k-1.8,classes[k])
        }
        if(plot_to_file){
            dev.off()
            outmessage<-paste("Classification results were",
                "plotted to the file", file,
                "on the current directory.",sep=" ")
            cat(outmessage,"\n")
        }
        colnames(Freqs)<-paste(testsamples,result_final,sep="-")
        return(list(class=result_final,freq=prob,allfreqs=Freqs))
    }
)
