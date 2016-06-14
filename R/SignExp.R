library(R.methodsS3)
library(R.oo)
library(class)
library(plot3D)
library(nloptr)

#### Class SignExp: pair of tensors with signatures and exposures ####
setConstructorS3("SignExp", function(Ps=NA,Es=NA,samplenames=NA,mutnames=NA){
  if(!(is.na(Ps[1]) | is.na(Es[1]))){
    if(is.list(Ps)){
      cat("Ps eh lista.\n")
      if(all(unlist(lapply(Ps,is.matrix)))){
        nrows<-unlist(lapply(Ps,NROW))
        ncols<-unlist(lapply(Ps,NCOL))
        if(all(nrows==nrows[1]) & all(ncols==ncols[1])){
          TEN<-to.tensor(0,c('i'=nrows[1],'n'=ncols[1],'r'=length(Ps)))
          for(k in 1:length(Ps)) TEN[,,'r'=k]<-Ps[[k]]
          Ps<-TEN
        }
      }
    }
    if(is.list(Es)){
      cat("Es eh lista.\n")
      if(all(unlist(lapply(Es,is.matrix)))){
        nrows<-unlist(lapply(Es,NROW))
        ncols<-unlist(lapply(Es,NCOL))
        if(all(nrows==nrows[1]) & all(ncols==ncols[1])){
          TEN<-to.tensor(0,c('n'=nrows[1],'j'=ncols[1],'r'=length(Es)))
          for(k in 1:length(Es)) TEN[,,'r'=k]<-Es[[k]]
          Es<-TEN
        }
      }
    }
    if(!(is.tensor(Ps) & is.tensor(Es))) stop("Signatures and exposures must be tensors or lists of matrices") 
    dp<-dim(Ps) #[i,n,r]
    de<-dim(Es) #[n,j,r]
    if(!(dp[[3]]==de[[3]] & dp[[2]]==de[[1]])) stop("Signatures and exposures are not compatible")
    sigSums<-matrix(0,de[[1]],de[[3]]) #[n,r]
    for(k in 1:de[[3]]){
      P<-Ps[,,k,drop=FALSE]
      sigSums[,k]<-colSums(P)
    }
    if(all(is.na(samplenames))) samplenames<-paste("sample",1:de[[2]],sep="_")
    if(all(is.na(mutnames))){
      bases=c("A","C","G","T")
      mutnames<-paste(rep(bases,each=4,6),rep(c("C","T"),each=48),rep(bases,24),">",c(rep(bases[-2],each=16),rep(bases[-4],each=16)),sep="")
    }
  }else{
    sigSums <- NA
    samplenames <- NA
    mutnames <- NA
  }
  extend(Object(), "SignExp",
         .Sign=Ps,
         .Exp=Es,
         .sigSums=sigSums,
         .samples=samplenames,
         .mutations=mutnames,
         .normalized=FALSE)
})

setMethodS3("Normalize",class="SignExp",function(this,...){
  #Normalize signatures.
  if(!this$.normalized){
    Ps<-this$.Sign
    Es<-this$.Exp
    dp<-dim(this$.Sign)
    de<-dim(this$.Exp)
    i<-dp[[1]]
    n<-dp[[2]]
    j<-de[[2]]
    r<-de[[3]]
    for(k in 1:r){
      P<-this$.Sign[,,k]
      E<-this$.Exp[,,k]
      vm<-this$.sigSums[,k,drop=TRUE]
      this$.Sign[,,k] <- t(t(matrix(as.vector(P),i,n))/vm)
      this$.Exp[,,k] <- matrix(as.vector(E),n,j)*vm
    }
    this$.normalized=TRUE
  }
})

setMethodS3("Reorder",class="SignExp",function(this,ord,...){
  #Change signatures order.
  if(length(ord)==dim(this$.Sign)[[2]]){
    this$.Sign<-this$.Sign[,ord,]
    this$.Exp<-this$.Exp[ord,,]
    if(!all(is.na(this$.sigSums))){
      this$.sigSums<-this$.sigSums[ord]
    }
  }else stop("'ord' needs to be a vector of length equal to signatures number.")
})

setMethodS3("Average_sign",class="SignExp",function(this,normalize=TRUE,...){
  if(normalize & !this$.normalized)  Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  i<-dp[[1]]; n<-dp[[2]] #r<-dp[[3]]
  Phat<-matrix(as.vector(mean.tensor(this$.Sign,along='r')),i,n)
  return(Phat)
})

setMethodS3("Median_sign",class="SignExp",function(this,normalize=TRUE,...){
  if(normalize & !this$.normalized)  Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  i<-dp[[1]]; n<-dp[[2]]; r<-dp[[3]]
  Phat<-matrix(0,i,n)
  for(k in 1:n){
    P<-matrix(as.vector(this$.Sign[,k,]),i,r)
    Phat[,k]<-apply(P,1,median)
  }
  return(Phat)
})

setMethodS3("Average_exp",class="SignExp",function(this,normalize=TRUE,...){
  if(normalize & !this$.normalized)  Normalize(this)
  de <- dim(this$.Exp) #[n,j,r]
  n<-de[[1]]; j<-de[[2]] #r<-de[[3]]
  Ehat<-matrix(as.vector(mean.tensor(this$.Exp,along='r')),n,j)
  return(Ehat)
})  

setMethodS3("Median_exp",class="SignExp",function(this,normalize=TRUE,...){
  if(normalize & !this$.normalized)  Normalize(this)
  dp <- dim(this$.Exp) #[n,j,r]
  n<-dp[[1]]; j<-dp[[2]]; r<-dp[[3]]
  Ehat<-matrix(0,n,j)
  for(k in 1:n){
    E<-matrix(as.vector(this$.Exp[k,,]),j,r)
    Ehat[k,]<-apply(E,1,median)
  }
  return(Ehat)
})

setMethodS3("Paths",class="SignExp",function(this,file_sufix="plot.png",plots_per_page=4,...){
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  if(length(grep("\\.pdf$",file_sufix))){
    filetype <- "pdf"
  }else{
    if(length(grep("\\.png$",file_sufix))==0) file_sufix<-paste(file_sufix,"png",sep=".") 
    filetype <- "png"
  }
  #P matrix plot
  cores<-gg2.col(i)
  if(filetype == "pdf"){
    pdf(file=paste("Signature_paths",file_sufix,sep="_"),width=8,height=2*plots_per_page)
    par(mfrow=c(plots_per_page,1))
  }else{
    if(n<=plots_per_page){
      png(filename=paste("Signature_paths",file_sufix,sep="_"),width=1000,height=200*n)
      par(mfrow=c(n,1))
    }else{
      png(filename=paste("Signature_paths","%03d",file_sufix,sep="_"),width=1000,height=200*plots_per_page)
      par(mfrow=c(plots_per_page,1))
    }
  }
  for (k in 1:n){
    if(k %% plots_per_page == 1){
      maintitle <-"P matrix, signatures"
    }else{
      maintitle <-""
    }
    y.max <- max(as.vector(this$.Sign[,k,]))
    y.min <- min(as.vector(this$.Sign[,k,]))
    plot(c(0), c(0), type="l", xlim=c(0,r), ylim=c(y.min,y.max),
         xlab="Gibbs sampler iterations", ylab="",
         main=maintitle)
    for (s in 1:i) lines(as.vector(this$.Sign[s, k, ]), col=cores[s], ...)
  }
  dev.off()
  #E matrix plot
  cores<-gg2.col(j)
  if(filetype == "pdf"){
    pdf(file=paste("Exposure_paths",file_sufix,sep="_"),width=8,height=2*plots_per_page)
    par(mfrow=c(plots_per_page,1))
  }else{
    if(n<=plots_per_page){
      png(filename=paste("Exposure_paths",file_sufix,sep="_"),width=1000,height=200*n)
      par(mfrow=c(n,1))
    }else{
      png(filename=paste("Exposure_paths","%03d",file_sufix,sep="_"),width=1000,height=200*plots_per_page)
      par(mfrow=c(plots_per_page,1))
    }
  }
  for (k in 1:n){
    if(k %% plots_per_page == 1){
      maintitle <-"E matrix, exposures"
    }else{
      maintitle <-""
    }
    y.max <- max(as.vector(this$.Exp[k,,]))
    y.min <- min(as.vector(this$.Exp[k,,]))
    plot(c(0), c(0), type="l", xlim=c(0,r), ylim=c(y.min,y.max),
         xlab="Gibbs sampler iterations", ylab="",
         main=maintitle)
    for (g in 1:j) lines(as.vector(this$.Exp[k,g, ]), col=cores[g], ...)
  }
  dev.off()
})

setMethodS3("SignPlot",class="SignExp",function(this,file="Signature_plot.pdf",pal='bcr1',
                                                threshold=0,plots_per_page=4,gap=1,reord=NA,...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  x <- c(1,16,32,48,64,80,96)
  if(is.null(this$.mutations)){
    rename<-TRUE
  }else if (all(is.na(this$.mutations))){
    rename<-TRUE
  }else rename<-FALSE
  if(rename){
    x.names <- c(rep(c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT",
                       "GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT"),3),
                 rep(c("ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT",
                       "GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"),3))
    muttypes<- c("C>A","C>G","C>T","T>A","T>C","T>G")
    mutord<-1:96
  }else{
    rmut<-read.mutation(this$.mutations)
    mutord<-order(rmut[[2]],rmut[[1]])
    x.names <- rmut[[1]][mutord]
    muttypes <- unique(rmut[[2]][mutord])
  }
  if(is.na(reord[1])){
    reord<-1:n
  }
  # define color palette
  if (pal == 'web') {
    bar.col <- c(rep("#009fffff",16),rep("#f14ab4ff",16),
                 rep("#005c30ff",16),rep("#ffe000ff",16),
                 rep("#4f429fff",16),rep("#ff7a00ff",16))
    border <- bar.col
    xcolores <- c("#009fffff","#f14ab4ff","#005c30ff","#ffe000ff",
                  "#4f429fff","#ff7a00ff")
    top_border<-"white"
  } else if (pal == 'brew1') {
    bar.col <- c(rep("#a6cee3",16),rep("#1f78b4",16),rep("#b2df8a",16),
                 rep("#33a02c",16),rep("#fb9a99",16),rep("#e31a1c",16))
    border <- bar.col
    xcolores <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                  "#fb9a99","#e31a1c")
    top_border<-"white"
  } else if (pal == 'brew2') {
    bar.col <- c(rep("#7fc97f",16),rep("#beaed4",16),rep("#fdc086",16),
                 rep("#ffff99",16),rep("#386cb0",16),rep("#f0027f",16)) 
    border <- bar.col
    xcolores <- c('#7fc97f','#beaed4','#fdc086','#ffff99',
                  '#386cb0','#f0027f')
    top_border<-"white"
  } else if (pal == 'brew3') {
    bar.col <- c(rep('#b2182b',16),rep('#ef8a62',16),rep('#fddbc7',16),
                 rep('#d1e5f0',16),rep('#67a9cf',16),rep('#2166ac',16))
    border <- bar.col
    xcolores <- c('#b2182b','#ef8a62','#fddbc7','#d1e5f0',
                  '#67a9cf','#2166ac')
    top_border<-"white"
  } else if (pal == 'alx') {
    bar.col <- c(rep('#29b4f4ff',16),rep('#000000ff',16),
                 rep('#f41d09ff',16),rep('#bfc0bfff',16),
                 rep('#76d248ff',16),rep('#f9b6b7ff',16))
    border <- bar.col
    xcolores <- c('#29b4f4ff','#000000ff','#f41d09ff','#bfc0bfff',
                  '#76d248ff','#f9b6b7ff')
    top_border<-"white"
  } else if (pal == 'gg2') {
    bar.col <- c(rep(gg2.col(6)[1],16),rep(gg2.col(6)[2],16),
                 rep(gg2.col(6)[3],16),rep(gg2.col(6)[4],16),
                 rep(gg2.col(6)[5],16),rep(gg2.col(6)[6],16))
    border <- bar.col
    xcolores <- gg2.col(6)   
    top_border<-"white"
  } else if (pal == 'gg') {
    bar.col <- c(rep(gg.col(6)[1],16),rep(gg.col(6)[2],16),
                 rep(gg.col(6)[3],16),rep(gg.col(6)[4],16),
                 rep(gg.col(6)[5],16),rep(gg.col(6)[6],16))
    border <- bar.col
    xcolores <- gg.col(6)   
    top_border<-"white"
  } else if (pal == 'bw') {
    bar.col <- c(rep("#949494FF",16), rep("#4D4D4DFF",16),
                 rep("#FFFFFFFF",16), rep("#B7B7B7FF",16),
                 rep("#DBDBDBFF",16), rep("#707070FF",16))
    xcolores <- c("#949494FF", "#4D4D4DFF", "#FFFFFFFF",
                  "#B7B7B7FF", "#DBDBDBFF", "#707070FF")
    border <- "black"
    top_border<-"black"
  } else if (pal == 'bcr1') {
    xcolores <- c("#006040FF", #c>a
                  "#BF7C40FF", #c>g
                  "#BF0040FF", #c>t
                  "#406000FF", #t>a
                  "#4000BFFF", #t>c
                  "#FF7C00FF") #t>g
    bar.col <- rep(xcolores,each=16)
    border <- bar.col #c(rep("#0000FFFF",48),rep("#CC0000FF",48))
    top_border<-"white"
  } else if (pal == 'bcr2') {
    xcolores <- c("#006040FF", #c>a
                  "#BF7C40FF", #c>g
                  "#FF0000FF", #c>t
                  "#406000FF", #t>a
                  "#0000FFFF", #t>c
                  "#FF7C00FF") #t>g
    bar.col <- rep(xcolores,each=16)
    border <- bar.col #c(rep("#0000FFFF",48),rep("#CC0000FF",48))
    top_border<-"white"
  }
  #Plot
  if(length(grep("\\.pdf$",file))){
    pdf(file,width=7,height=1.5*plots_per_page)
    par(cex= 1, cex.axis=1, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA)
  }else{
    if(length(grep("\\.png$",file))==0) file<-paste(file,"png",sep=".")
    if(n<=plots_per_page){
      png(file,width=1800,height=400*n)
      par(cex= 1.8, cex.axis=1.5, mfrow=c(n,1), mar=c(4.5,2,2,1), mgp=c(4,0.35,0), xpd=NA)
    }else{
      png(sub("\\.png$","_%03d.png",file),width=1800,height=400*plots_per_page)
      par(cex= 5, cex.axis=5, mfrow=c(plots_per_page,1), mar=c(4.5,2,2,1), mgp=c(4,2,1), xpd=NA)
    }
  }
  hg<-rep(1,plots_per_page)
  hg[1]<-1.15
  hg[plots_per_page]<-1.15
  layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
  for(k in 1:n){
    P<-matrix(as.vector(this$.Sign[mutord,reord[k],]),i,r)
    y.max <- max(P)
    y.width <- y.max+y.max*.05
    medians<-apply(P,1,median)
    q05<-apply(P,1,quantile,0.05)
    q25<-apply(P,1,quantile,0.25)
    q75<-apply(P,1,quantile,0.75)
    q95<-apply(P,1,quantile,0.95)
    medians[medians<threshold]<-0
    q05[q05<threshold]<-0
    q25[q25<threshold]<-0
    q75[q75<threshold]<-0
    q95[q95<threshold]<-0
    bottom <- (k %% plots_per_page==0 | k==n)
    top <- k %% plots_per_page==1
    topmar<-1
    bottommar<-1
    this.names<-NULL
    if(bottom){
      this.names<-x.names
      bottommar<-3
    }
    if(top){
      topmar<-3
      lastk <- n-plots_per_page*floor(k/plots_per_page)
      if (lastk<plots_per_page){
        hg[plots_per_page]<-1
        hg[lastk]<-1.15
        layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
      }
    }
    par(mar=c(bottommar,2,topmar,1),cex=0.6,cex.axis=0.6)
    mp <- barplot(medians, names=NULL, axes=FALSE, cex.names=1, cex.axis=0.8,las=2, 
                  col=bar.col, border=border, ylim=c(0,y.max), xlim=c(0,100*(1+gap)),
                  space=c(0,rep(gap,95)),tcl=NA,family="mono",font=2)
    MP <- mp + (mp[2,] - mp[1,])/2 #positions
    if(top){
      #big text mutations
      adj<-16*(1+gap)
      text(c(MP[8]+adj*(0:5)), c(rep(y.max+y.max*0.15, 6)), 
           labels=muttypes, 
           cex=1, col=c(rep("black",6)),family="mono",font=2)
    }
    axis(side=2,pos=-2,cex.axis=1,font=2)
    if(bottom){
      #trinucleotides
      nameletters<-strsplit(this.names,"")
      for (s in 1:length(this.names)){
        letters<-nameletters[[s]]
        if(pal=="bw"){
          middlecolor<-"black"
        }else{
          middlecolor<-xcolores[floor((s-1)/16)+1]
        }
        mtext(letters[1],side=1,line=0,at=mp[s],cex=0.55,las=2,col="grey30",adj=3.8,family="mono")
        mtext(letters[2],side=1,line=0,at=mp[s],cex=0.55,las=2,col=middlecolor,font=2,adj=2.7,family="mono")
        mtext(letters[3],side=1,line=0,at=mp[s],cex=0.55,las=2,col="grey30",adj=1.6,family="mono")
      }
    }
    #vertical lines
    segments(mp[1] - (mp[2,] - mp[1,])/2, 0, mp[1] - (mp[2,] - mp[1,])/2, 
             y.max, lwd=0.4)
    for (j in 2:7) segments(MP[x[j]], 0, MP[x[j]], y.max, lwd=0.4)
    #top bars
    rect(mp[1] - (mp[2,] - mp[1,])/2, y.max, MP[x[2]], y.width, 
         col=xcolores[1], border=top_border)
    for (j in 2:6) rect(MP[x[j]], y.max, MP[x[j+1]], y.width, 
                        col=xcolores[j], border=top_border)
    #Error bars
    segments(x0=mp,y0=q05,x1=mp,y1=q95,col="grey50",lwd=0.2)
    cap <- (mp[2,] - mp[1,])/6
    segments(x0=mp-(cap/3),y0=q05,x1=mp+(cap/3),y1=q05,col="grey40",lwd=0.2)
    segments(x0=mp-cap,y0=q25,x1=mp+cap,y1=q25,col="grey40",lwd=0.2)
    segments(x0=mp-cap,y0=q75,x1=mp+cap,y1=q75,col="grey40",lwd=0.2)
    segments(x0=mp-(cap/3),y0=q95,x1=mp+(cap/3),y1=q95,col="grey40",lwd=0.2)
    #side bars
    rect(97*(1+gap),0.01*y.max,100*(1+gap),0.99*y.max,col="grey90",border="grey10")
    text(98.5*(1+gap),y.max/2,labels=paste("S",k,sep=""),cex=0.9)
  }
  dev.off()
})

setMethodS3("SignBoxplot",class="SignExp",function(this,file="Signature_boxplot.pdf",pal='brew1',
                                                   threshold=0,plots_per_page=4,...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  x <- c(1,16,32,48,64,80,96)
  if(is.null(this$.mutations) | all(is.na(this$.mutations))){
    x.names <- c(rep(c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT",
                       "GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT"),3),
                 rep(c("ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT",
                       "GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"),3))
    muttypes<- c("C>A","C>G","C>T","T>A","T>C","T>G")
  }else{
    rmut<-read.mutation(this$.mutations)
    x.names <- rmut[[1]]
    muttypes <- unique(rmut[[2]])
  }
  # define color palette
  if (pal == 'web') {
    bar.col <- c(rep("#009fffff",16),rep("#f14ab4ff",16),
                 rep("#005c30ff",16),rep("#ffe000ff",16),
                 rep("#4f429fff",16),rep("#ff7a00ff",16))
    border <- bar.col
    xcolores <- c("#009fffff","#f14ab4ff","#005c30ff","#ffe000ff",
                  "#4f429fff","#ff7a00ff")
  } else if (pal == 'brew1') {
    bar.col <- c(rep("#a6cee3",16),rep("#1f78b4",16),rep("#b2df8a",16),
                 rep("#33a02c",16),rep("#fb9a99",16),rep("#e31a1c",16))
    border <- bar.col
    xcolores <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                  "#fb9a99","#e31a1c")
  } else if (pal == 'brew2') {
    bar.col <- c(rep("#7fc97f",16),rep("#beaed4",16),rep("#fdc086",16),
                 rep("#ffff99",16),rep("#386cb0",16),rep("#f0027f",16)) 
    border <- bar.col
    xcolores <- c('#7fc97f','#beaed4','#fdc086','#ffff99',
                  '#386cb0','#f0027f')
  } else if (pal == 'brew3') {
    bar.col <- c(rep('#b2182b',16),rep('#ef8a62',16),rep('#fddbc7',16),
                 rep('#d1e5f0',16),rep('#67a9cf',16),rep('#2166ac',16))
    border <- bar.col
    xcolores <- c('#b2182b','#ef8a62','#fddbc7','#d1e5f0',
                  '#67a9cf','#2166ac')
  } else if (pal == 'alx') {
    bar.col <- c(rep('#29b4f4ff',16),rep('#000000ff',16),
                 rep('#f41d09ff',16),rep('#bfc0bfff',16),
                 rep('#76d248ff',16),rep('#f9b6b7ff',16))
    border <- bar.col
    xcolores <- c('#29b4f4ff','#000000ff','#f41d09ff','#bfc0bfff',
                  '#76d248ff','#f9b6b7ff')
  } else if (pal == 'gg2') {
    bar.col <- c(rep(gg2.col(6)[1],16),rep(gg2.col(6)[2],16),
                 rep(gg2.col(6)[3],16),rep(gg2.col(6)[4],16),
                 rep(gg2.col(6)[5],16),rep(gg2.col(6)[6],16))
    border <- bar.col
    xcolores <- gg2.col(6)   
  } else if (pal == 'gg') {
    bar.col <- c(rep(gg.col(6)[1],16),rep(gg.col(6)[2],16),
                 rep(gg.col(6)[3],16),rep(gg.col(6)[4],16),
                 rep(gg.col(6)[5],16),rep(gg.col(6)[6],16))
    border <- bar.col
    xcolores <- gg.col(6)   
  } else if (pal == 'bw') {
    bar.col <- c(rep("white",6))
    xcolores <- "white"
    border <- "black"
  }
  y.max <- max(this$.Sign)
  y.width <- y.max+y.max*.05
  #Plot
  if(length(grep("\\.pdf$",file))){
    pdf(file)#,width=8,height=2*plots_per_page)
    par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2, bty="n")
  }else{
    if(length(grep("\\.png$",file))==0) file<-paste(file,"png",sep=".")
    if(n<=plots_per_page){
      png(file,width=1800,height=360*n)
      par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2, bty="n")
    }else{
      png(sub("\\.png$","_%03d.png",file),width=1800,height=360*plots_per_page)
      par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2, bty="n")
    }
  }
  for(k in 1:n){
    P<- data.frame(t(matrix(as.vector(this$.Sign[,k,]),i,r)))
    medians<-apply(P,2,median)
    colnames(P)<-x.names
    plot(1:96,medians,ylim=c(0,y.max),xlim=c(0,100),type="n",main="",xlab="",ylab="",xaxt="n")
    #vertical lines
    mp<-c(0.5,(1:6)*16+0.5)
    segments(mp, 0, mp,y.max, lwd=0.4)
    #top bars
    rect(mp[-7] , y.max, mp[-1], y.width, 
         col=xcolores, border="white")
    boxplot(P, names=x.names, col=bar.col, border="black", add=TRUE, cex=0.2, at=1:96,xaxt="s",lwd=0.2) 
    #big text mutations
    text(mp[1:6]+8, rep(y.max*1.15, 6), 
         labels=muttypes, 
         cex=0.8, col=c(rep("black",6)))
    rect(97,0.01*y.max,100,0.99*y.max,col="grey90",border="grey10")
    text(98.5,y.max/2,labels=paste("S",k,sep=""),cex=0.9)
  }
  dev.off()
})

setMethodS3("SignHeat",class="SignExp",function(this,file="Signature_heatmap.pdf",nbins=20,...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  signat<-paste("S",1:n,sep="")
  m.names <- paste(rep(c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T",
                         "G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T"),6),
                   rep(c("C>A","C>G","C>T","T>A","T>C","T>G"),each=16),sep=":")
  Phat<-Average_sign(this,normalize=TRUE)
  colors<-rainbow(nbins,start=0.5,end=0.8)
  minval<-min(Phat)
  maxval<-max(Phat)
  lims<-seq(minval,maxval,length=nbins+1)
  if(length(grep("\\.pdf$",file))){
    pdf(file)
  }else{
    if(length(grep("\\.png$",file))==0) file<-paste(file,"png",sep=".")
    png(file,width=1800,height=3000)
  }
  layout(matrix(1:2,1,2),widths=c(4,1),heights=2)
  par(mar=c(2,1.2,2,1.2))
  plot(1:n,xlim=c(-0.3,n),ylim=c(-1,96),type="n",axes=FALSE,main="signeR Heatmap",xlab="Signatures",ylab="")
  for(s in 1:i){
    for(r in 1:n){
      cor<-colors[max(sum(lims<Phat[97-s,r]),1)]
      rect(r-1,s-1,r,s,col=cor,border=cor)
      if(s==1){
        text(r-0.5,-1,signat[r],cex=0.6)
      }
    }
    text(-0.3,s+0.1,m.names[97-s],cex=0.4) 
  }
  par(mar=c(5,1,5,1))
  plot(rep(1,nbins+1),lims,xlim=c(0,2),ylim=c(0,nbins+1),type="n",axes=FALSE,xlab="",ylab="")
  for(i in 1:(nbins+1)){
    cor<-colors[i]
    rect(1,i-1,2,i,col=cor,border=cor)
    text(0.5,i-1,round(lims[i],2),cex=0.5)
  }
  dev.off()
})

setMethodS3("Classify",class="SignExp",function(this,labels,method="knn",k=3,plotfile="Classification_barplot.pdf",...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  knearest<-k
  if(!length(labels)==j){
    stop("Labels should be provided for each sample. Samples to be classified should be labeled with 'NA'.")
  }
  totest<-is.na(labels)
  totrain<-!totest
  ntest<-sum(totest)
  testsamples<-this$.samples[totest]
  cl<-labels[totrain]
  classes<-as.vector(levels(as.factor(cl)))
  nclass<-length(classes)
  Freqs<-matrix(0,nclass,ntest)
  rownames(Freqs)<-classes
  for (smp in 1:r){ # generates classifications
    D<-matrix(as.vector(this$.Exp[,,'r'=smp,drop=FALSE]),n,j)
    Train<-t(D[,!is.na(labels)])
    Test<-t(D[,is.na(labels)])
    if(length(grep("knn",method,ignore.case = TRUE))){ # kNN classification
      Classific<-knn(Train,Test,cl,k=knearest)
      thisclass<-as.vector(Classific)
      for (k in 1:ntest){
        Freqs[thisclass[k],k]<-Freqs[thisclass[k],k]+1
      }
    }else if(length(grep("rf",method,ignore.case = TRUE))){ # Random Forest classification
      RF<-randomForest(x=Train, y=as.factor(cl), ntree=1000,importance=TRUE, keep.forest=TRUE)
      Classific<-predict(RF,newdata=Test,type="vote",norm.votes=FALSE,predict.all=TRUE)
      Freqs<-Freqs+Classific$aggregate[,classes]
    }else{
      Classific<-method(Train,Test,cl,...)
      thisclass<-as.vector(Classific)
      for (k in 1:ntest){
        Freqs[thisclass[k],k]<-Freqs[thisclass[k],k]+1
      }
    }
  }
  result<-rep("",ntest)
  prob<-rep(0,ntest)
  for (k in 1:ntest){
    counts<-Freqs[,k]
    result[k] <- classes[which(counts==max(counts))][1]
    prob[k]   <- max(counts)/sum(counts)
  }
  colnames(Freqs)<-paste(testsamples,result,sep="\n")
  names(result)<-testsamples
  names(prob)<-testsamples
  cores<-rainbow(nclass,start=0.5,end=0.8)
  if(!is.na(plotfile)){
    if(length(grep("\\.pdf$",plotfile))){
      pdf(file=plotfile)
    }else{
      if(length(grep("\\.png$",plotfile))==0) plotfile<-paste(plotfile,"png",sep=".") 
      png(filename=plotfile,width=800,height=800)
    }
    layout(matrix(1:2,1,2),widths=c(4,1),heights=1)
    bpl<-barplot(Freqs,main="Sample Classification", ylab="Frequencies",col=cores)
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
    plot(1:4,1:4,xlim=c(0,4),ylim=c(0,2*nclass),axes=FALSE,type="n",main="",xlab="",ylab="")
    for (k in 1:nclass){
      rect(0,nclass+k-2,1,nclass+k-1,col=cores[k])
      text(2,nclass+k-1.8,classes[k])
    }
    dev.off()
  }
  colnames(Freqs)<-paste(testsamples,result,sep="-")
  return(list(class=result,freq=prob,allfreqs=Freqs))
})

setMethodS3("DiffExp",class="SignExp",function(this,labels,method="kw",contrast="all",quant=0.5,
                                               cutoff=0.05,plotfile="Diffexp_boxplot.pdf",pal="color",...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  cl<-labels[!is.na(labels)]
  if(all(contrast=="all")){
    classes <- as.vector(levels(as.factor(cl)))
  }else{
    classes <- contrast
  }
  used <- labels %in% classes
  used_labels<-as.factor(as.vector(labels[used]))
  if(pal=="color"){
    col1<-"darkgreen"
    col2<-"red"
    col3<-"blue"
  }else{
    col1<-"black"
    col2<-"black"
    col3<-"black"
  }
  Pval<-matrix(NA,n,r)
  MoreExp<-data.frame(matrix(NA,r,n))
  for (k in 1:r){
    Exposure <- matrix(as.vector(this$.Exp[,,'r'=k,drop=FALSE]),n,j)
    Pval[,k]<-sapply(1:n,function(m){
      if (method=="kw"){
        Test<-kruskal.test(as.vector(Exposure[m,used]),used_labels)
        return(Test$p.value)
      }else{
        Test<-method(as.vector(Exposure[m,used]),used_labels,...)
        return(Test$p.value)
      }
    })
    MoreExp[k,]<-sapply(1:n,function(m){
      vet<-as.vector(Exposure[m,used])
      totexp<-sapply(classes,function(cl){
        sum(vet[used_labels==cl])
      })
      maxexp <- max(totexp)
      return(classes[which(totexp==maxexp)[1]])
    })
  }
  Lpval <- -1*log(Pval) #n x r
  rownames(Lpval)<-paste("S",1:n,sep="")
  y.min<-min(Lpval)
  y.max<-max(Lpval)
  lcut <- -1*log(cutoff)
  invquant<-1-quant
  Lpmed<-apply(Lpval,1,quantile,invquant)
  cor<-rep("black",n)
  cor[Lpmed>=lcut]<-col1 #########
  bigexp<-rep(NA,n)
  boxnames<-paste("S",1:n,sep="")
  boxlines<-rep(0.5,n)
  for(k in which(Lpmed>=lcut)){
    morevet<-MoreExp[,k]
    freqcl<-sapply(classes,function(cl){
      sum(morevet==cl)
    })
    maxfreq <- max(freqcl)
    bigexp[k]<-classes[which(freqcl==maxfreq)[1]]
    boxnames[k]<-paste(boxnames[k],bigexp[k],sep="\n")
    boxlines[k]<-1.5
  }
  if(!is.na(plotfile)){
    if(length(grep("\\.pdf$",plotfile))){
      pdf(file=plotfile,width=5,height=5)
    }else{
      if(length(grep("\\.png$",plotfile))==0) plotfile<-paste(plotfile,"png",sep=".") 
      png(filename=plotfile,width=800,height=800)
    }
    par(mar=c(3.1,4.2,2,2))
    plot(1:n,rep(lcut,n),type="n",main="",xlab="",
         ylab="-log(pvalue)",xlim=c(0.5,n+0.5),ylim=c(y.min,y.max),xaxt="n",cex.lab=1.2)
    boxplot(data.frame(t(Lpval)),at=1:n,add=TRUE,border=cor,names=rep("",n),pch=45)
    mtext(boxnames,side=1,line=boxlines,at=1:n,cex=1)
    abline(h=-1*log(c(0.05,0.01,0.001)),lty=2)
    #abline(h=lcut,col=col2)
    #mtext("-log(0.01)",side=4,at=-log(0.01),cex=0.7,las=1)
    #mtext("-log(0.05)",side=4,at=-log(0.05),cex=0.7,las=1)
    #mtext("-log(0.001)",side=4,at=-log(0.001),cex=0.7,las=1)
    for (k in 1:n){
      segments(k-0.39,Lpmed[k],k+0.39,Lpmed[k],col=col3,lwd=2)
    }
    dev.off()
  }
  Pmed<-apply(Pval,1,quantile,quant)
  signif <- Pmed<=cutoff
  mainresult<-data.frame(matrix(signif,1,n))
  colnames(mainresult)<-paste("S",1:n,sep="")
  return(list(result=mainresult,Pvquant=Pmed,Pvalues=Pval,MostExposed=bigexp))
})

## support function

read.mutation<-function(stvet){
  triplets<-c()
  mutations<-c()
  for (st in stvet){
    if (grepl("[ACGT]{3,3}",st,ignore.case=TRUE,perl=TRUE)){
      strt<-regexpr("[ACGT]{3,3}",st)[[1]]
      trivet<-strsplit(st,"")[[1]][strt:(strt+2)]
      triplet<-paste(trivet,collapse="")
      st<-sub(triplet,"",st,ignore.case=TRUE)
      st<-sub(trivet[2],"",st,ignore.case=TRUE)
      strt<-regexpr("[ACGT]",st,ignore.case=TRUE)[[1]]
      newbase<-strsplit(st,"")[[1]][strt]
      mut<-paste(trivet[2],newbase,sep=">")
    }else if(grepl("[ACTG]>[ACTG]",st,ignore.case=TRUE,perl=TRUE)){
      strt<-regexpr("[ACTG]>[ACTG]",st,ignore.case=TRUE)[[1]]
      mutvet<-strsplit(st,"")[[1]][strt:(strt+2)]
      mut<-paste(mutvet,collapse="")
      st<-sub(mut,"",st,ignore.case=TRUE)
      splitvet<-strsplit(st,"")[[1]]
      splitvet<-splitvet[splitvet %in% c("A","C","G","T","a","c","g","t")]
      firstbase<-splitvet[1]
      lastbase<-splitvet[2]
      triplet<-paste(c(firstbase,mutvet[1],lastbase),collapse="")
    }else if (grepl("[ACTG]\\.[ACTG]",st,ignore.case=TRUE,perl=TRUE)){
      strt<-regexpr("[ACTG]\\.[ACTG]",st,ignore.case=TRUE)[[1]]
      trivet<-strsplit(st,"")[[1]][strt:(strt+2)]
      st<-sub(paste(trivet,collapse=""),"",st,ignore.case=TRUE)
      splitvet<-strsplit(st,"")[[1]]
      splitvet<-splitvet[splitvet %in% c("A","C","G","T","a","c","g","t")]
      oldbase<-splitvet[1]
      newbase<-splitvet[2]
      trivet[2]<-oldbase
      triplet<-paste(trivet,collapse="")
      mut<-paste(oldbase,newbase,sep=">")
    }else{
      triplet<-st
      mut<-""
    }
    triplets<-c(triplets,triplet)
    mutations<-c(mutations,mut)
  }
  return(list(triplets,mutations))
}

BICboxplot<-function(signeRout,file="Model_selection_BICs.pdf"){
  nopt <- signeRout$Nsign
  tested_n <- signeRout$tested_n
  Test_BICs <- signeRout$Test_BICs
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
  if(length(grep("\\.pdf$",file))){
    pdf(file,width=7,height=7)
  }else{
    if(length(grep("\\.png$",file))==0) file<-paste(file,"png",sep=".") 
    png(file,width=800,height=700)
  }
  par(mar=c(5,6,3,2))
  plot(tested_n,seq(ymin,ymax,length=length(tested_n)),xlim=c(xmin,xmax),ylim=c(ymin,ymax),type='n',main="Bayesian Information Criterion",
       xlab="Number of processes",ylab="",xaxt="n",yaxt="n")
  axis(side=2,at=ticks,labels=parse(text=labels),tick=TRUE,las=1,cex=0.5)
  boxplot(Test_BICs,add=TRUE,at=tested_n,names=tested_n,xaxt="s",yaxt="n",cex=0.7)
  abline(v=nopt,lty=3)
  dev.off()
}

setMethodS3("ExposureBoxplot",class="SignExp",function(this,file="Exposure_boxplot.pdf", pal='terrain',
                                                   threshold=0,plots_per_page=4,...){
  if(!this$.normalized) Normalize(this)
  dp <- dim(this$.Sign) #[i,n,r]
  de <- dim(this$.Exp) #[n,j,r]
  i<-dp[[1]]; n<-dp[[2]]; j<-de[[2]]; r<-de[[3]]
  if(is.null(this$.samples) | all(is.na(this$.samples))){
    samplenames<-paste("Sample",1:j,sep="_")
  }else{
    samplenames<-this$.samples
  }
  # define color palette
  if (pal == 'rainbow') {
    bar.col = rainbow(j, s = 1, v = 1, start = 0, end = (j - 1)/j, alpha = 1)
  } else if (pal == 'heat') {
    bar.col = heat.colors(j, alpha = 1)
  } else if (pal == 'terrain') {
    bar.col = terrain.colors(j, alpha = 1)
  } else if (pal == 'topo') {
    bar.col = topo.colors(j, alpha = 1)
  } else if (pal == 'cm') {
    bar.col = cm.colors(j, alpha = 1)
  } else bar.col = pal
  #Plot
  if(length(grep("\\.pdf$",file))){
    pdf(file)#,width=8,height=2*plots_per_page)
    par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2)#, bty="n")
  }else{
    if(length(grep("\\.png$",file))==0) file<-paste(file,"png",sep=".")
    if(n<=plots_per_page){
      png(file,width=1800,height=360*n)
      par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2)#, bty="n")
    }else{
      png(sub("\\.png$","_%03d.png",file),width=1800,height=360*plots_per_page)
      par(cex= 1, cex.axis=0.5, mfrow=c(plots_per_page,1), mar=c(3,2,2,1), mgp=c(4,0.35,0), xpd=NA, las=2)#, bty="n")
    }
  }
  hg<-rep(1,plots_per_page)
  hg[plots_per_page]<-1.4
  layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
  for(k in 1:n){
    E<- data.frame(t(matrix(as.vector(this$.Exp[k,,]),j,r)))
    y.max <- max(E)
    y.width <- y.max+y.max*.05
    colnames(E)<-samplenames
    medians<-apply(E,2,median)
    q05<-apply(E,2,quantile,0.05)
    q25<-apply(E,2,quantile,0.25)
    q75<-apply(E,2,quantile,0.75)
    q95<-apply(E,2,quantile,0.95)
    medians[medians<threshold]<-0
    q05[q05<threshold]<-0
    q25[q25<threshold]<-0
    q75[q75<threshold]<-0
    q95[q95<threshold]<-0
    bottom <- (k %% plots_per_page==0 | k==n)
    top <- k %% plots_per_page==1
    topmar<-1
    bottommar<-1
    this.names<-NULL
    if(bottom){
      this.names<-samplenames
      bottommar<-5
    }
    if(top){
      #topmar<-2
      lastk <- n-plots_per_page*floor(k/plots_per_page)
      if (lastk<plots_per_page){
        hg[plots_per_page]<-1
        hg[lastk]<-1.4
        layout(mat=matrix(1:plots_per_page,plots_per_page,1),heights=hg)
      }
    }
    par(mar=c(bottommar,3.5,topmar,1),cex=0.7,cex.axis=1)
    if (bottom){
      plot(1:j,medians,ylim=c(0,y.max),xlim=c(0.8,1.05*j+0.5),type="n",main="",xlab="",ylab="",xaxt="n",font=2)
      boxplot(E, col=bar.col, border="black", add=TRUE, cex=0.5, at=1:j,xaxt="n",lwd=0.2,pch=45) 
      axis(side=1,at=1:j,labels=this.names,font=2)
    }else{
      plot(1:j,medians,ylim=c(0,y.max),xlim=c(0.8,1.05*j+0.5),type="n",main="",xlab="",ylab="",xaxt="n",font=2)
      boxplot(E, col=bar.col, border="black", add=TRUE, cex=0.5, at=1:j,xaxt="n",lwd=0.2,pch=45) 
    }
    #Error bars
    mp<-1:j
    segments(x0=mp,y0=q05,x1=mp,y1=q95,col="grey50",lwd=0.2)
    cap <- 1/6
    segments(x0=mp-(cap/3),y0=q05,x1=mp+(cap/3),y1=q05,col="grey40",lwd=0.2)
    segments(x0=mp-cap,y0=q25,x1=mp+cap,y1=q25,col="grey40",lwd=0.2)
    segments(x0=mp-cap,y0=q75,x1=mp+cap,y1=q75,col="grey40",lwd=0.2)
    segments(x0=mp-(cap/3),y0=q95,x1=mp+(cap/3),y1=q95,col="grey40",lwd=0.2)
    #side bars
    rect(j+1,0.01*y.max,1.05*j+1,0.99*y.max,col="grey90",border="grey10")
    text(1.02*j+1,y.max/2,labels=paste("S",k,sep=""),cex=0.9)
  }
  dev.off()
})

