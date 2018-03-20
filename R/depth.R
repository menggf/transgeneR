#' Predict fragment usage
#'
#' Predict fragment usage
#'
#' @param output.dir a directory for all the output. It is better to be a empty or non-existing directory;
#' @param min.counts the minimum couts cutoff
#' @param seq.depth the sequencing depth
#' @param chr the chromosome to estimate the sequencing depth
#' @param homozygote homozygote or heterozygote?
#' @param plot.width plot width, equal to "width" in pdf()
#' @param plot.height plot height, equal to "height" in pdf()
#'
#' @author Guofeng Meng
#' @references To be updated later
#'
#' @import BiocParallel
#' @import BH
#' @import Rcpp
#' @import nnls
#' @import graphics
#' @import grDevices
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @importFrom stats cor median sd window
#' @importFrom utils head write.table
#' @importFrom data.table fread
#' @importFrom seqinr read.fasta write.fasta
#' @importFrom ShortRead readFastq srduplicated writeFastq
#' @importFrom plyr ldply
#'
#'
#' @return a plot in pdf file: plot_fragment.pdf
#'
#' @examples
#' test.data=system.file("extdata", "test.zip", package = "transgeneR")
#' temp <- tempdir()
#' unzip(test.data, exdir=temp)
#' output.dir=paste(temp,"/test",sep="")
#' fragment.estimation(output.dir)
#' @export

fragment.estimation<-function(output.dir, min.counts=0, seq.depth=NULL,  chr="chr1", homozygote= FALSE, plot.width=10, plot.height=7){
	if(is.null(seq.depth)){
		chr.file=paste(output.dir,"/temp_files/read_nonsplit_", chr,".txt",sep="");
		leftchr=fread(chr.file,sep="\t", showProgress=FALSE);
		x=GRanges(seqnames=Rle("chr"),
			 	ranges = IRanges(start=as.vector(leftchr$V2), end=as.vector(leftchr$V3)),
			 	strand=Rle("+"))
		rm(leftchr);
		chrcv=coverage(x);
		rm(x);
		chrwin <- as.vector(window(chrcv[["chr"]], 1, length(chrcv[["chr"]])))
		rm(chrcv)
		chrwin=chrwin[chrwin!=0]
		nn=length(chrwin)
		seq.depth=round(median(chrwin[round(nn/4):round(nn*3/4)], na.rm=TRUE))
		print(paste("Estimated sequence depth is:",round(seq.depth, digits=1) ,sep=" "))
		leftchr1=1;leftchr=1;x=1;chrwin=1;
	}
	
	insert.file=paste(output.dir,"/temp_files/insert.fa",sep="");
    sq=read.fasta(insert.file, as.string = TRUE, forceDNAtolower = FALSE, seqonly = TRUE)[[1]];
    n.seq=nchar(sq)
	#left.file=paste(output.dir,"/temp_files/read_left_insert.txt",sep="");
	range.file=paste(output.dir,"/temp_files/range_split.txt",sep="");
	both.file=paste(output.dir,"/temp_files/read_both.txt",sep="");
	nonsplit.file=paste(output.dir,"/temp_files/read_nonsplit_insert.txt",sep="");
	split.file=paste(output.dir,"/temp_files/read_split.txt",sep="");
	#left=fread(left.file, header=FALSE,sep="\t")
	both=fread(both.file, header=TRUE,sep="\t")
	rang=fread(range.file, header=TRUE,sep="\t")
	split=fread(split.file, header=TRUE,sep="\t")
	nonsplit=fread(nonsplit.file, header=F,sep="\t")
	names(nonsplit) <- c("tag","from","to")
    nonsplit$chr<- "insert"
	#names(left)=c("tag","id","chr","from","to")
	temp=both[,c("chr1","from","to")]
	names(temp)<-c("chr","from","to")
	temp1=split[,c("chr1","from1","to1")];
	temp2=split[,c("chr2","from2","to2")];
	temp3=split[,c("chr3","from3","to3")];
	temp4=split[,c("chr4","from4","to4")];
	names(temp1)<-c("chr","from","to")
	names(temp2)<-c("chr","from","to")
	names(temp3)<-c("chr","from","to")
	names(temp4)<-c("chr","from","to")
          
	tmp1=rang[,c("chr1","from1","to1")]
	tmp2=rang[,c("chr2","from2","to2")];
	names(tmp1)<-c("chr","from","to")
	names(tmp2)<-c("chr","from","to")
	input=unique(rbind( nonsplit[,c("chr","from","to")], temp1,temp2,temp3,temp4,temp,tmp1,tmp2))
	rm(temp1, temp2, temp3, temp4, tmp1, tmp2, both, split, nonsplit );
	new.input=input[input$to > input$from & input$chr=="insert",]
	x1=GRanges(seqnames=Rle("insert"), 
  		ranges = IRanges(start=as.vector(new.input$from), end=as.vector(new.input$to)),
  		strand=Rle("+"));
  	cv1=coverage(x1);
	pred.file=paste(output.dir,"/report.txt",sep="");
    min.depth=seq.depth/4;
    if(homozygote)
        min.depth=seq.depth/2;
	pred=subset(fread(pred.file, header=TRUE,sep="\t"), count > min.counts);
	
    sub.sites1=subset(pred, chr.from != "insert" & chr.to == "insert");
      sub.sites2=subset(pred, chr.from == "insert");
      
      from1.chr=as.vector(sub.sites1$chr.from);
      from1.pos=as.vector(sub.sites1$pos.to);
      names(from1.chr)<-from1.pos    
      
      dd1=as.vector(sub.sites1$direction)
      dd2=as.vector(sub.sites2$direction)
      p1=as.vector(sub.sites1$pos.to)
      p3=as.vector(sub.sites2$pos.from)
      p4=as.vector(sub.sites2$pos.to)
      from1=vector();
      to1=vector();
      for(i in seq_along(dd1)){
         if(dd1[i]=="ff" | dd1[i]=="rf")
           from1=append(from1, p1[i]);
         if(dd1[i]=="rr" | dd1[i]=="fr")
           to1=append(to1, p1[i]);
      }
      from1=sort(unique(from1));
      to1=sort(unique(to1));
      from2=vector();
      to2=vector();
      for(i in seq_along(dd2)){
         if(dd2[i]=="ff" ){
           from2=append(from2, p4[i]);
           to2=append(to2, p3[i]);
         }
         if(dd2[i]=="rr" ){
           from2=append(from2, p3[i]);
           to2=append(to2, p4[i]);
         }
         if( dd2[i]=="fr"){
           to2=append(to2, p3[i]);
           to2=append(to2, p4[i]);
         }
         if( dd2[i]=="rf"){
           from2=append(from2, p3[i]);
           from2=append(from2, p4[i]);
         }
      }
      from2=sort(from2)
      to2=sort(to2)
      from=sort(c(from1, from2))
      to=sort(c(to1, to2))
      cmb.from=vector();
      cmb.to=vector();
      cmb.dd=vector()
      lab1=vector()
      lab2=vector()
      lab =vector()
      ann=matrix(0, ncol=length(to),nrow=length(from))
      for(i in seq_along(from)){
          for(j in seq_along(to)){
            if(from[i] < to[j]){
              cmb.from=append(cmb.from, from[i])
              cmb.to=append(cmb.to, to[j])
              cmb.dd=append(cmb.dd, 1)
              lab1=append(lab1, i)
              lab2=append(lab2, j)
              lab=append(lab, paste(i,"/",j,sep=""))
              ann[i,j]=1
            }
          }
      }
  
    win=200;
    step=20
    w1=seq(1, n.seq-win, by=step)
    w2=w1+win
    mm=length(cmb.from)
    nn=length(w1)
    mx=matrix(nrow=mm, ncol=nn);
    label=vector();
    for(i in 1:mm){
      for(j in 1:nn){
        mx[i,j]=length(intersect(cmb.from[i]:cmb.to[i], w1[j]:(w2[j]-1))) *cmb.dd[i]
      }
    }
    
    
    cvs=vector()
    for(j in 1:nn){
      cvs=append(cvs,  sum(as.vector(window(cv1[["insert"]], w1[j], w2[j]))));
    }
    
    fit=nnls(t(mx),cvs/seq.depth);
    cof=fit$x
    uniq.lab1=unique(lab1)
    have=vector()
    used=vector();
    for(x in uniq.lab1){
        wh=lab1==x
        test=cof;
        test[!wh]=0;
        max.cof=which.max(test);
        if(cof[max.cof] > 0.1){
          have=append(have,max.cof);
          used=append(used, x)
        }
    }
    mx2=matrix(nrow=length(have), ncol=nn);
    for(i in 1:length(have) ){
      for(j in 1:nn){
        mx2[i,j]=length(intersect(cmb.from[have[i]]:cmb.to[have[i]], w1[j]:(w2[j]-1))) *cmb.dd[have[i]]
      }
    }
    est.cvs=colSums(mx2);
    myr=cor(est.cvs, cvs)
    left=uniq.lab1[!uniq.lab1 %in% used];
    final.rr=myr
    while(length(left) >0){
      which.left=which(lab1%in%left);
      rr=sapply(which.left, function(y){
        ad=sapply(1:nn, function(z) length(intersect(cmb.from[y]:cmb.to[y], w1[z]:(w2[z]-1))) *cmb.dd[y])
        cor(cvs, ad + est.cvs)
      })
      wh.rr=which.max(rr);
      if(rr[wh.rr] < final.rr){
        break;
      }
      have=append(have, which.left[wh.rr])
      left=left[left!=lab1[which.left[wh.rr]]]
      final.rr=rr[wh.rr]
      #print(final.rr)
      mx3=matrix(nrow=length(have), ncol=nn);
      label=vector();
      for(i in 1:length(have)){
        for(j in 1:nn){
          mx3[i,j]=length(intersect(cmb.from[have[i]]:cmb.to[have[i]], w1[j]:(w2[j]-1))) *cmb.dd[have[i]]
        }
      }
      est.cvs=colSums(mx3);
      myr2=cor(est.cvs, cvs)
      #print(myr2)
    }
    df=data.frame(from=cmb.from, to=cmb.to, dd=cmb.dd, cof=fit$x)
     sub.df=df[have,]
      
    labs=vector()
    
    for(j in 1:nn){
      labs[j]=(w1[j]+w2[j])/2
    }
    est.cvs=colSums(mx3);
    pdf(paste(output.dir,"/plot_fragment.pdf",sep=""), width=10,height=14)
      par(mfrow=c(3,1))
      .plotCoverage(cv1,"insert")
      
      lb=seq(0, n.seq, by=300);
      axis(3, at=seq(0, n.seq, by=300),label=lb,cex.axis=0.6)
      tt=mean(cvs)/mean(est.cvs)
      plot(labs, cvs/tt, col="red",type="l",xlab="transgene", ylab="Depth",lwd=4)
      lines(labs, est.cvs,col="blue",lty=3,lwd=4)
      legend("topleft",c("Sequencing depth","Predicted depth"), lwd=4,lty=c(1,3), col=c("red","blue"))
      
      mm=length(have)
      loc=-1 * max(floor(mm/50),1);
      input=data.frame(from=cmb.from[have], to=cmb.to[have])
      input=input[order(input$from),]
      ann.chr1=from1.chr[as.character(input$from)]
      ann.chr2=from1.chr[as.character(input$to)]
      plot(c(1, n.seq), c(1, -2*length(have)+loc), type="n",xaxt="n",yaxt="n",ylab="Fragments",xlab="",axes = 0);
      axis(3, at=seq(0, n.seq, by=300),label=seq(0, n.seq, by=300),cex.axis=0.6)
      for(i in 1:mm){
        rect(input[i,1],-2*i+1,input[i,2],-2*i+0.5, border="green")
        if(!is.na(ann.chr1[i])){
          text(input[i,1]-30, -2*i+1, ann.chr1[i], cex=0.6)
        }
        if(!is.na(ann.chr2[i])){
          text(input[i,2]+30, -2*i+1, ann.chr1[i], cex=0.6)
        }
      }
  
  dev.off()
  
}    
   
	
.plotCoverage <-function(x, chrom, start=1, end=length(x[[chrom]]), col="blue", xlab="Index", ylab="Coverage"){
	xWindow <- as.vector(window(x[[chrom]], start, end))
	x <- start:end
	xlim <- c(start, end)
	ylim <- c(min(xWindow), min(max(xWindow),50000))
	plot(x = start, y = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "", type = "n")
	polygon(c(start, x, end), c(0, xWindow, 0), col = col)
}

.clustervector<-function(x, rg, min.split.reads=1){
	cl=rep(0,length(x));
	zz=0;
	center=vector();
	while(length(cl[cl==0])>0){
		zz=zz+1;
		tx=table(x[cl==0]);
		cc=as.numeric(names(which.max(tx)))
		center=append(center, cc);
		tag=cl==0 & x > cc-rg & x < cc+rg;
    if(length(tag[tag]) < min.split.reads)
      break;
		cl[tag]=zz;
	}
	names(center)=1:zz
	return(list(cl=cl, center=center))
}


