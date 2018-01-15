#' Predict fragment usage
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
.fragment.estimation<-function(output.dir, seq.depth=NULL, chr="chr1", homozygote= FALSE, plot.width=10, plot.height=7){
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
	pred.file=paste(output.dir,"/report.txt",sep="");
  min.depth=seq.depth/4;
  if(homozygote)
    min.depth=seq.depth/2;
	pred=fread(pred.file, header=TRUE,sep="\t");

	locs=c(as.vector(as.matrix(subset(pred[,c("chr.from", "pos.from")],
									chr.from=="insert","pos.from"))),
									as.vector(as.matrix(subset(pred[,c("chr.to", "pos.to")],
									chr.to=="insert","pos.to"))))

	bk=sort(clustercpp(c(50, locs))$center)
	n=length(bk)
	ip=vector();
	for(i in seq_len(n-1)){
		x=GRanges(seqnames=Rle("insert"),
			 			ranges = IRanges(start=bk[i], end=bk[i+1]),
			 			strand=Rle("+"))
		ip[i]=sum(coverage(x))
	}

	input=input[!is.na(input$from),]
	tag=as.vector(input$to)-as.vector(input$from) > 0
	input=input[tag,]
	gr=GRanges(seqnames=Rle(as.vector(input$chr)),
			 ranges = IRanges(start=as.vector(input$from), end=as.vector(input$to)),
			 strand=Rle(rep("+",dim(input)[1])))
	cv=coverage(gr)
	xwin <- sapply(1:(n-1), function(x) sum(window(cv[["insert"]], bk[x], bk[x+1])))
	p=n;
	q=n*(n-1)/2
	mx=matrix(0,ncol=n-1, nrow=q)
	z=1;
	rd1=vector()
	rd2=vector();
	for(i in 1:(n-1)){
		for(j in i:(n-1)){
			for(x in i:j){
				mx[z,x]=1
			}
			rd1[z]=bk[i];
			rd2[z]=bk[j+1];
			z=z+1;
		}
	}
	cc= 2 * xwin/(seq.depth*ip)
  if(homozygote)
     cc = xwin/(seq.depth*ip)
	fit=nnls(t(mx),cc)
	x=round(fit$x)
	ww=which(x>0);
	sm=rep(0,n-1)
	for(i in 1:q){
		sm=sm+x[i]*t(mx)[,i]
	}
	rr=cor(sm,cc)
	print(paste("Consistency of predicted transgenic fragment:", round(rr,digits=3),sep=" "))
	plot.file=paste(output.dir,"/plot_fragment.pdf",sep="");
	pdf(plot.file,width=plot.width, height=plot.height)
	par(mfrow=c(2,1),mar=c(3,4,4,2)+0.1)
	plot(c(1,length(cv[["insert"]])), c(1,sum(x)+2),type="n",xlab="",ylab="transgenic Fragments",yaxt="n",main="Transgene Sequence" )
	tt=2;
	random.col=rainbow(length(ww))
	for(p in 1:length(ww)){
		for(i in 1:x[ww[p]]){
			lines(c(rd1[ww[p]],rd2[ww[p]]), c(tt,tt), lwd=2,col=random.col[p])
			tt=tt+1
		}
	}
	#par(mar=c(2,4,2,2)+0.1)
	#draw.bg()
	par(mar=c(5,4,2,2)+0.1)
	.plotCoverage(cv,"insert")
	dev.off();
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


#draw.bg=function(){
#	fr=c(42,27,27,27,1,27,10,10010,1,27, 1950, 20,42,4, 20, 20, 20)
#	to=c(10217,10217,6228,4328,10217,7664,10010,10217,10109,9010,10217,9976,10217,10217, 3691, 3691, 3691)
#	dsf=data.frame(fr=fr, to=to);
#	od=order(fr);
#	fr=fr[od]
#	to=to[od]
#	gb=GRanges(seqnames=Rle("insert"),
#			 ranges = IRanges(start=fr, end=to),
#			 strand=Rle(rep("+",length(fr))))
#	cb=coverage(gb)
#	dd=as.vector(window(cb[["insert"]], 1, length(cb[["insert"]])))
#	plot(c(1,length(cv[["insert"]])), c(1,length(fr)+2),type="n",xlab="",ylab="transgenes Fragments",yaxt="n",main="" )
#	for(i in 1:length(fr)){
#		lines(c(fr[i],to[i]), c(i+1,i+1), lwd=2,col="green")
#	}
#}

