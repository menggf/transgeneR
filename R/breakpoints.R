#' Predict the breakpoints
#' @import BiocParallel
#' @import BH
#' @import Rcpp
#' @import nnls
#' @import graphics
#' @import grDevices
#' @importFrom stats cor median sd window
#' @importFrom utils head write.table
#' @importFrom data.table fread
#' @importFrom seqinr read.fasta write.fasta
#' @importFrom ShortRead readFastq srduplicated writeFastq
#' @importFrom plyr ldply
#'
.breakpoints<-function(st, output.dir,sites,region,reads,  span=250, col.left="lightblue", col.right="lightgreen",max.show=1000){
	output.dir=sub("/$","",output.dir, perl=T);
	sub.site=subset(sites, sites==st, c("chr.from","pos.from", "chr.to", "pos.to", "direction", "gap","count"))
	direct=as.vector(as.matrix(sub.site[1,"direction"]))
	gap=as.numeric(as.matrix(sub.site[1,"gap"]))
	used.reads=reads[[st]];

	if(direct=="ff" | direct =="rr"){
	  chr1 =as.vector(as.matrix(sub.site[1,"chr.from"]))
	  chr2 =as.vector(as.matrix(sub.site[1,"chr.to"]))
	  pos1.from=as.numeric(as.matrix(sub.site[1,"pos.from"]))-span
	  pos2.from=as.numeric(as.matrix(sub.site[1,"pos.to"]))
    sub.region=region[region$id %in% used.reads & (region$direction=="ff" | region$direction=="rr"),]
	  sub.region=unique(sub.region[order(sub.region$from2),c(-1,-2)])
	  if(direct =="rr"){
	  	chr1 =as.vector(as.matrix(sub.site[1,"chr.to"]))
	  	chr2 =as.vector(as.matrix(sub.site[1,"chr.from"]))
	  	pos1.from=as.numeric(as.matrix(sub.site[1,"pos.to"]))-span
	  	pos2.from=as.numeric(as.matrix(sub.site[1,"pos.from"]))
      sub.region=unique(sub.region[order(sub.region$from3),])
	  }

	  if(gap<0)
	  	    gap=0;
	  dd=as.vector(sub.region$direction);
	  c1=as.vector(sub.region$chr1);
	  f1=as.vector(sub.region$from1);
	  t1=as.vector(sub.region$to1);
	  c2=as.vector(sub.region$chr2);
	  f2=as.vector(sub.region$from2);
	  t2=as.vector(sub.region$to2);
	  c3=as.vector(sub.region$chr3);
	  f3=as.vector(sub.region$from3);
	  t3=as.vector(sub.region$to3);
	  c4=as.vector(sub.region$chr4);
	  f4=as.vector(sub.region$from4);
	  t4=as.vector(sub.region$to4);
	  mm=nrow(sub.region)
	  used=1:mm;
	  if(mm> max.show){
	  	used=sort(sample(used, max.show))
	  	mm=max.show
	  }

	  loc=-1*max(floor(mm/50),1);
	  par(mar=c(4,3,2,1),mgp=c(0,0.5,0))

	  plot(c(1, 2*span + gap), c(1, -2*mm+loc), type="n",xaxt="n",yaxt="n",ylab="Reads",xlab="", axes = 0);
	  text((1+span)/2, 0, chr1,pos=3)
	  text((2*gap+3*span)/2, 0, chr2,pos=3)

	  arrows(0, -1, span, -1, lwd=par("lwd")*2, length=0.1, col=col.left, lend=2)
	  arrows(span+gap, -1, 2*span+gap, -1, lwd=par("lwd")*2, length=0.1, col=col.right, lend=2)
	  if(gap> 0){
	  	#lines(c(span, span+gap),c(-1,-1),lty=3, col="pink")
	  	text(span+gap/2, -1, paste(gap,"bp",sep=""),pos=3, cex=0.8)
	  }
	  pp=1;
	  for(i in used){
	  	if(dd[i] == "ff"){
	  		if(t2[i]-pos1.from < 0 | f2[i]-pos1.from > span)
	  			next();
	  		if(t3[i]-pos2.from < 0 | f3[i]-pos2.from > span)
	  			next();
		  	if(!is.na(f1[i]))
				.mylines(c(f1[i]-pos1.from, t1[i]-pos1.from),c(loc - 2* pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
			   	.mylines(c(f2[i]-pos1.from, t2[i]-pos1.from),c(loc - 2* pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
			   	.mylines(c(f3[i]-pos2.from+span+gap, t3[i]-pos2.from+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f4[i]))
			  	.mylines(c(f4[i]-pos2.from+span+gap, t4[i]-pos2.from+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
			pp=pp+1;
		}
		if(dd[i] == "rr"){
			if(t3[i]-pos1.from < 0 | f3[i]-pos1.from > span)
				next();
			if(t2[i]-pos2.from < 0 | f2[i]-pos2.from > span)
				next();
			if(!is.na(f4[i]))
				.mylines(c(f4[i]-pos1.from, t4[i]-pos1.from),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
			   	.mylines(c(f3[i]-pos1.from, t3[i]-pos1.from),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
			   	.mylines(c(f2[i]-pos2.from+span+gap, t2[i]-pos2.from+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f1[i]))
			  	.mylines(c(f1[i]-pos2.from+span+gap, t1[i]-pos2.from+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
			pp=pp+1;
		}
	  }
       lines(c(span, span), c(loc, -2*mm-2), col="yellow",lty=3)
	     lines(c(gap+span, gap+span), c(loc, -2*mm-2), col="yellow",lty=3)
	     axis(1, at=c(0,span/2, span ),labels=c(pos1.from, pos1.from+0.5*span, pos1.from+span),cex.axis=0.6)
	     axis(3, at=c(span+gap,1.5*span+gap,2*span+gap),labels=c( pos2.from,pos2.from+0.5*span, pos2.from+span),cex.axis=0.6)
   }

   if(direct=="fr"){
	  chr1 =as.vector(as.matrix(sub.site[1,"chr.from"]))
	  chr2 =as.vector(as.matrix(sub.site[1,"chr.to"]))
	  pos1.from=as.numeric(as.matrix(sub.site[1,"pos.from"]))-span
	  pos2.from=as.numeric(as.matrix(sub.site[1,"pos.to"]))
	  sub.region=region[region$id %in% used.reads & (region$direction=="fr" | region$direction=="rf"),]
    sub.region=unique(sub.region[order(sub.region$from2),c(-1,-2)])
	  if(gap<0)
	  	gap=0;
	  dd=as.vector(sub.region$direction);
	  c1=as.vector(sub.region$chr1);
	  f1=as.vector(sub.region$from1);
	  t1=as.vector(sub.region$to1);
	  c2=as.vector(sub.region$chr2);
	  f2=as.vector(sub.region$from2);
	  t2=as.vector(sub.region$to2);
	  c3=as.vector(sub.region$chr3);
	  f3=as.vector(sub.region$from3);
	  t3=as.vector(sub.region$to3);
	  c4=as.vector(sub.region$chr4);
	  f4=as.vector(sub.region$from4);
	  t4=as.vector(sub.region$to4);
	  mm=dim(sub.region)[1]
	  used=1:mm;
	  if(mm> max.show){
	  	used=sort(sample(used, max.show))
	  	mm=max.show
	  }
	  loc=-1*max(floor(mm/50),2);
	  par(mar=c(4,3,2,1),mgp=c(0,0.5,0))
	  plot(c(1, 2*span + gap), c(1, -2*mm+loc), type="n",xaxt="n",yaxt="n",ylab="Reads",xlab="", axes = 0);
	  text((1+span)/2, 0, chr1,pos=3)
	  text((2*gap+3*span)/2, 0, chr2,pos=3)
	  arrows(0, -1, span, -1, lwd=par("lwd")*2, length=0.1, col=col.left, lend=2)
	  arrows(2*span+gap, -1, span+gap, -1, lwd=par("lwd")*2, length=0.1, col=col.right, lend=2, code=2)
	  if(gap> 0){
	  	#lines(c(span, span+gap),c(-1,-1),lty=3, col="pink")
	  	text(span+gap/2, -1, paste(gap,"bp",sep=""),pos=3, cex=0.8)
	  }
	  pp=1;
	  for(i in used){
	  	if(direct == dd[i]){
	  		if(t2[i]-pos1.from < 0 | f2[i]-pos1.from > span)
	  			next()
	  		if(pos2.from-f3[i] < 0 | pos2.from-t3[i] > span)
	  			next()
	  		if(!is.na(f1[i]))
		    	.mylines(c(f1[i]-pos1.from, t1[i]-pos1.from),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
		 	  	.mylines(c(f2[i]-pos1.from, t2[i]-pos1.from),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
		   		.mylines(c(pos2.from-t3[i]+span+gap,pos2.from-f3[i]+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f4[i]))
		  		.mylines(c(pos2.from-t4[i]+span+gap,pos2.from-f4[i]+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
		  	pp=pp+1
		}
		if(direct != dd[i]){
			if(t3[i]-pos1.from < 0 | f3[i]-pos1.from > span )
				next();
			if(pos2.from-f2[i] < 0 | pos2.from-t2[i] >  span )
				next();
			if(!is.na(f4[i]))
		    	.mylines(c(f4[i]-pos1.from, t4[i]-pos1.from),c(loc - 2*pp, loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
		 	  	.mylines(c(f3[i]-pos1.from, t3[i]-pos1.from),c(loc - 2*pp, loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
		   		.mylines(c(pos2.from-t2[i]+span+gap, pos2.from-f2[i]+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f1[i]))
		  		.mylines(c(pos2.from-t1[i]+span+gap,pos2.from-f1[i]+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
		  	pp=pp+1
		}
	  }
	  lines(c(span, span), c(loc, -2*mm-2), col="yellow",lty=3)
	  lines(c(gap+span, gap+span), c(loc, -2*mm-2), col="yellow",lty=3)
	  axis(1, at=c(0,span/2, span ),labels=c(pos1.from, pos1.from+0.5*span, pos1.from+span),cex.axis=0.6)
	  axis(3, at=c(span+gap,1.5*span+gap,2*span+gap),labels=c( pos2.from,pos2.from+0.5*span, pos2.from+span),cex.axis=0.6)
   }
   if(direct=="rf"){
	  chr1 =as.vector(as.matrix(sub.site[1,"chr.from"]))
	  chr2 =as.vector(as.matrix(sub.site[1,"chr.to"]))
	  pos1.from=as.numeric(as.matrix(sub.site[1,"pos.from"]))+span
	  pos2.from=as.numeric(as.matrix(sub.site[1,"pos.to"]))
	  sub.region=region[region$id %in% used.reads & (region$direction=="fr" | region$direction=="rf"),]
    sub.region=unique(sub.region[order(sub.region$from2),c(-1,-2)])
	  if(gap<0)
	  	gap=0;
	  dd=as.vector(sub.region$direction);
	  c1=as.vector(sub.region$chr1);
	  f1=as.vector(sub.region$from1);
	  t1=as.vector(sub.region$to1);
	  c2=as.vector(sub.region$chr2);
	  f2=as.vector(sub.region$from2);
	  t2=as.vector(sub.region$to2);
	  c3=as.vector(sub.region$chr3);
	  f3=as.vector(sub.region$from3);
	  t3=as.vector(sub.region$to3);
	  c4=as.vector(sub.region$chr4);
	  f4=as.vector(sub.region$from4);
	  t4=as.vector(sub.region$to4);
	  mm=dim(sub.region)[1]
	  used=1:mm;
	  if(mm> max.show){
	  	used=sort(sample(used, max.show))
	  	mm=max.show
	  }
	  loc=-1*max(floor(mm/50),2);
	  par(mar=c(4,3,2,1),mgp=c(0,0.5,0))
	  plot(c(1, 2*span + gap), c(1, -2*mm+loc), type="n",xaxt="n",yaxt="n",ylab="Reads",xlab="",axes = 0);
	  text((1+span)/2, 0, chr1,pos=3)
	  text((2*gap+3*span)/2, 0, chr2,pos=3)
	  arrows(span, -1, 0, -1, lwd=par("lwd")*2, length=0.1, col=col.left, lend=2)
	  arrows(span+gap, -1, 2*span+gap, -1, lwd=par("lwd")*2, length=0.1, col=col.right, lend=2, code=2)
	  if(gap> 0){
	  	#lines(c(span, span+gap),c(-1,-1),lty=3, col="pink")
	  	text(span+gap/2, -1, paste(gap,"bp",sep=""),pos=3, cex=0.8)
	  }
	  pp=1;
	  for(i in used){
	  	if(direct == dd[i]){
	  		if(pos1.from-f2[i] <0 | pos1.from-t2[i] > span)
	  			next();
	  		if(t3[i]-pos2.from <0 | f3[i]-pos2.from > span)
	  			next();
	  		if(!is.na(f1[i]))
		    	.mylines(c(pos1.from-t1[i], pos1.from-f1[i]),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
		 	  	.mylines(c(pos1.from-t2[i], pos1.from-f2[i]),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
			   	.mylines(c(f3[i]-pos2.from + span+gap, t3[i]-pos2.from+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f4[i]))
			  	.mylines(c(f4[i]-pos2.from + span+gap, t4[i]-pos2.from+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
			 pp=pp+1
		}
		if(direct != dd[i]){
			if(pos1.from-f3[i] < 0 | pos1.from-t3[i] > span)
				next();
			if(t2[i]-pos2.from < 0 | f2[i]-pos2.from > span)
				next();
			if(!is.na(f4[i]))
		    	.mylines(c(pos1.from-t4[i], pos1.from-f4[i]),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f3[i]))
		 	  	.mylines(c(pos1.from-t3[i], pos1.from-f3[i]),c(loc - 2*pp,loc-1 - 2*pp),col=col.left)
			if(!is.na(f2[i]))
			   	.mylines(c(f2[i]-pos2.from + span+gap, t2[i]-pos2.from+span+gap),c(loc - 2*pp, loc-1 - 2*pp),col=col.right)
			if(!is.na(f1[i]))
			  	.mylines(c(f1[i]-pos2.from + span+gap, t1[i]-pos2.from+span+gap),c(loc - 2*pp,loc-1 - 2*pp),col=col.right)
			pp=pp+1
		}
	  }
	  lines(c(span, span), c(loc, -2*mm-2), col="yellow",lty=3)
	  lines(c(gap+span, gap+span), c(loc, -2*mm-2), col="yellow",lty=3)
	  axis(1, at=c(0,span/2, span ),labels=c(pos1.from, pos1.from+0.5*span, pos1.from+span),cex.axis=0.6)
	  axis(3, at=c(span+gap,1.5*span+gap,2*span+gap),labels=c( pos2.from,pos2.from+0.5*span, pos2.from+span),cex.axis=0.6)
   }

}

.mylines<-function(x,y,col){
	rect(x[1],y[2],x[2],y[1], border=col)
}
