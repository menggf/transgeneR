#' Predict fragment usage
#' @import BiocParallel
#' @import BH
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
.predict.recombination.sites<-function(input, ranges, splts, region=10,min.split.reads=2, cores=1){
	labels=paste(as.vector(input$chr1),as.vector(input$chr2), as.vector(input$direction),sep="#")
	label=table(labels);
	splts.ids = as.vector(splts$id)
	splts.c.from = as.vector(splts$chr2)
	splts.c.to = as.vector(splts$chr3)
	res=bplapply(names(sort(label[label >= min.split.reads ], decreasing=T)), function(x){
		wh=labels==x;
		sub.input=input[wh,];
		nn=length(wh[wh]);
		if(nn < min.split.reads){
			return(NULL)
		}
		cl1 =clustercpp(c(region, as.numeric(sub.input$from)))
		cl2 =clustercpp(c(region, as.numeric(sub.input$to)))
		cl=paste(cl1$cl,"/", cl2$cl,sep="")
		chr.from=vector();
		pos.from=vector();
		chr.to=vector();
		pos.to=vector();
		types=vector();
		dirs=vector();
		gap=vector();
		gap.sd=vector();
		count=vector()
		score=vector();
		reads=vector();
		bal=vector();
     bal2=vector();
		cross=vector();
		for(cc in unique(cl)){
			wtag=cl==cc;
       ccs=strsplit(cc, "/")[[1]];
       if(ccs[1]=="0" | ccs[2]=="0")
         next();
			mm=length(which(wtag));
			if(mm < min.split.reads)
				next();
			sub.sub.input=unique(sub.input[wtag,])
			ids=as.vector(sub.sub.input$id)

			c.from=as.vector(sub.sub.input$chr1)[1];
			c.to=as.vector(sub.sub.input$chr2)[1];
			sub.splts=splts[splts.ids %in% ids & splts.c.from == c.from & splts.c.to == c.to,]
			xx=nrow(sub.splts)
			if(xx < min.split.reads)
			    next();
			dirs=append(dirs,as.vector(sub.sub.input$direction)[1]);
			chr.from=append(chr.from, c.from);
			chr.to=append(chr.to, c.to);
			left.dd =apply(sub.splts[, c(5,6,8,9)],    1, function(y) max(y, na.rm = TRUE)- min(y[y!=0], na.rm = TRUE))
			right.dd=apply(sub.splts[,c(11,12,14,15)], 1, function(y) max(y, na.rm = TRUE)- min(y[y!=0], na.rm = TRUE))
			count=append(count, xx)
			p.from = cl1$center[as.numeric(ccs[1])];
			p.to =   cl2$center[as.numeric(ccs[2])];
			dd = min( max(left.dd), max(right.dd ), na.rm =TRUE);
			bal = append(bal, dd)

      bal2 = append(bal2, length(unique(sub.sub.input$assign)))
			pos.from = append(pos.from, p.from)
			pos.to = append(pos.to,  p.to)

			if(c.from == c.to){
			  tg= ranges$chr1 == c.from & ranges$chr2==c.to;
				sub.range=ranges[tg,]
				tg=((abs(sub.range$from1 - p.from) < 50 | abs(sub.range$to1 - p.from) < 50)
				    & (abs(sub.range$from2 - p.to) < 50 | abs(sub.range$to2 - p.to) < 50)) |
				  ((abs(sub.range$from2 - p.from) < 50 | abs(sub.range$to2 - p.from) < 50)
				   & (abs(sub.range$from1 - p.to) < 50 | abs(sub.range$to1 - p.to) < 50))
				sub.sub.range=sub.range[tg,]
				cross=append(cross, dim(sub.sub.range)[1]);
			}
			if(c.from != c.to){
        tg= ranges$chr1==c.from & ranges$chr2==c.to
				sub.range1=ranges[tg,]
				tg=((abs(sub.range1$from1 - p.from) < 50 | abs(sub.range1$to1 - p.from) <50)
				  & (abs(sub.range1$from2 - p.to) < 50 | abs(sub.range1$to2 - p.to) <50))
				sub.sub.range1=sub.range1[tg,]
				tg=ranges$chr2==c.from & ranges$chr1==c.to
				sub.range2=ranges[tg,]
				tg=((abs(sub.range2$from2 - p.from) < 50 | abs(sub.range2$to2 - p.from) <50)
				  & (abs(sub.range2$from1 - p.to) < 50 | abs(sub.range2$to1 - p.to) <50))
				sub.sub.range2=sub.range2[tg,]
				cross=append(cross, dim(sub.sub.range1)[1] + dim(sub.sub.range2)[1]);
			}
			score=append(score, mean(as.numeric(sub.sub.input$score)))
			gap=append(gap, median(as.numeric(sub.sub.input$gap)))
			reads=append(reads, paste(as.vector(sub.sub.input$id), collapse=","))
			gap.sd=append(gap.sd, sd(as.numeric(sub.sub.input$gap)))
		}
		if(length(chr.from)==0)
			 return(NULL)
		resu=data.frame(chr.from=chr.from, pos.from=pos.from, chr.to=chr.to,
                     pos.to=pos.to, direction=dirs, gap=gap,
                     gap.sd= round(gap.sd,digits=1), count=count,
                     cross=cross, score=round(score,digits=0),balance.len=bal, balance.assign=bal2,
                     reads=reads)
		return(resu)
	}, BPPARAM= MulticoreParam( workers= cores))
	res[sapply(res, is.null)] <- NULL
	dif <- ldply(res, data.frame)
  if(nrow(dif) > 1)
	dif=dif[order(dif$count, decreasing=T),]
	mm=dim(dif)[1];
	sites=paste("site",1:mm, sep="")
	output=cbind(sites,dif)
	return(output)
}

#' Refine the sites
.refine.sites <-function(output.dir, sites, min.clip=30, cores=1){
  chr.from = as.vector(sites$chr.from)
  chr.to   = as.vector(sites$chr.to)
  pos.from = as.vector(sites$pos.from)
  pos.to   = as.vector(sites$pos.to)
  used.chr=unique(c(chr.from, chr.to))

  dd=min(min.clip, 40);
  nonsplit.left=rep(0,length(chr.from))
  nonsplit.right=rep(0,length(chr.from))
  for(chr in used.chr){
    #print(chr)
    rin=fread(paste(output.dir,"/temp_files/read_nonsplit_",chr,".txt",sep="" ),header=FALSE,sep="\t", showProgress=FALSE)
    rin2=rin[order(rin$V2),]
    from = as.vector(rin2$V2) + dd
    to =as.vector(rin2$V3) - dd
    wh.from=chr.from==chr
    wh.to=chr.to==chr
    #print("ok")
    nonsplit.left[wh.from]=cmpitcpp(c(cores, length(wh.from[wh.from]), length(from), pos.from[wh.from], from,to) )
    nonsplit.right[wh.to]= cmpitcpp(c(cores, length(wh.to[wh.to]),     length(from), pos.to[wh.to],   from,to) )
  }
  temp=data.frame(ns.left=nonsplit.left,ns.right=nonsplit.right)
  return(cbind(sites, temp)[,c("sites","chr.from","pos.from", "chr.to", "pos.to", "direction", "gap", "gap.sd", "count", "cross", "ns.left", "ns.right","score", "balance.len", "balance.assign" ,"reads")])
}

#' cluster the breakpoints
.clustervector<-function(x, rg, min.split.reads=1){
	cl=rep(0,length(x));
	zz=0;
	center=vector();
	while(length(cl[cl==0])>0){
		zz=zz+1;
		tx=table(x[cl==0]);
		cc=as.numeric(names(which.max(tx)))
		center=append(center, cc);
		tag=cl==0 & x >= cc-rg & x <= cc+rg;
    if(length(tag[tag]) < min.split.reads)
      break;
		cl[tag]=zz;
	}
	names(center)=1:zz
	return(list(cl=cl, center=center))
}
