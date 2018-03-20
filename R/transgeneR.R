#' A integrated tool for transgenic integration discovery using sequencing data
#'
#' A integrated tool for transgenic integration discovery using sequencing data
#'
#' @param output.dir a directory for all the output. It is better to be a empty or non-existing directory;
#' @param genome.ref the bowtie2 genome reference prefix. It should be "bt2_index_base" of the bowtie2-build;
#' @param homozygote the transgene sequence is integrated as homozygote? Default: False;
#' @param mate1 the first fastq file of pair-end reads;
#' @param mate2 the second fastq file of pair-end reads;
#' @param insert.seq the transgenic sequence in fasta format;
#' @param estimate.copy to estimate copy number or not. This option can be set TRUE only when whole genome sequencing data is available
#' @param calculate.nonsplit.reads to estimate number of nonsplit reads crossing the integration or rearrangement sites
#' @param seq.depth the sequencing depth of the whole genome sequencing data. It is used to estimated copy number of transgenic sequences. This option is used only when estimate.copy is TRUE
#' @param backgroud.chrom the chromosome used to estimate the sequencing information, e.g. coverage and copy number
#' @param remove.duplicate to remove the duplicated reads or not. Recommend to set it TURE if PCR-based sequencing data are used
#' @param min.clip the minimum clip size for reads split
#' @param max.tail the maximum tail unmapped length
#' @param min.split.reads the minimum number of reads for a validated integration site
#' @param cores the thread number
#' @param max.fragment equal to the setting of "-X" of bowtie2
#' @param allow.gap the maximum gap between two part of insertion
#' @param homo.windows the windows size during homologous sequence searching
#' @param homo.step the steps during homologous sequence searching
#' @param strict bool value, to use discordanantly mapped reads or not
#' @param region the tolarated region size of integration sites
#' @param span the spanned size to plot the integration sites
#' @param col.left the color of left part of integration site plot
#' @param col.right the color of right part of integration site plot
#' @param max.show the maximum number of read to display in integration site plot
#' @param gz bool value, to compressed the alignment output or not. Not recommend to set it TRUE
#' @param block.size the number of reads to analysis each run; it is used to control ROM usage
#' @param reads.dup.score the minimum alignment score difference for unique mapping
#'
#' @useDynLib transgeneR
#'
#' @author Guofeng Meng
#' @references To be updated later
#'
#' @import BiocParallel
#' @import BH
#' @import Rcpp
#' @import GenomicRanges
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
#'
#' @details This a one-stop analysis pipeline to find the transgenic intergation sites and transgenic recombination. When the whole genome sequencing data are available, it can estimate the copy number of transgenic sequences in animal genome.  To do the whole analysis, it has following steps:
#' \itemize{
#' \item Build the bowtie2 reference for transgenic sequence. The output is store in a directory "insert_ref/";
#' \item Map the homologous regions of transgenic sequences in genome. Output is a file "homo.txt";
#' \item reads local alignment to both genome and transgenic sequences: in this step, users need to pre-install bowtie2 and build the genome reference of studied animals. This step will generate two file: aln_genome.sam and aln_insert.sam.
#' \item Assign the reads to genome, transgenic sequence or both and collect the clipping parts of read for second-round alignments. The output will be store in a dictory "temp_files/";
#' \item Second-round alignments. The output are stored as "temp_files/fragment_genome.sam" and "temp_files/fragment_insert.sam";
#' \item Connect the break sites in either transgenic sequence,  genome or between transgenic sequence and genome;
#' \item Make the plot for the split reads  in the integration sites. Figures are created in "sites/*.pdf";
#' \item If whole genome sequencing data are used, it can estimated both the full and imcomplete integration; and calculate their copy numbers. One figure is drawn as "plot_fragment.pdf"
#' }
#'
#' In same case that users wish to re-run part of the analysis, user can just delete the output in mentioned steps and this will make it to skip the steps with outputs.
#'
#'
#'
#' @return The output are files located in "output.dir". They have a structure of:
#' output.dir:/
#' \itemize{
#'         \item aln_genome.sam (or aln_genome.sam.gz)
#'         \item aln_insert.sam (or aln_insert.sam.gz)
#'         \item assign.txt : read assignment in genome and transgenic sequence
#'         \item homo.txt : homologous annotation of transgeneic sequence
#'         \item mapping_summary.txt : reads alignment information
#'         \item report.txt : the predicted results
#'         \item warning.txt : the warning information
#'         \item plot_fragment.pdf the copy information of transgenic sequence
#'         \item insert_ref/ : the bowtie2 reference for transgenic sequence
#'         \item temp_files/: some temparory files
#'         \item sites/: the transgenic integration or recombination information
#'        \enumerate{
#'             \item site1.pdf: the plot for first site
#'             \item site2.pdf: the plot for second site
#'      }
#' }
#'
#' @examples
#' #\dontrun{
#' test.data=system.file("extdata", "test.zip", package = "transgeneR")
#' temp <- tempdir()
#' unzip(test.data, exdir=temp)
#' output.dir=paste(temp,"/test",sep="")
#' transgene.sq=paste(temp,"/test/temp_files/insert.fa",sep="")
#' transgeneR(output.dir,  insert.seq=transgene.sq, estimate.copy=FALSE)
#' #}
#' @export

transgeneR<-function(output.dir, genome.ref, homozygote= FALSE, mate1=NULL, mate2=NULL, insert.seq=NULL, estimate.copy=FALSE,calculate.nonsplit.reads=TRUE, seq.depth=NULL, backgroud.chrom="chr10", remove.duplicate=FALSE, min.clip=30, max.tail=10, min.split.reads=2, cores=1, max.fragment=600, allow.gap=100, homo.windows=100, homo.step=50, strict=1, region=10, span=250, col.left="lightblue", col.right="lightgreen",max.show=1000, gz=F, block.size=2000000, reads.dup.score=50){
    # judge bowtie2 installed or not
    if(Sys.info()['sysname']=="Linux" | Sys.info()['sysname']=="Darwin"){
        info=system("which bowtie2",intern=TRUE);
        if(length(grep("no bowtie2", info))!=0)
            print("Warning: bowtie2 should be installed or exported in $PATH")
    }
  if(Sys.info()['sysname']=="Windows"){
    info=system("where bowtie2",intern=TRUE, ignore.stdout=TRUE,ignore.stderr=TRUE,show.output.on.console =FALSE);
    if(length(grep("Could not find", info))!=0)
      print("Warning: bowtie2 should be installed or exported in $PATH")
  }
    # creat directory
    output.dir=sub("/$","",output.dir, perl=T);
    ref.dir=paste(output.dir,"/insert_ref",sep="");
    temp.dir=paste(output.dir,"/temp_files",sep="");
    if(!file.exists(output.dir)){
        dir.create(output.dir);
        dir.create(ref.dir);
        dir.create(paste(output.dir,"/sites",sep=""));

    }
    if(!file.exists(ref.dir))
        dir.create(ref.dir);
    if(!file.exists(temp.dir))
        dir.create(temp.dir);

    #remove duplate reads

    if(remove.duplicate ){
        if(!file.exists(paste(output.dir,"/temp_files/input_R1.fastq.gz",sep=""))){
            fastq1=readFastq(mate1)
            fastq2=readFastq(mate2)
            dup1=srduplicated(fastq1)
            dup2=srduplicated(fastq2)
            dup=dup1 & dup2;
            mate1=paste(output.dir,"/temp_files/input_R1.fastq.gz",sep="")
            mate2=paste(output.dir,"/temp_files/input_R2.fastq.gz",sep="")
            writeFastq(fastq1[!dup], mate1, compress=TRUE,mode="w")
            writeFastq(fastq2[!dup], mate2, compress=TRUE, mode="w")
            rm(fastq1)
            rm(fastq2)
        }    else{
            mate1=paste(output.dir,"/temp_files/input_R1.fastq.gz",sep="")
            mate2=paste(output.dir,"/temp_files/input_R2.fastq.gz",sep="")
        }
    }

    # re-write insert file
    insert.file=paste(output.dir,"/temp_files/insert.fa",sep="");
    if(!file.exists(insert.file)){
        if(is.null(insert.seq))
            stop("Error: insert.seq should be given as fasta file");
        write.fasta(read.fasta(insert.seq, as.string = TRUE, forceDNAtolower = FALSE, seqonly = TRUE)[[1]], "insert", file.out=insert.file);
    }

    #building the insert refrerence
    insert.ref=paste(output.dir,"/insert_ref/insert",sep="");
    if(!file.exists(paste(insert.ref,".1.bt2",sep=""))){
        cmd=paste("bowtie2-build --threads", cores, "-q", insert.file, insert.ref, sep=" ");
        system(cmd);
    }

    # homo searching
    homo.file=paste(output.dir,"/homo.txt",sep="");
    if(!file.exists(homo.file)){
        dna.seq=read.fasta(insert.file, as.string = TRUE, forceDNAtolower = FALSE, seqonly = TRUE)[[1]]
        from=1;
        num=nchar(dna.seq)
        temp=vector();
        while(from < num- homo.step){
            sub.seq=substr(dna.seq, from, from+homo.windows-1 );
            temp=append(temp,c(paste(">",from,sep=""), sub.seq));
            from=from+homo.step;
        }
        rm(dna.seq);
        temp.seq=paste(output.dir,"/temp_homo.fa",sep="");
        write.table(temp,temp.seq,row.names=F,col.names=F,quote=F)

        temp.output=paste(output.dir,"/temp_homo.sam",sep="")
        cmd=paste("bowtie2 -x", genome.ref, "--local  --quiet -f --no-head --reorder -U", temp.seq, "-S",temp.output,sep=" ");
        system(cmd);
        findhomo(output.dir, c(homo.windows, homo.step));
        tag=file.remove(temp.seq)
        tag=file.remove(temp.output)
    }

    # do alignment
    aln.genome=paste(output.dir,"/aln_genome.sam",sep="");
    aln.genome2=paste(output.dir,"/aln_genome.sam.gz",sep="");
    if(!file.exists(aln.genome) & !file.exists(aln.genome2)){
        if(is.null(genome.ref))
            stop("Error: genome.ref should be given as a prefix");
        if(is.null(mate1) | is.null(mate2))
            stop("Error: mate1 and mate2 should be given as fastq files");
        cmd=paste("bowtie2 -x", genome.ref, "--local --reorder --dovetail -X", max.fragment, "-p", cores, "-1", mate1, "-2", mate2,"-S",aln.genome, sep=" ");
        print("Start bowtie2 alignment in genome...");
        system(cmd);
    }

    aln.insert=paste(output.dir,"/aln_insert.sam",sep="");
    aln.insert2=paste(output.dir,"/aln_insert.sam.gz",sep="");
    if(!file.exists(aln.insert) & !file.exists(aln.insert2)){
        if(is.null(mate1) | is.null(mate2))
            stop("Error: mate1 and mate2 should be given as fastq files");

        cmd=paste("bowtie2 -x", insert.ref, "--local --reorder --dovetail -X", max.fragment, "-p", cores, "-1", mate1, "-2", mate2,"-S",aln.insert, sep=" ");
        print("Start bowtie2 alignment in insert...");
        system(cmd);
    }

    if(gz & !file.exists(aln.genome2)){
        system(paste("gzip", aln.genome, sep=" "));
    }
    if(gz & !file.exists(aln.insert2)){
        system(paste("gzip", aln.insert, sep=" "));
    }
    #assign reads
    assign.file=paste(output.dir,"/assign.txt",sep="")
    if(!file.exists(assign.file)){
        print("begin assign reads...");
        assignreads(c(output.dir,backgroud.chrom), c(min.clip,cores,strict,max.tail, block.size))
    }
    # do fragment alignment
    temp.align.genome=paste(output.dir,"/temp_files/fragment_genome.sam",sep="")
    temp.align.insert=paste(output.dir,"/temp_files/fragment_insert.sam",sep="")
    if(!file.exists(temp.align.genome) | !file.exists(temp.align.insert)){
        temp.fragment=paste(output.dir,"/temp_files/temp_fragment.txt",sep="");
        #cmd1=paste("bowtie2 -x", genome.ref, "--quiet --local -f -p", cores, "-U", temp.fragment, "-S", temp.align.genome ,sep=" ");
        #cmd2=paste("bowtie2 -x", insert.ref, "--quiet --local -f -p", cores, "-U", temp.fragment, "-S", temp.align.insert ,sep=" ");
        cmd1=paste("bowtie2 -x", genome.ref, "--quiet --local -f -U", temp.fragment, "-S", temp.align.genome ,sep=" ");
        cmd2=paste("bowtie2 -x", insert.ref, "--quiet --local -f -U", temp.fragment, "-S", temp.align.insert ,sep=" ");
        print("begin second-round assign reads...");
        system(cmd1)
        system(cmd2);
        reassignreads(output.dir, c(allow.gap, block.size, reads.dup.score))
    }

    #file.remove(temp.fragment)
    #file.remove(paste(output.dir,"/temp_mark.txt",sep=""));
    #file.remove(paste(output.dir,"/temp_seq.txt",sep=""));
    #file.remove(paste(output.dir,"/temp_id.txt",sep=""));

    #find the integration sites
    print("Finding the integration sites...");
    split.file=paste(output.dir,"/temp_files/site_split.txt",sep="")
    range.file=paste(output.dir,"/temp_files/range_split.txt",sep="")
    assign.file2=paste(output.dir,"/assign2.txt",sep="")
    if(!file.exists(split.file)){
        print("Error: '/temp_files/site_split.txt' does not exist!")
    }

    splits=unique(fread(split.file, sep="\t",header=TRUE, showProgress=FALSE));
    #assigns=fread(assign.file2,header=TRUE,sep="\t", showProgress=FALSE)
    ranges=as.data.frame(fread(range.file,header=TRUE,sep="\t", showProgress=FALSE))
    splts<-unique(fread(paste(output.dir,"/temp_files/read_split.txt",sep=""), header=T,sep="\t", showProgress=FALSE))

    sites=.predict.recombination.sites(splits, ranges, splts, region=region, min.split.reads=min.split.reads,cores=cores);
    new.sites=sites;

    if(calculate.nonsplit.reads | estimate.copy){
        print("Finding the non-split reads... ")
        ipt=paste("my $dir=\"", output.dir,"\";",sep="")
        cmd=paste(
          ipt,
          "open DATA,\"$dir/temp_files/read_nonsplit.txt\" or die;",
          "my %chrs;",
          "<DATA>;",
          "while(my $str=<DATA>){",
          "    my @str=split(/\t/,$str);",
          "    if(defined($chrs{$str[3]})){",
          "        print {$chrs{$str[3]}} \"$str[0]\t$str[4]\t$str[5]\";",
          "    }",
          "    else{",
          "      open(my $fileout,\">\",\"${dir}/temp_files/read_nonsplit_$str[3].txt\");",
          "      print $fileout \"$str[0]\t$str[4]\t$str[5]\";",
          "      $chrs{$str[3]}=$fileout;",
          "    }",
          "}",
          "close DATA;",

          "foreach (keys %chrs){",
          " close($chrs{$_});",
          "}",
          sep="\n"
        )
        script.file=paste(output.dir,"/temp_files/perlscript.pl",sep="")
        writeLines(cmd, script.file)
        system(paste("perl  ",script.file,sep=""))
        new.sites=.refine.sites(output.dir, sites, min.clip=10, cores=cores);
    }
    predict.file=paste(output.dir,"/report.txt",sep="")
    write.table(new.sites[,names(new.sites)!="reads"], predict.file, row.names=F,quote=F,sep="\t")
    reads=lapply(seq_len(nrow(new.sites)), function(x){
      y=strsplit(as.vector(as.matrix(new.sites[x,"reads"])),",")[[1]];
        return(y)
    })
    names(reads)<-as.vector(new.sites$sites)
    #predict.reads.file=paste(output.dir,"/temp_files/.predict_recombination_reads.rda",sep="")
    #save(reads,file=predict.reads.file)
    sub.sites=subset(new.sites, count >= min.split.reads & score > 80 & balance.len > 50 );

    #predict the integration sites
    print("Begin integration sites prediction...");
    sites.dir=paste(output.dir, "/sites",sep="");
    if(!file.exists(sites.dir))
        dir.create(sites.dir);

    print("Begin drawing read plots...");
    for(st in head(as.vector(sub.sites$sites),500)){
          pdf(paste(output.dir,"/sites/",st,".pdf",sep=""));
          .breakpoints(st, output.dir, new.sites, splts, reads,  span=span, col.left=col.left, col.right=col.right,max.show=max.show)
          dev.off()
    }
    # estimate copy number
    if(estimate.copy){
        print("Begin estimate copy...");
        fragment.estimation(output.dir,  seq.depth=seq.depth, homozygote= FALSE);
    }
}
