# transgeneR
This a a package designed for trangenic integration and rearrangement sites discovery using sequencing data

## Install TransgeneR
library(devtools)

install_github("menggf/transgeneR")

## prerequisite

TransgeneR uses bowtie2 as the alignment tool. Therefore, before usage, bowtie2 must have been installed
in the linux mechine and its location has been added to $PATH variable. Meanwhile, the genome reference
index has been constructed using "bowtie2-build" command.

## usage
test.data=system.file("extdata", "test.zip", package = "transgeneR")

temp <- tempdir()

unzip(test.data, exdir=temp)

output.dir=paste(temp,"/test",sep="")

transgene.sq=paste(temp,"/test/temp_files/insert.fa",sep="")

transgeneR(output.dir,  insert.seq=transgene.sq)

## Q & A
1. error: /usr/bin/ld: cannot find -lboost_iostreams

Please install boost library in your system using the command like:
yum install boost-devel
