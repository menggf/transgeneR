# transgeneR
This a a package designed for trangenic integration and rearrangement site discovery using sequencing data

## Install TransgeneR
library(devtools)

install_github("menggf/transgeneR")

## usage
test.data=system.file("extdata", "test.zip", package = "transgeneR")

temp <- tempdir()

unzip(test.data, exdir=temp)

output.dir=paste(temp,"/test",sep="")

transgene.sq=paste(temp,"/test/temp_files/insert.fa",sep="")

transgeneR(output.dir,  insert.seq=transgene.sq)
