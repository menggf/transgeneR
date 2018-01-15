---
title: 'TransgeneR: A tool for transgenic integration and recombination sites discovery'
author: "Guofeng Meng"
date: '2017-06-24'
output:
  pdf_document:
    keep_tex: yes
    number_sections: yes
    toc: yes
    toc_depth: 2
  html_document:
    toc: yes
    toc_depth: '2'
vignette: "\\VignetteIndexEntry{TransgeneR: A tool for transgenic integration and
  recombination sites discovery} \n\\VignetteEncoding{UTF-8} \n\\VignetteEngine{knitr::rmarkdown}\n"
---

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{TransgeneR: A tool for transgenic integration and recombination sites discovery}
-->

## Introduction

TransgeneR is designed to find the transgenic integration information in the animal genome using the whole genome sequencing data or PCR-based sequencing data. In many case, the transgenic sequences can have multiple integration sites and even recombination between transgenic sequences. Therefore, transgeneR is supposed to answer following question:
* where are the transgenic sequences integrated on the genome?
* Is there transgenic recombination?
* How many transgenic sequences are integrated on the genome?

To use transgeneR:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("transgeneR")
```

Please note that "bowtie2" should be installed in the users' computers. And the bowtie2 geneome reference has been built the studied animals before using transgeneR. 

TransgeneR is an one-stop analysis pipeline. 

The output of transgeneR are stored in a directory set by "output.dir". 
