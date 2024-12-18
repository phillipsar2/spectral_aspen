## Title: nQuack
## Author: Alyssa Phillips
## Date: 12/16/24

library(nQuack)
library(dplyr)
library(kableExtra)


## Prepare many samples
inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/filtered/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/"
filelist <- list.files(path = inpath, pattern = "*.bam" )
filelist <- gsub(".bam", "", filelist)

for( i in 1:length(filelist)){
  prepare_data(filelist[i], inpath, outpath)
}
