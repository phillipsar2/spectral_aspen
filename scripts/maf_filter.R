### Title: MAF filter
### Author: Alyssa Phillips
### Date: 9/24/2025

library(vcfR)
library(dplyr)
library(stringr)
library(tidyverse)

vcf <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.Chr19.vcf.gz",
                 verbose = FALSE)
vcf

# Extract genotypes
gt <- extract.gt(vcf, element = "GT", as.numeric = F)
gt[1:5,1:5]

# Convert genotypes to numbers
## "1/1/1" is AAA (ref homozygous)
str <- list()
str$genotypes <- c( "0/0",   "0/1" , "1/1", "0/0/0", "0/0/1", "0/1/1", "1/1/1")
str$num_genos <- c("0","1","2","0","1","2","3")

gt_num <- apply(gt, 2, 
                function(x){as.character(factor(x, levels = str$genotypes, labels = str$num_genos))})
gt_numeric <- apply(gt_num, c(1,2), as.numeric)
# dim(gt_numeric)
# rbind(gt[2,], gt_num[2,]) # verify it worked

# Get ploidy of each genotype
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL.ploidycalls.2025-04-17.csv")
colnames(gt_numeric) == meta$sample

ploidy_sorted <- c()
for (p in 1:length(colnames(gt_numeric))){
  ploidy_sorted[p] <- meta$ploidy_call[colnames(gt_numeric)[p] == meta$sample]
}

sum(is.na(ploidy_sorted)) # should be zero
ploidy_num <- ifelse(ploidy_sorted == "diploid", 2, 3)

# Calculate MAF per row
allele_freq <- function(row){
  max_alleles <- sum(ploidy_num[!is.na(row)]) # take individuals that have data at the site and sum to get max alleles observed
  p_count <- sum(row, na.rm = T) # sum allele observations
  p_count/max_alleles
}
# allele_freq(gt_numeric[8,])

p <- c()
p <- apply(gt_numeric, 1, allele_freq)
length(p)
hist(p)

p_df <- data.frame("SNP" = rownames(gt), 
           "p" = p,
           "Chr" = NA,
           "POS" = NA)
p_df[,3:4] <- str_split(p_df$SNP, "_", simplify = T)
head(p_df)

# filter
hist(p[p > 0.05 & p < 0.95])
length(p[p > 0.05 & p < 0.95])

p_filt <- p_df[p_df$p > 0.05 & p_df$p < 0.95,]
hist(p_filt$p, breaks = 20)

# export sites
write.table(p_filt[,3:4], "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.Chr19.mafSNPs.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

