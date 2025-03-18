### Title: gbs2ploidy
### Author: Alyssa
### Date: 3/4/25

library('rjags')
library(gbs2ploidy)
library(vcfR)
# devtools::install_github("DevonDeRaad/SNPfiltR")
library(SNPfiltR)


# (1) Load vcf
vcf <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz",
                 verbose = FALSE )
vcf

# (2) Additional VCF filters

## minimum coverage 2x per genotype
vcf <- hard_filter(vcf = vcf, depth = 8)

## minimum 8 reads per alt allele at a site 
# vcf <- min_mac(vcf, min.mac = 8)

# (3) Extract read depths for alleles at heterozygous sites

## Extract allele depth
ad <- extract.gt(vcf, element = "AD", as.numeric = F)

## MAF filter of 0.1
maf <- maf(vcf)
hist(maf[,4])

dim(ad)
ad <- ad[maf[,4] > 0.1,]
dim(ad) # final number of SNPs

## Create input matrix of depth of reference read
refmat <- masplit(ad, record = 1, sort = 0)

## Create input matrix of alt read depth
altmat <- masplit(ad, record = 2, sort = 0)

## Plot allele ratio
ar <- altmat/(altmat + refmat)
dim(ar)
ar[1:5,1:5]

par(mfrow=c(3,3))
for(i in 1:9){ # change number from 1 to 9
  hist(ar[,i], 
       # axes=FALSE, 
       xlab="Observed allelic ratio", 
       ylab="Count",
       main = colnames(ar)[i],
       # xlim = c(0.1,0.9), 
       breaks = 20)
}

## missing data?
miss <- colSums(is.na(ar))
miss[1:9]/dim(ar)[1] * 100
hist(miss/dim(ar)[1])

## Average depth per genotype
dp <- extract.gt(vcf, element = "DP", as.numeric = T) %>% colMeans(na.rm = T)
hist(dp, breaks = 20, xlim=c(0,35))

# (4) Infer ploidy 
##  Bayesian inference of allelic proportions
## began at 8:39
propOut <- estprops(cov1 = altmat[,1:3], # [SNPs, ind], non-ref
                    cov2 = refmat[,1:3], # ref allele
                    props = c(0.33, 0.5, 0.66),
                    mcmc.nchain = 3, 
                    mcmc.steps = 1000, mcmc.burnin = 100, mcmc.thin = 2)



## allelic proportions for the first nine individuals
geno_names <- colnames(vcf@gt)[-1]
pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/nquack/gbs2ploidy_alleleratios.pdf")
par(mfrow=c(3,3))
for(i in 1:3){ # change number from 1 to 9
  plot(propOut[[i]][,5], ylim=c(0,1), axes=FALSE, 
       xlab="allelic ratios", 
       # ylab="75th quantile of the posterior distribuion of allelic proportions",
       ylab="proportions",
       main = geno_names[i])
  axis(1, at = 1:3,c("1:2","1:1","2:1"))
  axis(2)
  box()
  segments(1:5, propOut[[i]][,1], 1:5,propOut[[i]][,5])
  # title(main=paste("true ploidy =",dat[[3]][i]))
}
dev.off()

## PCA & DA to assign ploidies
pout <- estploidy(alphas = propOut, 
                  het = rep(0.15, 9), 
                  depth = dp[1:9], 
                  nclasses = 2, # number of cytotypes expected
                  ids = geno_names[1:9], pcs = 1:2)

plot(pout$pcscrs[,1] ~ pout$pcscrs[,2])
