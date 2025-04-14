### Title: gbs2ploidy
### Author: Alyssa
### Date: 3/4/25

library('rjags') # v4.3.2
library(gbs2ploidy)
library(vcfR) #v1.15
# devtools::install_github("DevonDeRaad/SNPfiltR")
library(SNPfiltR) # v.1.0.2
library(ggplot2)
library(dplyr)
library(tidyr)
library(argparser)

# Argument name
ap <- arg_parser("gbs2ploidy")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "gzipped vcf subset to one sample")
# ap <- add_argument(ap, "--tmp", help = "Path to the temp directory")
ap <- add_argument(ap, "--out", help = "Outfile name")

# parse arguments
argv <- parse_args(ap)

vcffile <- as.character(argv$vcf)
out <- as.character(argv$out)

# (1) Load vcf
# vcf <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz",
                 # verbose = FALSE )
# vcf <- read.vcfR("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.1dp30.per0.25.vcf.gz",
                 # verbose =F)
vcf <- read.vcfR(vcffile, verbose = F)
vcf

# Load metadata for flow cytometry samples
# flow_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/known_flow_with_accnum.csv")
# head(flow_meta)

# (2) Extract heterozygous sites from vcf
gt <- extract.gt(vcf)
hets <- is_het(gt)
# Censor non-heterozygous positions by changing them to NAs
is.na(vcf@gt[,-1][!hets]) <- TRUE

# (3) Additional VCF filters

## Minimum coverage of 6x per site per genotype
vcf <- hard_filter(vcf = vcf, depth = 6)
vcf

# (3) Extract read depths for alleles
## Extract allele depth
ad <- extract.gt(vcf, element = "AD", as.numeric = F)
dim(ad)

# ad <- ad[,1]

## Subset flow cytometry samples for testing
# flow_samp <- colnames(ad) %in% flow_meta$accession
# ad <- ad[, flow_samp]
# dim(ad)

## Exclude genotypes with less than 10K sites
# sites <- colSums(!is.na(ad))
# hist(sites, xlab = "sites with data", breaks = 10)

## Write samples that failed filtering to a file
if ( sum(!is.na(ad)) < 10000 ){
  write.table("fail", out,
              row.names = F, col.names = F,)
} else {
  refmat <- masplit(as.matrix(ad), record = 1, sort = 0)
  altmat <- masplit(as.matrix(ad), record = 2, sort = 0)
  
  propOut <- estprops(cov1 = as.matrix(altmat), # [SNPs, ind], non-ref
                      cov2 = as.matrix(refmat), # ref allele
                      props = c(0.33, 0.5, 0.66),
                      mcmc.nchain = 3, 
                      mcmc.steps = 1000, mcmc.burnin = 100, mcmc.thin = 2)
  
  names(propOut) <- colnames(ad)
  propOut_long <- lapply(propOut, function(x){ as.data.frame(x) %>%
      tibble::rownames_to_column(var = "allelic_ratio") %>%
      pivot_longer(cols = `2.5%`:`97.5%`, names_to ="quantile", values_to = "proportion") } )
  
  propOut_df <- reshape2::melt(propOut_long) %>%
    dplyr::select(quantile, value, L1, allelic_ratio)
  colnames(propOut_df) <- c("quantile", "proportion", "sample", "allelic_ratio")
  write.csv(propOut_df, out, row.names = F)
}

# fail <- colnames(ad)[sites < 10000]
# write.table(fail, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/rad_aspen_predictions.mindp6.10ksites.droppedsamples.csv", 
#           row.names = F, sep = ",", col.names = F,)

## Apply filter
# ad <- ad[,sites >= 10000]
# dim(ad)

## Create input matrix of depth of reference read
# refmat <- masplit(ad, record = 1, sort = 0)
# avg_ref_depth <- colMeans(refmat, na.rm = T)

## Create input matrix of alt read depth
# altmat <- masplit(ad, record = 2, sort = 0)
# avg_alt_depth <- colMeans(altmat, na.rm = T)

## Plot allele ratio
# ar <- altmat/(altmat + refmat)

# pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/rad_aspen.mindp6.10ksites.raw_allelicratio.plots.pdf")
# par(mfrow=c(4,4))
# for(i in 1:dim(ar)[2]){ # change number from 1 to 9
#   hist(ar[,i], 
#        # axes=FALSE, 
#        xlab="Observed allelic ratio", 
#        ylab="Count",
#        main = colnames(ar)[i],
#        # xlim = c(0.1,0.9), 
#        breaks = 20)
# }
# dev.off()


## Average depth per genotype
# dp <- extract.gt(vcf, element = "DP", as.numeric = T) %>% 
#   colMeans(na.rm = T) 
# dp <- dp[flow_samp]
# hist(dp, breaks = 20, xlim=c(0,35))

## Combine all dataframe stats
# qual_stats <- data.frame(dp, sites, miss, avg_alt_depth, avg_ref_depth)
# qual_stats$sample <- rownames(qual_stats)

# (4) Infer ploidy ----
##  Bayesian inference of allelic proportions
# propOut <- estprops(cov1 = as.matrix(altmat[,1]), # [SNPs, ind], non-ref
#                     cov2 = as.matrix(refmat[,1]), # ref allele
#                     props = c(0.33, 0.5, 0.66),
#                     mcmc.nchain = 3, 
#                     mcmc.steps = 1000, mcmc.burnin = 100, mcmc.thin = 2)

## Export propOut object. Need to convert list to table
# names(propOut) <- colnames(ad)
# propOut_long <- lapply(propOut, function(x){ as.data.frame(x) %>%
#     tibble::rownames_to_column(var = "allelic_ratio") %>%
#     pivot_longer(cols = `2.5%`:`97.5%`, names_to ="quantile", values_to = "proportion") } )
# 
# propOut_df <- reshape2::melt(propOut_long) %>%
#   dplyr::select(quantile, value, L1, allelic_ratio)
# colnames(propOut_df) <- c("quantile", "proportion", "sample", "allelic_ratio")
# write.csv(propOut_df, out, row.names = F)


## Import propOut values ----
# propOut_df <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/flowcyt_predictions.allsnps.propOut.csv")
# 
# str(propOut_df)
# dim(propOut_df)
# 
# pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/flowcyt_predictions.allelicratio.allsnps.mind6.minsites10k.plots.pdf")
# propOut_df %>%
#   dplyr::filter(quantile == '75%') %>%
#   ggplot(aes(x = as.factor(allelic_ratio), y = proportion)) +
#   geom_point() +
#   facet_wrap( ~ sample, nrow = 9) +
#   theme_bw()+
#   xlab("Allelic ratio") +
#   ylab("75th quantile of the posterior distribuion of allelic proportions")
# dev.off()
# 
# 
# ## allelic proportions for the first nine individuals
# # geno_names <- colnames(ad)
# # pdf("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/nquack/gbs2ploidy_alleleratios.pdf")
# # par(mfrow=c(3,3))
# # for(i in 1:9){ # change number from 1 to 9
# #   plot(propOut[[i]][,5], ylim=c(0,1), axes=FALSE, 
# #        xlab="allelic ratios", 
# #        # ylab="75th quantile of the posterior distribuion of allelic proportions",
# #        ylab="proportions",
# #        main = geno_names[i])
# #   axis(1, at = 1:3,c("1:2","1:1","2:1"))
# #   axis(2)
# #   box()
# #   segments(1:5, propOut[[i]][,1], 1:5,propOut[[i]][,5])
# #   # title(main=paste("true ploidy =",dat[[3]][i]))
# # }
# # dev.off()
# 
# 
# ## What is the winner
# winner_prop <- propOut_df %>%
#   group_by(sample, quantile) %>%
#   filter(proportion == max(proportion)) %>%
#   arrange(sample, quantile)
# 
# winner_num <- winner_prop %>%
#   group_by(sample) %>%
#   count(allelic_ratio, name = 'Count')
# 
# winner_num$ploidy_call <- ifelse(winner_num$allelic_ratio == '0.5', yes = 'diploid', no = 'triploid')
# head(winner_num)
# 
# accuracy <- c()
# for (i in 1:dim(winner_num)[1]){
#   true_ploidy <- flow_meta$flow_ploidy[which(flow_meta$accession == winner_num$sample[i])]
#   accuracy[i] <- true_ploidy == winner_num$ploidy_call[i]
# }
# sum(accuracy)/length(accuracy) * 100
# 
# propcalls <- merge(winner_num, flow_meta, by.x = 'sample', by.y = 'accession') 
# dim(propcalls)so
# write.csv(propcalls, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/flowcyt_predict.propOut.allsnps.mindp6.minsites10k.accuracy.csv",
#           row.names = F)
# 
# # Assess quality stats
# incorrect <- propcalls[!propcalls$ploidy_call == propcalls$flow_ploidy,]
# 
# qual_stats$answer <- ifelse(qual_stats$sample %in% incorrect$sample, yes = "wrong", "right")
# plot(qual_stats$sites, qual_stats$avg_alt_depth,
#      col = as.factor(qual_stats$answer),
#      pch = propcalls$flow_ploidy,
#      xlab = "Number of sites",
#      ylab = "Mean alternate read depth")
# 
# plot(qual_stats$sites, qual_stats$dp,
#      col = as.factor(qual_stats$answer),
#      pch = propcalls$flow_ploidy,
#      xlab = "Number of sites",
#      ylab = "Mean read depth")
# 
# plot(qual_stats$avg_alt_depth, qual_stats$avg_ref_depth,
#      col = as.factor(qual_stats$answer),
#      pch = propcalls$flow_ploidy,
#      xlab = "Mean alternate read depth",
#      ylab = "Mean reference read depth")
# 
# ## PCA & DA to assign ploidies
# pout <- estploidy(alphas = propOut, 
#                   het = rep(0.15, 9), 
#                   depth = dp[1:9], 
#                   nclasses = 2, # number of cytotypes expected
#                   ids = geno_names[1:9], pcs = 1:2)
# 
# plot(pout$pcscrs[,1] ~ pout$pcscrs[,2])
