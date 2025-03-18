# Title: SNP Quality Assessment - RADSEQ
# Author: Alyssa Phillips
# Date: 11/27/2024

library(dplyr)
library(ggplot2)

##############
### Raw snps
##############

# Load quality table

dir <- as.character("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/")
files <- list.files(path = dir, pattern = ".table")

# extract_data_from_file <- function(file) { 
#   data <- read.table(file, header = TRUE, 
#                    stringsAsFactors = FALSE)
#   return(data) 
# } 
# 
# qual <- lapply(paste0(dir,files), 
#                     extract_data_from_file) 
# 


qual <- read.table(paste0(dir, files[1]), header = T)
dim(qual)
str(qual)

pdf(paste0(dir, "alignment_quality.chr1.",Sys.Date(),".pdf"))

qual %>%
  ggplot(aes(x=QUAL)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  xlim(c(0,150)) +
  theme_bw() +
  geom_vline(xintercept = 20, color = "black")+
  geom_vline(xintercept = 30, color = "black") +
  xlab("Base quality (QUAL)")

qual %>%
  ggplot(aes(x=MQ)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  xlim(c(0,65)) + # bwa-mem doesn't go above MQ of 60
  theme_bw() +
  geom_vline(xintercept = 20, color = "black")+
  geom_vline(xintercept = 30, color = "black") +
  xlab("Mapping quality (MQ)")

qual %>%
  ggplot(aes(x=DP)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef") +
  xlim(c(0,50)) +
  theme_bw() +
  xlab("Total depth at each site (DP)")

dev.off()

# Playing around with filter thresholds
dim(qual)

dplyr::filter(qual, QUAL > 20, MQ > 20) %>%
  dim()

dplyr::filter(qual, QUAL > 30, MQ > 50) %>%
  dim()

dplyr::filter(qual, QUAL > 30, MQ > 50, DP < 50, DP > 10) %>%
  dim()

dplyr::filter(qual, QUAL > 30, MQ > 50, DP < 50, DP > 20) %>%
  dim()

##############
### Filtering for Depth
##############
# Load quality table
dir <- as.character("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/depth/")
files <- list.files(path = dir, pattern = ".table")

qual <- read.table(paste0(dir, files[1]), header = T)
# qual <- read.table("~/aspen/radseq/erincar/rad_aspen.all.filtered.nocall.table" , header = T, skip = 1)
head(qual)
str(qual)
dim(qual)

qual[1:10,1:10]

# Exclude HSYDC_448_18.DP
# qual <- qual[,colnames(qual) != "HSYDC_448_18.DP"]

genotypes <- dim(qual)[2] - 2

# Prep table
# qual <- qual[qual$CHROM != "CHROM",] # shouldn't be necessary if files are properly joined

qual[,3:dim(qual)[2]] <- lapply(qual[,3:dim(qual)[2]], as.numeric)
str(qual)
dim(qual)

# Replace 0s with NA
# qual[qual == 0] <- "NA"
# dim(qual)

# Estimate mean genotype depth for sites with coverage
geno_dp <- colMeans(qual[3:dim(qual)[2]], na.rm = T)

## Plot mean genotype depth
p_gdp <- reshape2::melt(geno_dp) %>%
  ggplot(aes(x=value)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  xlim(c(0,max(geno_dp, na.rm = T))) +
  theme_bw() +
  geom_vline(xintercept = mean(geno_dp, na.rm = T), color = "black") +
  labs(x = "Average genotype depth", 
       title = paste0("Average genotype depth for ", genotypes, " genotypes"))

ggsave(plot = p_gdp, filename = paste0(dir, "avg_geno_dp.", Sys.Date(),".jpg"),
       width = 5, height = 4, units = "in")

# Plot genotype depth distributions
qlist <- matrix(nrow = genotypes, ncol = 3) # define number of samples (10 samples here)
qlist <- data.frame(qlist, row.names=colnames(qual)[-c(1:2)])

pdf(paste0(dir, "genotype_depth_distributions.",Sys.Date(),".pdf"))
par(mfrow=c(4,3)) # mfrow sets max number of plots (rows by columns), automatically makes new pages if PDF

for (i in 3:dim(qual)[2]) {
  qlist[i-2,] <- quantile(qual[,i], c(.05, .1, .99), na.rm=T)
  d <- density(qual[,i], from=0, to=100, bw=1, na.rm = T)
  plot(d, xlim = c(0,50), main=rownames(qlist)[i-2], col="blue", lwd=2, xlab = NULL)
  abline(v=qlist[i-2,c(1,3)], col='red', lwd=3)
}

dev.off()

# # Test out some depth + missingness filters
# max = 30
# min = 6
# miss = 0.1
# 
# qual[qual < 6] <- "NA"
# qual[qual > max] <- "NA"
# 
# genos_with_data <- rowSums(is.na(qual[,3:dim(qual)[2]])) # number of genotypes to be removed at each site
# hist(genos_with_data)

# Depth per site
site_dp <- rowSums(qual[3:dim(qual)[2]], na.rm = T)


# hist(site_dp, xlim = c(0, 400), breaks = 50)
p_sdp <- as.data.frame(site_dp) %>% 
  ggplot(aes(x=site_dp)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.9, ) +
  xlim(c(0,50)) +
  theme_bw() +
  geom_vline(xintercept = mean(site_dp, na.rm = T), color = "black") +
  labs(x = "Average genotype depth at each site", 
       title = paste0("Average genotype depth at each site for ", length(site_dp), " sites"))

ggsave(plot = p_sdp, filename = paste0(dir, "avg_geno_dp.", Sys.Date(),".jpg"),
       width = 5, height = 4, units = "in")

# Genotypes without data (not zero)
miss <- rowSums(qual[,3:dim(qual)[2]] == 0)
hist(miss/genotypes * 100,
     xlab = "Missing data per site", ylab = "Number of sites",
     breaks = 20)

sum(miss <= 10)
