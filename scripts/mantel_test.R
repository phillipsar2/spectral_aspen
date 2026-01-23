### Title: Mantel test
### Author: Alyssa Phillips
### Date: 14/4/25

library(vcfR)
library(stringr)
library(vegan)

# Load data ----

prjdir <- "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/"

## Genetics ----
### VCF used in gusrelate for kinmat (thinned and ploidy only data)
# vcf.thin <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/vcf/RMBL.genowithploidy.vcf.gz")
# vcf.thin

file <- "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz"
vcf <- read.vcfR(file)
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL.ploidycalls.2025-04-17.csv")
sum(colnames(vcf@gt) %in% meta$sample)
samples_keep <- colnames(vcf@gt) %in% meta$sample
samples_keep[1] <- TRUE # keep FORMAT column
vcf.sub <- vcf[,samples_keep]
vcf.sub

# Extract vcf values including GT and AD (from Shelly Gaynor)
ccdf1 <- vcfR::vcfR2tidy(vcf.sub,
                         info_fields = c("DP"),
                         format_fields = c("GT", "AD", "DP"))
ccdf1gt <- ccdf1$gt
rm(ccdf1)

## Extract and format Total Depth
df1_dp <- data.frame(POS=paste0(ccdf1gt$ChromKey,".", ccdf1gt$POS),
                     IND=ccdf1gt$Indiv,
                     DP = as.numeric(ccdf1gt$gt_DP))

df1_dp <- df1_dp %>%
  tidyr::pivot_wider(names_from = POS, values_from = DP)
names <- df1_dp$IND
df1_dp <- df1_dp[,-1]
df1_dp <- as.matrix(df1_dp)
row.names(df1_dp) <- names

## Extract and format Allele Depth
ADset  <- data.frame(do.call(rbind, strsplit(ccdf1gt$gt_AD,",")))
ADset$X1 <- as.numeric(ADset$X1)
ADset$X2 <- as.numeric(ADset$X2)

df1_ad <- data.frame(POS=paste0(ccdf1gt$ChromKey,".", ccdf1gt$POS),
                     IND=ccdf1gt$Indiv, 
                     AD =  ADset$X1)

df1_ad <- df1_ad %>%
  tidyr::pivot_wider(names_from = POS, values_from = AD)
names <- df1_ad$IND
df1_ad <- df1_ad[,-1]
df1_ad <- as.matrix(df1_ad)
row.names(df1_ad) <- names

rm(ccdf1gt)

## Divide allele depth by total depth
dfc <- df1_ad/df1_dp
dfc[df1_ad == 0 | is.na(df1_ad)] <- 0
dfc[is.na(df1_dp)] <- NA

## Select one read per genotype
dfc_sr <- apply(dfc, c(1,2), function(x){rbinom(n = 1, size = 1, prob = x)})
dfc_sr[is.na(dfc)] <- NA

## Calculate IBS
# IBS == 0 if 0,1
# IBS == 1 if 1,1 or 0,1
# sum across positions then divide by the number of positions with data for both ind

combn(colnames(dfc_sr), 2) %>% head()

### Mixed-ploidy vcfs from updog
# files <- list.files(path = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/", pattern  = "\\.vcf.gz")
# files <- paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/", files)
# 
# gt <- lapply(files, function(x){read.vcfR(x) %>% extract.gt(element = "GT")}) # extract genotypes
# gt_mat <- do.call(rbind, gt)
# dim(gt_mat)
# rm(gt)
# 
# ad <- lapply(files, function(x){read.vcfR(x) %>% extract.gt(element = "AD")}) # extract allele depth
# ad_mat <- do.call(rbind, ad)
# dim(ad_mat)
# rm(ad)



## grm ----
grm <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/RMBLgrm.filt.2025-05-12.csv", 
                header = T, row.names = 1)
grm[1:5,1:5]
dim(grm)

# sum(colnames(grm) %in% genotypes[,1])

## extract genotypes list
genotypes <- colnames(grm) %>% as.data.frame()
genotypes$sitecode <- genotypes[,1] %>% 
  str_split(pattern = "_", simplify = T) %>%
  as.data.frame() %>%
  dplyr::select(V1)
dim(genotypes)
genotypes <- apply(genotypes,2, as.factor)

## rename grm
colnames(grm) <- genotypes[,2]
rownames(grm) <- genotypes[,2]
grm[1:5,1:5]

## Spectral distance ----
sv_dist <- read.csv(paste0(prjdir, "data/spectral_variation/spectral_distance/is/is_pairwise_spectral_distance_CR.csv"),
                    row.names = 1)
dim(sv_dist)
sv_dist[1:5,1:5]

sum(colnames(sv_dist) %in% genotypes[,2])

sv_dist_tmp <- sv_dist[,colnames(sv_dist) %in%  genotypes[,2]]
sv_dist_sub <- sv_dist_tmp[rownames(sv_dist_tmp) %in%  genotypes[,2],]
dim(sv_dist_sub)

colnames(sv_dist_sub)[!colnames(sv_dist_sub) %in% genotypes[,2]]
genotypes[!genotypes[,2] %in% colnames(sv_dist),]

## Spatial distance matrix ----
geo <- read.csv(paste0(prjdir, "data/spectral_variation/spectral_distance/rmbl_pairwise_spatial_distance_m.csv"),
                row.names = 1)
dim(geo)
geo[1:5,1:5]

sum(colnames(geo) %in% colnames(sv_dist_sub))

geo_tmp <- geo[,colnames(geo) %in%  colnames(sv_dist_sub)]
geo_sub <- geo_tmp[rownames(geo_tmp) %in%  colnames(sv_dist_sub),]
dim(geo_sub)


### Subset grm
grm_tmp <- grm[,colnames(grm) %in%  colnames(sv_dist_sub)]
grm_sub <- grm_tmp[rownames(grm_tmp) %in%  colnames(sv_dist_sub),]
dim(grm_sub)

### Genotypes lacking spectral data due to few aspen pixels:
### "ASIBO"   "BS"      "JCSAT"   "SGBTF01" "VBHJU"

### Genotypes lacking genetic data or unknown ploidy call:
### [1] "AEBVG"   "AQXKL"   "AQZOX"   "CCZV"    "CZUO"    "FYSUQ"   "GJBOV"   "GKPW"    "GPQAO"   "HSYDC"  
### [11] "HVDLX"   "HYXUD"   "KGWC"    "KLZG"    "KNCV"    "MEMRS01" "MGWO"    "NISHT"   "QGVH"    "RJHT"   
### [21] "RNUV"    "SDIP"    "SGYC"    "XJOM"    "YIBK" 

### Final n = 457

# Estimate genetic distance ----
## Average Euclidean distance between continuous genotypes

### Filter missing data
miss <- rowSums(is.na(gt_mat))/dim(gt_mat)[2]
hist(miss*100)

gt_miss <- gt_mat[miss < 0.10,]
gt_miss[1:5,1:5]
dim(gt_miss)
rm(gt_mat)
rm(miss)

### Convert to continuous genotypes
unique(gt_miss[3,])

gt_cont <- apply(gt_miss, c(1,2), function(x) {
  if (is.na(x)){      # is.na needs to be the first test
    "<NA>"              
  } else if (x == "0/0"| x== "0/0/0"){
    "0"
  } else if (x== "0/0/1"){
    "0.33"
  } else if (x == "0/1" ){
    "0.5"
  } else if (x == "0/1/1" ){
    "0.66"
  } else if (x == "1/1" | x== "1/1/1"){
    "1"}})

dim(gt_cont)
gt_cont[1:5, 1:5]

### Thin by LD
gt_num <- apply(gt_cont, 2, as.numeric)
gt_num[1:5,1:5]
rownames(gt_num) <- rownames(gt_cont)
rm(gt_miss)
rm(gt_cont)

#### Calculate MAF and filter
q <- rowMeans(gt_num, na.rm = T)
hist(q)

gt_maf <- gt_num[q > 0.05 & q < 0.95,] # MAF filter of 0.05

#### Correlation among q
library(ldsep)
ld_tri <- mldest(
  gt_maf,
  2,
  nc = 1,
  type = "comp",
  se = F # don't calculate standard errors
)
dim(ld_tri)

# Mantel test ----
## genetic distance vs geographic distance
grm_geo <- mantel(grm_sub, geo_sub, 
                  method = "pearson", 
                  permutations = 9999)

grm_geo # Pearson r = -0.1978, p-value = 1

## genetic distance vs spectral distance
grm_sv <- mantel(grm_sub, sv_dist_sub, 
                  method = "pearson", 
                  permutations = 9999)

grm_sv # pearson r = -0.0698, p-value = 1

# Plot
library(reshape)
library(tidyr)
library(dplyr)
library(data.table)

grm_sub[upper.tri(grm_sub)] <- NA
grm_df = melt(data.table(grm_sub, keep.rownames = TRUE))
colnames(grm_df) <- c('geno1', 'geno2', 'R')
grm_df <- filter(grm_df, geno1 != geno2) %>%
  filter(R != 'NA')
dim(grm_df)
head(grm_df)

sv_dist_sub[upper.tri(sv_dist_sub)] <- NA
sv_dist_df = melt(data.table(sv_dist_sub, keep.rownames = TRUE))
dim(sv_dist_df)
colnames(sv_dist_df) <- c('geno1', 'geno2', 'sv')
sv_dist_df <- filter(sv_dist_df, geno1 != geno2) %>%
  filter(sv != 'NA')
dim(sv_dist_df)
head(sv_dist_df)

grm_df$sv_dist <- sv_dist_df$sv
plot(grm_df$R, grm_df$sv_dist)

library(ggplot2)
ggplot(grm_df, aes(x = R, y = sv_dist)) +
  geom_point() +
  geom_smooth(method = "lm")
