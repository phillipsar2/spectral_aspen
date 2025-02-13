# Title: RMBL PCA
# Author: Alyssa Phillips
# Date: 1/28/25

library(vcfR)
library(stringr)
library(ggplot2)

# (1) Load genotype matrices from updog ----
dip_genotypes <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.2025-01-27.txt")
dim(dip_genotypes)

trip_genotypes <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.2025-01-27.txt")
dim(trip_genotypes)
trip_genotypes[1:5,1:5]

# (2) Load in metadata ----
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv")
str(meta)

# (3) Prep the data ----
## Merge dataframes to shared sites
sum(rownames(dip_genotypes) %in% rownames(trip_genotypes)) # number of sites
dip_sub <- dip_genotypes[rownames(dip_genotypes) %in% rownames(trip_genotypes),]
trip_sub <- trip_genotypes[rownames(trip_genotypes) %in% rownames(dip_sub),]

dip_sub[1:5,1:5]
trip_sub[1:5,1:5]

genos <- cbind(dip_sub, trip_sub) %>% as.matrix()
dim(genos)

## Missing data across sites
is.na(genos) %>% sum

hist(as.matrix(dip_sub))
hist(as.matrix(trip_sub))

# (4) Run PCA ----
pca_out <- prcomp(genos, scale. = T, center = TRUE)
pve <- summary(pca_out)$importance[2,1:10] *100


# (4) Plot PCA ----
sites <- str_split(colnames(genos), pattern = "_", simplify = T)[,1] %>% 
  as.data.frame()
meta_sub <- merge(x = sites, y = meta, by.y = "Site_Code", by.x = ".", all.x = TRUE, sort = F)
meta_sub$ind <- colnames(genos)

pca_df <- data.frame(ind = meta_sub$ind,
                     ploidy = meta_sub$Ploidy_level,
                     geno = meta_sub$Genotype,
                     pca_out$rotation[,1:10])

ggplot(pca_df, aes(x = PC1, y = PC2, color = geno, shape = ploidy)) +
  geom_point() +
  theme_bw() +
  guides(color="none") +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)"))

ggsave("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/pca/all.bygenotypes.byploidy.pca.jpeg",
       width = 6, height = 4, unit = "in")

ggplot(pca_df, aes(x = PC1, y = PC2, color = ploidy)) +
  geom_point() +
  theme_bw() +
  guides(color="none") +
  xlab(paste0("PC1 (", round(pve[1],2), "%)")) +
  ylab(paste0("PC2 (", round(pve[2],2), "%)"))
