# Title: Updog genotyping
# Author: Alyssa Phillips
# Date: 1/15/2025

devtools::install_github("dcgerard/updog")
library("updog")
devtools::install_github(repo="knausb/vcfR")
library(vcfR)
library(dplyr)
library(stringr)
library("argparser")
library(ldsep)

# Argument name
ap <- arg_parser("Updog genotyping")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "Input VCF")

# add additional arguments
# ap <- add_argument(ap, "--meta", help = "Metadata file")
ap <- add_argument(ap, "--ploidy", help = "choose 'diploid' or 'triploid'")
ap <- add_argument(ap, "--cores", help = "Number of available cores")
ap <- add_argument(ap, "--outdir", help = "Output directory")

# parse arguments
argv <- parse_args(ap)

ploidy_level <- as.character(argv$ploidy)
cores <- as.numeric(argv$cores)
outdir <- as.character(argv$outdir)

# Load VCF
# vcf <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz",
  # verbose = FALSE )
vcf <- read.vcfR(as.character(argv$vcf), verbose = FALSE )

vcf

# Extract allele depth
ad <- extract.gt(vcf, element = "AD", as.numeric = F)

# Create input matrix of depth of reference read
refmat <- masplit(ad, record = 1, sort = 0)

# Create input matrix of total read depth
alt <- masplit(ad, record = 2, sort = 0)

sizemat <- refmat + alt

# Load metadata
# meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020 (1).csv")
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL.ploidycalls.2025-04-17.csv")
meta <- arrange(meta, sample) 
# meta <- read.csv(as.character(argv$meta))

# Split input dataframes into diploids and triploids
# sites <- str_split(colnames(refmat), pattern = "_", simplify = T)[,1]

# key <- cbind(colnames(refmat), sites) %>% as.data.frame()
# colnames(key) <- c("tree", "site")

# key_ploidy <- merge(x = key, y = meta, by.x = "site", by.y = "Site_Code") %>%
  # dplyr::select(site, tree, Ploidy_level)

# Genotype

if (ploidy_level == 'diploid'){
  # dips <- key_ploidy$tree[key_ploidy$Ploidy_level == "Diploid"]
  # print(paste0("Number of diploids: ", length(dips[complete.cases(dips)]) ))
  dips <- meta$sample[meta$ploidy_call == "diploid"]
  print(paste0("Number of diploids: ", length(dips[complete.cases(dips)]) ))
  
  refmat_dips <- refmat[ , colnames(refmat) %in% dips]
  sizemat_dips <- sizemat[,colnames(sizemat) %in% dips]
  
  mout <- multidog(refmat = refmat_dips,
                       sizemat = sizemat_dips,
                       ploidy = 2,
                       model = "norm",
                       nc = cores)

  # write.table(genomat, file=paste0(outdir, "/updog.genomat.diploid.", Sys.Date(),".txt"), quote = F)

  
} else if (ploidy_level == 'triploid'){
  # trips <- key_ploidy$tree[key_ploidy$Ploidy_level == "Triploid"]
  # print(paste0("Number of Triploids: ", length(trips[complete.cases(trips)]) ))
  trips <- meta$sample[meta$ploidy_call == "triploid"]
  print(paste0("Number of triploids: ", length(trips[complete.cases(trips)]) ))
  
  refmat_trips <- refmat[,colnames(refmat) %in% trips]
  sizemat_trips <- sizemat[,colnames(sizemat) %in% trips]
  
  mout <- multidog(refmat = refmat_trips, 
                        sizemat = sizemat_trips, 
                        ploidy = 3, 
                        model = "norm",
                        nc = cores)
  
} else {
  print("ERROR: --ploidy can only be diploid or triploid")
}

# Examine genotype and SNP Quality
pdf(paste0(outdir,"/genotype_depth_distributions.",ploidy_level,".",Sys.Date(),".pdf"))
## The (posterior) proportion of individuals mis-genotyped at each site
hist(mout$snpdf$prop_mis, main = "Proportion of ind mis-genotyped at each SNP")

## Overdispersion of each snp - simulations suggest dropping > 0.05
hist(mout$snpdf$od, main = "Overdispersion of each SNP")

## Bias - simulations suggest filtering 0.5 < x > 2
hist(mout$snpdf$bias, main = "Bias at each SNP")

dev.off()

## Filter SNPs based on updog recommendations
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > 0.5 & bias < 2)

# Extract genotype matrix
genomat <- format_multidog(mout_cleaned, varname = "geno")

# Save filtered genotype matrix
write.table(genomat, 
            file = paste0(outdir, "/updog.genomat.", ploidy_level, ".", Sys.Date(),".txt"), 
            quote = F)

## Estimate LD
ploidy <- ifelse(ploidy_level == "diploid", 2, 3) 
varnames <- paste0("logL_", 0:ploidy)
larray <- format_multidog(x = mout_cleaned, varname = varnames) # get genotype log-likelihoods

like_ld <- mldest(geno = larray, K = ploidy, type = "comp") # estimate LD

ldmat <- format_lddf(obj = like_ld, element = "r2") # make into matrix

write.table(ldmat, 
            file = paste0(outdir, "/updog.LDmat.", ploidy_level, ".", Sys.Date(),".txt"), 
            quote = F)


##### Test code below
# dips <- key_ploidy$tree[key_ploidy$Ploidy_level == "Diploid"]
# paste0("Number of diploids: ", length(dips[complete.cases(dips)]) )
# trips <- key_ploidy$tree[key_ploidy$Ploidy_level == "Triploid"]
# paste0("Number of triploids: ", length(trips[complete.cases(trips)]) )

# refmat_dips <- refmat[ , colnames(refmat) %in% dips]
# refmat_trips <- refmat[,colnames(refmat) %in% trips]

# sizemat_dips <- sizemat[,colnames(sizemat) %in% dips]
# sizemat_trips <- sizemat[,colnames(sizemat) %in% trips]

# Run multidog - Diploids
# mout_dip <- multidog(refmat = refmat_dips, 
#                      sizemat = sizemat_dips, 
#                      ploidy = 2, 
#                      model = "norm",
#                      nc = 4)

# write.table(genomat, file=paste0("~/aspen/radseq/erincar/updog/updog.genomat.diploids.", Sys.Date(),".txt"), quote = F)

# Run multidog - tripoids
# mout_trip <- multidog(refmat = refmat_trips, 
#                       sizemat = sizemat_trips, 
#                       ploidy = 3, 
#                       model = "norm",
#                       nc = 4)

# write.table(genomat, file=paste0("~/aspen/radseq/erincar/updog/updog.genomat.diploids.", Sys.Date(),".txt"), quote = F)

# Test triploids with one snp
# test_trip <- flexdog(refvec  = refmat_trips[2,], 
#                      sizevec = sizemat_trips[2,], 
#                      ploidy  = 3, 
#                      model   = "norm")

# 2.270 seconds per SNP
# system.time(flexdog(refvec  = refmat_trips[1,], 
#                     sizevec = sizemat_trips[1,], 
#                     ploidy  = 3, 
#                     model   = "norm"))

# plot_geno(refvec = refmat_trips[2,], sizevec = sizemat_trips[2,], ploidy = 3,)
# plot(test_trip)


# Visualize output
# plot(mout_dip, indices = c(500, 5, 100))
# 
# str(mout$snpdf)
# 
# str(mout$inddf)
