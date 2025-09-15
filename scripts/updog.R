# Title: Updog genotyping
# Author: Alyssa Phillips
# Date: 1/15/2025

# devtools::install_github("dcgerard/updog")
library("updog")
#devtools::install_github(repo="knausb/vcfR")
library(vcfR)
library(dplyr)
library(stringr)
library("argparser")
library(ldsep)
library(tidyr)

# Argument name
ap <- arg_parser("Updog genotyping")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "Input VCF")

# add additional arguments
# ap <- add_argument(ap, "--meta", help = "Metadata file")
ap <- add_argument(ap, "--ploidy", help = "choose 'diploid' or 'triploid'")
ap <- add_argument(ap, "--cores", help = "Number of available cores")
ap <- add_argument(ap, "--outdir", help = "Output directory")
ap <- add_argument(ap, "--chr", help = "Chromosome")

# parse arguments
argv <- parse_args(ap)

ploidy_level <- as.character(argv$ploidy)
cores <- as.numeric(argv$cores)
outdir <- as.character(argv$outdir)
chr <- as.character(argv$chr)

# Load VCF ----
# vcf <- read.vcfR("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.Chr09.RMBL.filtered.gvcf.gz",
#                  verbose = FALSE, nrows = 10000 )
vcf <- read.vcfR(as.character(argv$vcf), verbose = FALSE )

vcf

# Extract allele depth
ad <- extract.gt(vcf, element = "AD", as.numeric = F)

# Create input matrix of depth of reference read
refmat <- masplit(ad, record = 1, sort = 0)

# Create input matrix of total read depth
alt <- masplit(ad, record = 2, sort = 0)

sizemat <- refmat + alt

# Load metadata ----
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

# Genotype ----

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

# Examine genotype and SNP Quality ----
pdf(paste0(outdir,"/genotype_depth_distributions.",ploidy_level,".",chr,".pdf"))
## The (posterior) proportion of individuals mis-genotyped at each site
hist(mout$snpdf$prop_mis, main = "Proportion of ind mis-genotyped at each SNP")

## Overdispersion of each snp - simulations suggest dropping > 0.05
hist(mout$snpdf$od, main = "Overdispersion of each SNP")

## Bias - simulations suggest filtering 0.5 < x > 2
hist(mout$snpdf$bias, main = "Bias at each SNP")

dev.off()

## Filter SNPs based on updog recommendations ----
# prop_miss = posterior proportion of ind mis-genotyped
# bias = allelic bias parameter
mout_cleaned <- filter_snp(mout, prop_mis < 0.05 & bias > 0.5 & bias < 2)

# Extract genotype matrix
genomat <- format_multidog(mout_cleaned, varname = "geno")

# Save filtered genotype matrix
write.table(genomat, 
            file = paste0(outdir, "/updog.genomat.", ploidy_level, ".", chr,".txt"), 
            quote = F)

# Export as VCF ----
library(BIGr)
# multidog.object <- mout
# output.file <- "text.vcf"
updog2vcf <- function (multidog.object, output.file, updog_version = NULL, 
          compress = TRUE) {
  mout <- multidog.object
  ploidy <- as.numeric(unique(multidog.object$snpdf$ploidy))
  if (!grepl(".vcf", output.file)) 
    output.file <- paste0(output.file, ".vcf")
  model_selected <- unique(multidog.object$snpdf$model)
  updog_meta <- paste0("##UpdogCommandLine.multidog=<ID=Multidog,Version=\"", 
                       updog_version, "\",CommandLine=\"> multidog(refmat = matrices$ref_matrix, sizemat = matrices$size_matrix, ploidy = ", 
                       ploidy, ", model = ", model_selected, ")\">")
  bigr_meta <- paste0("##BIGrCommandLine.updog2vcf=<ID=updog2vcf,Version=\"", 
                      packageVersion("BIGr"), "\",Data=\"", Sys.time(), "\", CommandLine=\"> updog2vcf(", 
                      deparse(substitute(multidog.object)), ",", output.file, 
                      ",", updog_version, ")\">")
  vcf_header <- c("##fileformat=VCFv4.3", "##reference=NA", 
                  "##contig=<ID=NA,length=NA>", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", 
                  "##INFO=<ID=ADS,Number=R,Type=Integer,Description=\"Depths for the ref and each alt allele in the order listed\">", 
                  "##INFO=<ID=BIAS,Number=1,Type=Float,Description=\"The estimated allele bias of the SNP from updog\">", 
                  "##INFO=<ID=OD,Number=1,Type=Float,Description=\"The estimated overdispersion parameter of the SNP from updog\">", 
                  "##INFO=<ID=PMC,Number=1,Type=Float,Description=\"The estimated proportion of individuals misclassified in the SNP from updog\">", 
                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype, where 1 is the count of alternate alleles\">", 
                  "##FORMAT=<ID=UD,Number=1,Type=Integer,Description=\"Dosage count of reference alleles from updog, where 0 = homozygous alternate\">", 
                  "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", 
                  "##FORMAT=<ID=RA,Number=1,Type=Integer,Description=\"Reference allele read depth\">", 
                  "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", 
                  "##FORMAT=<ID=MPP,Number=1,Type=Float,Description=\"Maximum posterior probability for that dosage call from updog\">", 
                  updog_meta, bigr_meta)
  if (any(model_selected %in% c("f1", "f1pp", "s1", "s1pp"))) {
    mout$snpdf$maxpostprob <- apply(mout$snpdf[, grep("Pr_", 
                                                      colnames(mout$snpdf))], 1, max)
    if (model_selected == "f1" | model_selected == "f1pp") {
      sele_col <- c("snp", "p1ref", "p1size", "p2ref", 
                    "p2size", "p1geno", "p2geno", "maxpostprob")
      parents <- pivot_longer(mout$snpdf[, sele_col], cols = c("p1geno", 
                                                               "p2geno"), names_to = "ind", values_to = "geno")
      parents.depth <- apply(parents, 1, function(x) {
        if (grepl("p1", x[7])) 
          return(data.frame(ref = x[2], size = x[3], 
                            check.names = FALSE))
        else if (grepl("p2", x[7])) 
          return(data.frame(ref = x[4], size = x[5], 
                            check.names = FALSE))
        else return(data.frame(ref = NA, size = NA, check.names = FALSE))
      })
      parents.depth <- do.call(rbind, parents.depth)
      parents <- cbind(parents[, c(1, 6:8)], parents.depth)
      parents$ind <- gsub("p1geno", "parent1", parents$ind)
      parents$ind <- gsub("p2geno", "parent2", parents$ind)
    }
    else {
      sele_col <- c("snp", "maxpostprob", "pgeno", "p1ref", 
                    "p1size")
      parents <- data.frame(mout$snpdf[, sele_col[1:2]], 
                            ind = "parent", mout$snpdf[, sele_col[3:5]], 
                            check.names = FALSE)
      colnames(parents)[4:6] <- c("geno", "ref", "size")
    }
    inddf <- mout$inddf[, c(1, 7, 2, 3, 4, 5)]
    inddf <- rbind(parents, inddf)
    inddf$ref <- as.numeric(inddf$ref)
    inddf$size <- as.numeric(inddf$size)
  }
  else {
    inddf <- mout$inddf[, c(1, 7, 2, 3, 4, 5)]
  }
  depth_df <- inddf %>% group_by(snp) %>% summarize(total_ref = sum(ref), 
                                                    total_size = sum(size), total_alt = sum(size) - sum(ref))
  depth_df <- depth_df %>% arrange(match(snp, mout$snpdf$snp))
  new_df <- mout$snpdf %>% tidyr::separate(snp, into = c("CHROM", 
                                                  "POS", "other"), sep = "_") %>% select(CHROM, POS)
  new_df$POS <- sub("^0+", "", new_df$POS)
  vcf_df <- data.frame(CHROM = new_df$CHROM, POS = new_df$POS, 
                       # ID = mout$snpdf$snp,
                       ID = ".",
                       ### NEED TO CHANGE TO REAL  REF ALT ----
                       REF = "A", ALT = "B", QUAL = ".", 
                       FILTER = ".", INFO = NA, FORMAT = NA, check.names = FALSE)
  vcf_df$INFO <- paste0("DP=", depth_df$total_size, ";", "ADS=", 
                        depth_df$total_ref, ",", depth_df$total_alt, ";", "BIAS=", 
                        mout$snpdf$bias, ";", "OD=", mout$snpdf$od, ";", "PMC=", 
                        mout$snpdf$prop_mis)
  vcf_df$FORMAT <- paste("GT", "UD", "DP", "RA", "AD", "MPP", 
                         sep = ":")
  convert_dosage_to_genotype <- function(dosage, ploidy) {
    if (is.na(dosage)) {
      return("./.")
    }
    ref_count <- dosage
    alt_count <- ploidy - dosage
    genotype <- paste(c(rep("0", ref_count), rep("1", alt_count)), 
                      collapse = "/")
    return(genotype)
  }
  geno_df <- inddf[, c("snp", "geno")] %>% mutate(genotype = sapply(geno, 
                                                                    convert_dosage_to_genotype, ploidy = as.numeric(ploidy)))
  format_df <- data.frame(snp = inddf$snp, ind = inddf$ind, 
                          format = paste0(geno_df$genotype, ":", inddf$geno, ":", 
                                          inddf$size, ":", inddf$ref, ":", inddf$ref, ",", 
                                          (inddf$size - inddf$ref), ":", inddf$maxpostprob), 
                          check.names = FALSE)
  format_wide <- format_df %>% pivot_wider(names_from = ind, 
                                           values_from = format)
  vcf_df <- cbind(vcf_df, format_wide[, -1])
  colnames(vcf_df)[1] <- "#CHROM"
  if (!compress) {
    writeLines(vcf_header, con = output.file)
    suppressWarnings(write.table(vcf_df, file = output.file, 
                                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, 
                                 append = TRUE))
  }
  else {
    temp_loc <- tempfile()
    writeLines(vcf_header, con = temp_loc)
    suppressWarnings(write.table(vcf_df, file = temp_loc, 
                                 sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, 
                                 append = TRUE))
    suppressMessages(Rsamtools::bgzip(temp_loc, dest = paste0(output.file, 
                                                              ".gz"), overwrite = FALSE))
  }
}
updog2vcf(mout_cleaned, 
          paste0(outdir, "/updog.genomat.", ploidy_level, ".", chr,".vcf"),
          updog_version = '2.1.5', compress = F)

## Estimate LD
ploidy <- ifelse(ploidy_level == "diploid", 2, 3) 
varnames <- paste0("logL_", 0:ploidy)
larray <- format_multidog(x = mout_cleaned, varname = varnames) # get genotype log-likelihoods

like_ld <- mldest(geno = larray, K = ploidy, type = "comp") # estimate LD

ldmat <- format_lddf(obj = like_ld, element = "r2") # make into matrix

write.table(ldmat, 
            file = paste0(outdir, "/updog.LDmat.", ploidy_level, ".", chr,".txt"), 
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
