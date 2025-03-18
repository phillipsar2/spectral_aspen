# Title: BAMQC Stats for the radseq data
# Author: Alyssa
# Date: 12/5/24

library(stringr)
library(ggplot2)

# Load files and prep
stats <- read.table("/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/stats.bamqc.txt", header = F)
colnames(stats) <- c("bams", "numreads", "GC", "meanMQ", "perreadsmapped")
stats$bams <- str_split(stats$bams, pattern = "/", simplify = T)[,10]
# stats$sitename <- str_split(stats$bams, pattern = "_", simplify = T)[,1]
stats$samplename <- str_split(stats$bams, pattern = ".rg", simplify = T)[,1]

# stats[,c(2:4)] <- lapply(stats[,c(2:4)], function(x) gsub(",", "", x))
stats[,2] <- gsub(",", "", stats[,2])
stats[,3] <- gsub("%", "", stats[,3])
stats[,5] <- gsub("%)", "", stats[,5])
stats[,2:5] <- lapply(stats[,2:5], as.numeric)

head(stats)
str(stats)

# Load metadata and join
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv")
str(meta)

# stats_meta <- merge(stats, meta, by.x = "sitename", by.y = "Site_Code", all = T)

# Overview of data quality

pdf(paste0("/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/bamqc.histograms.",Sys.Date(),".pdf"))

stats %>%
  ggplot(aes(x = numreads/1000000)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9, adjust = 0.8) +
  # xlim(c(0,max(geno_dp, na.rm = T))) +
  theme_bw() +
  labs(x = "Total number of reads (millions)")

stats %>%
  ggplot(aes(x = GC)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Average GC content")
# GC content can be affected by duplicates, bacterial contamination, etc.

stats %>%
  ggplot(aes(x = meanMQ)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Mean mapping quality (MQ)")

stats %>%
  ggplot(aes(x = perreadsmapped)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  theme_bw() +
  labs(x = "Percent reads mapped")

dev.off()

# Comparing the quality statistics
pairs(stats[,c(2:5)], upper.panel = NULL)
cor(stats[,c(2:5)])

'Nothing is really standing out as too correlated. Not useful for figuring out what is happening with GC content.'

# Looking at characters of the submitted data
pdf(paste0("~/aspen/radseq/erincar/quality_check/bamqc.stats_by_ploidy.",Sys.Date(),".pdf"))

stats_meta %>%
  ggplot(aes(x = Ploidy_level, y = GC, fill = Ploidy_level)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "GC content", x = "Ploidy")

stats_meta %>%
  ggplot(aes(x = Ploidy_level, y = perreadsmapped, fill = Ploidy_level)) +
  geom_violin() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill="none") +
  labs(y = "% reads mapped", x = "Ploidy")

dev.off()       


# look at samples for nquack
nquack <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/known_ploidy_genos.txt")
names <- readxl::read_xlsx("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/final_plant-level_spectra_2022-2023_03-14-24.xlsx", sheet = 2)



stats$samplename[stats$samplename %in% names$`Our accession`]
stats$samplenames[stats$samplename %in% names$`SRA Accession`]


nquack_sub <- names$`SRA Accession`[names$ID_genotype_JGI[names$`SRA Accession` %in% stats$samplename] %in% nquack$V1]
temp <- names$`Our accession`[names$ID_genotype_JGI[names$`Our accession` %in% stats$samplename] %in% nquack$V1]

keep_names <- c(nquack_sub, temp) %>% sort()
keep_names <- keep_names[!keep_names == "NA"]
length(keep_names)

nquack_sub <- stats[stats$samplename %in% keep_names,]

write.csv(nquack_sub$samplename, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/test_cytotype_samplename.csv", row.names = F)
