### Title: Aspen MAR
### Author: Alyssa Phillips
### Date: 2/18/25

# library(devtools)
# install_github("meixilin/mar")
library(mar)
library(stringr)
library(dplyr)


# Prepare genotype data ----
## no missing data, bialleleic
# geno_dips <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.2025-01-27.txt", header = T)
# geno_trips <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.2025-01-27.txt", header = T)
# 
# geno_trips[1:5,1:5]

## Swap homozygous ref to 0
### Triploids
# geno_trips_flip <- apply(geno_trips, c(1,2), function(x) {
#   if (x == 3)
#     print("0") else
#   if (x == 2)
#     print("1") else
#   if (x == 1)
#     print("2") else
#   if(x == 0)
#     print("3")
#   }
#   )

# geno_trips[1:15,1:5]
# geno_trips_flip[1:15,1:5]

### Diploids
# geno_dips_flip <- apply(geno_dips, c(1,2), function(x) {
#       if (x == 2)
#         print("0") else
#           if (x == 1)
#             print("1") else
#               if(x == 0)
#                 print("2")
# }
# )

## Merge
# genos <- merge(geno_dips_flip, geno_trips_flip, by = 0, all = TRUE)
# dim(genos)
# 
# genos_comp <- genos[complete.cases(genos),]
# dim(genos_comp)
# genos_comp[1:5,1:5]

## Thin SNPs by one every 500 bp
# pos <- str_split(genos_comp$Row.names, pattern = "_", simplify = T) %>%
#   as.data.frame()
# colnames(pos) <- c("chr", "pos")
# 
# pos$int <- cut(as.numeric(pos[,2]), 
#     seq( from = min(as.numeric(pos[,2])), 
#          to = max(as.numeric(pos[,2])), 
#          by = 500 )
#     )
# 
# keep_pos <- pos %>%
#   group_by(chr, int) %>% # group by the chr and interval, 4743 intervals
#   slice_sample(n = 1) %>% # randomly select one
#   mutate(chr_pos = paste0(chr, "_", pos)) # remerge names
#   
# genos_filt <- filter(genos_comp, Row.names %in% keep_pos$chr_pos)
# paste0("Final number of SNPs: ",dim(genos_filt)[1])

## Write to file without rownames or colnames. Save rownames and colnames to other files.
# write.table(genos_filt[,-1],
          # "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/mar.genomat.all.thinned.tsv",
          # row.names = F, col.names = F, sep = "\t")
# write.table(colnames(genos_filt[,-1]),
          # "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/mar.genomat.colnames.tsv",
          # row.names = F, col.names = F, sep = "\t")
# write.table(genos_filt[,1],
#             "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/mar.genomat.pos.tsv",
#             row.names = F, col.names = F, sep = "\t")

# Prepare Geographic data ----
geno_names <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/mar.genomat.colnames.tsv")
site_names <- str_split(geno_names$V1, pattern = "_", simplify =T)[,1]

meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv")
str(meta)

## subset metadata to those in genotype file
meta_sub <- meta %>%
  filter(Site_Code %in% site_names) %>%
  select(Site_Code, Latitude, Longitude, Ploidy_level) # pick appropriate columns

## order metadata to match genotype file
meta_sub$Site_Code == site_names

meta_sub_order <- meta_sub[match(site_names, meta_sub$Site_Code),] 

meta_sub_order$Site_Code == site_names

## Change ploidy to integers
ploidy <- ifelse(meta_sub$Ploidy_level == "Diploid", yes = 2, no = 3)

## export
write.table(meta_sub_order[,1:3],
            "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/mar.latlon.tsv",
            row.names = F, col.names = T)
write.table(ploidy,
            "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/mar.ploidy.tsv",
            row.names = F, col.names = F)

# Run MAR pipeline ----
geno_file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/mar.genomat.all.thinned.tsv"
geo_file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/mar.latlon.tsv"
dir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/"
ploidy_file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/mar.ploidy.tsv"

MARPIPELINE(name = "RMBL_aspen", 
            workdir = dir, 
            genofile = geno_file, 
            lonlatfile = geo_file, 
            saveobj = TRUE,
            option_geno = list(ploidy = 3) )
