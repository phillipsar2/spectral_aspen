#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")
library(dplyr)
library(kableExtra)
library(ggplot2)

inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/model_inference/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/model_interpret/"

# samples <- c("MLG013", "MLG014", "MLG015")
samples <- list.files(path = inpath, pattern = "*.csv" )
samples <- lapply(samples, function(x) gsub(x = x, pattern = ".csv", replacement ="")) %>% unlist()

known_samp <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/known_flow_with_accnum.csv")
known_samp$file <- paste0(known_samp$accession, ".rg")
dim(known_samp)

samples_sub <- samples[samples %in% known_samp$file]

# Select models based on the BIC scores only considering 2x, 3x, and 4x models
# for(i in 1:length(samples)){
#   temp <- read.csv(paste0(inpath, samples[i], ".csv"))
#   summary <- quackit(model_out =  temp, 
#                      summary_statistic = "BIC", 
#                      mixtures = c("diploid", "triploid", "tetraploid"))
#   write.csv(summary, 
#             file = paste0(outpath, samples[i], ".csv"),
#             row.names = FALSE)
# }

for(i in 1:length(samples_sub)){
  temp <- read.csv(paste0(inpath, samples_sub[i], ".csv"))
  ## ADDED LINE HERE
  temp <- temp[which(temp$mixture %in% c("diploid", "triploid")), ]
  summary <- quackit(model_out =  temp, 
                     summary_statistic = "BIC", 
                     mixtures = c("diploid", "triploid"))
  write.csv(summary, 
            file = paste0(outpath, samples_sub[i], ".csv"),
            row.names = FALSE)
}

# Load metadata for evaluation
## gbs2ploidy
# meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/spectral_aspen.subset_ploidy.meta.csv")
# meta$sample <- paste0(meta$sample, ".rg")
# meta$ploidal.level <- lapply(meta$Ploidy, function(j){
#   if (j == 2){
#     "diploid"
#   } else if (j == 3){
#     "triploid"
#   } else {"NA"}
# }
# )

## flow samples
# names_key <- readxl::read_xlsx("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/final_plant-level_spectra_2022-2023_03-14-24.xlsx", sheet = 2)
# flow_data <- readxl::read_xlsx("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/test_cytotype_data.xlsx")

flow_data <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/metadata/known_flow_with_accnum.csv")
flow_data$file <- paste0(flow_data$accession, ".rg")

# meta <- merge(x = flow_data, y = names_key, by = "ID_genotype_JGI")

# meta_sub <- rbind(meta[!meta$`Our accession` == 'NA',], meta[!meta$`SRA Accession` == 'NA',])
# meta_sub <- meta_sub %>%
#   arrange('Our accession', 'SRA Accession') %>%
#   select(ID_genotype_JGI, flow_ploidy, 'SRA Accession', 'Our accession') 
# meta_sub$sample <- c(meta_sub$`Our accession`[1:5], meta_sub$`SRA Accession`[6:45])
# meta_sub$sample <- paste0(meta_sub$sample, ".rg")
# str(meta_sub)


# Read in quackit() output 
# dfs <- lapply(list.files(path = outpath, full.names = TRUE  ), read.csv)
dfs <- lapply(paste0(outpath, samples_sub, ".csv"), read.csv)
# alloutput <- do.call(rbind, dfs)
# alloutput <- bind_rows(dfs) %>% # binds even if different dimensions
#   select(-X, -total, -correct)
alloutput <- bind_rows(dfs)

# Combined
# alloutputcombo <- dplyr::left_join(alloutput, meta, keep = T) # join by 'sample'
# alloutputcombo <- merge(x = alloutput, y= known_samp, by.x='sample', by.y = 'file')
alloutputcombo <- merge(x = alloutput, y = flow_data, by.x='sample', by.y = 'file', all.y = T)
head(alloutputcombo)
tail(alloutputcombo)
dim(alloutputcombo)

alloutputcombo_sub <- alloutputcombo[complete.cases(alloutputcombo),]
dim(alloutputcombo_sub)
head(alloutputcombo_sub)

# Check the accuracy
alloutputcombo_sub <- alloutputcombo_sub %>%
  dplyr::mutate(accuracy = ifelse(winnerBIC == flow_ploidy, 1, 0))

## What distribution and model type should we use?
sumcheck <- alloutputcombo_sub %>% 
  group_by(Distribution, Type) %>% 
  summarize(total = n(), correct = sum(accuracy))

sumcheck_by_ploidy <- alloutputcombo_sub %>% 
  group_by(Distribution, Type, flow_ploidy) %>% 
  # filter(flow_ploidy == 'diploid') %>%
  summarize(total = n(), correct = sum(accuracy))

alloutputcombo_sub %>% 
  group_by(Distribution, Type, flow_ploidy) %>% 
  filter(Distribution == 'normal-uniform') %>%
  summarize(total = n(), correct = sum(accuracy)) 
  
sumcheck_by_sample <- alloutputcombo_sub %>%
  # group_by(Distribution, Type, sample, flow_ploidy) %>%
  group_by(sample, flow_ploidy) %>%
  summarize(total = n(), correct = sum(accuracy))

# dip_genos <- sumcheck_by_geno$Genotype[sumcheck_by_geno$ploidal.level == "diploid"] %>% unique()
# alloutputcombo$ID[alloutputcombo$Genotype %in% dip_genos] %>% unique

kbl(sumcheck) %>%
  kable_paper("hover", full_width = F) 

write.csv(sumcheck, paste0(outpath, "model_distribution_accuracy.flow_samples." ,Sys.Date(), ".csv"))
write.csv(sumcheck_by_ploidy, paste0(outpath, "model_distribution_accuracy_by_ploidy.flow_samples." ,Sys.Date(), ".csv"))
write.csv(sumcheck_by_sample, paste0(outpath, "model_distribution_accuracy_by_sample.flow_samples." ,Sys.Date(), ".csv"))

# Are the same genotypes assigned similar ploidy levels?


# Plot accuracy for diploids vs triploids under each model
sumcheck_by_ploidy$distribution_type  <- paste0(sumcheck_by_ploidy$Distribution, "_", sumcheck_by_ploidy$Type)
alloutputcombo_sub %>% 
  mutate(model = paste0(Distribution, "_",Type)) %>% 
  group_by(model, flow_ploidy) %>%
  summarize(total = n(), correct = sum(accuracy)) %>%
  ggplot( aes(x = model, y = correct, fill = flow_ploidy)) +
  geom_bar(stat="identity", width = 0.5) +
  # scale_fill_viridis(discrete=TRUE, name="") +
  # theme_ipsum() +
  ylab("Number correct") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept = 41) +
  ylim(c(0,45)) +
  guides()
