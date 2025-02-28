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

for(i in 1:length(samples)){
  temp <- read.csv(paste0(inpath, samples[i], ".csv"))
  ## ADDED LINE HERE
  temp <- temp[which(temp$mixture %in% c("diploid", "triploid")), ]
  summary <- quackit(model_out =  temp, 
                     summary_statistic = "BIC", 
                     mixtures = c("diploid", "triploid"))
  write.csv(summary, 
            file = paste0(outpath, samples[i], ".csv"),
            row.names = FALSE)
}

# Load metadata for evaluation
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/spectral_aspen.subset_ploidy.meta.csv")
meta$sample <- paste0(meta$ID, ".rg")
str(meta)

meta$ploidal.level <- lapply(meta$Ploidy, function(j){
  if (j == 2){
    "diploid"
    } else if (j == 3){
      "triploid"
    } else {"NA"}
  }
)

# Read in quackit() output 
## 98 has issues
dfs <- lapply(list.files(path = outpath, full.names = TRUE  ), read.csv)
alloutput <- do.call(rbind, dfs[1:97])
str(alloutput)

# Combined
alloutputcombo <- dplyr::left_join(alloutput, meta, keep = T) # join by 'sample'
head(alloutputcombo)
tail(alloutputcombo)

# Check the accuracy
alloutputcombo <- alloutputcombo %>%
  dplyr::mutate(accuracy = ifelse(winnerBIC == ploidal.level, 1, 0))

## What distribution and model type should we use?
sumcheck <- alloutputcombo %>% 
  group_by(Distribution, Type) %>% 
  summarize(total = n(), correct = sum(accuracy))

sumcheck_by_ploidy <- alloutputcombo %>% 
  group_by(Distribution, Type, ploidal.level) %>% 
  summarize(total = n(), correct = sum(accuracy))

sumcheck_by_geno <- alloutputcombo %>% 
  group_by(Distribution, Genotype, ploidal.level) %>% 
  summarize(total = n(), correct = sum(accuracy))

dip_genos <- sumcheck_by_geno$Genotype[sumcheck_by_geno$ploidal.level == "diploid"] %>% unique()
alloutputcombo$ID[alloutputcombo$Genotype %in% dip_genos] %>% unique

kbl(sumcheck) %>%
  kable_paper("hover", full_width = F) 

write.csv(sumcheck, paste0(outpath, "model_distribution_accuracy" ,Sys.Date(), ".csv"))

# Are the same genotypes assigned similar ploidy levels?


# Plot accuracy for diploids vs triploids under each model
ggplot(alloutputcombo, aes(x = ploidal.level, y = BICdif, color = ploidal.level)) +
  geom_point()
