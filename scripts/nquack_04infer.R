#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")
library(dplyr)

getwd()

# Argument name
ap <- arg_parser("nQuack")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "samp", help = "CSV file from step 3")

# parse arguments
argv <- parse_args(ap)

samp <- as.character(argv$samp)
print(samp)

inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/model_inference/"

# Sample List
# samples <- c("MLG013", "MLG014", "MLG015")
#samples <- list.files(path = inpath, pattern = "*.csv" ) 
#samples <- lapply(samples, function(x) gsub(x = x, pattern = ".csv", replacement ="")) %>% unlist()
#samp_sub <- samples[1:10]


# Run all models for subset of samples (100?)
#for(i in 1:length(samp_sub)){
#  temp <- as.matrix(read.csv(paste0(inpath, samples[i], ".csv")))
#  out1 <- quackNormal(xm = temp, samplename = samples[i], cores = 10, parallel = FALSE)
#  out2 <- quackBeta(xm = temp, samplename = samples[i], cores = 10, parallel = FALSE)
#  out3 <- quackBetaBinom(xm = temp, samplename = samples[i], cores = 10, parallel = FALSE)
#    allout <- rbind(out1, out2, out3)
#  write.csv(allout, 
#            file = paste0(outpath, samples[i], ".csv"),
#            row.names = FALSE)
#}

### One at a time
temp <- as.matrix(read.csv(paste0(inpath, samp, ".csv")))
out1 <- quackNormal(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
out2 <- quackBeta(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
out3 <- quackBetaBinom(xm = temp, samplename = samp, cores = 10, parallel = FALSE)
allout <- rbind(out1, out2, out3)

write.csv(allout, 
             file = paste0(outpath, samp, ".csv"),
             row.names = FALSE)
