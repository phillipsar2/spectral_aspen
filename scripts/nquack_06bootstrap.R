#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")
library(dplyr)

# Argument name
ap <- arg_parser("nQuack")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "samp", help = "CSV file from step 3")

# parse arguments
argv <- parse_args(ap)

samp <- as.character(argv$samp)
print(samp)

inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/bootstrap/"


bout <- c()

temp <- as.matrix(read.csv(paste0(inpath, samp, ".csv")))
# beta-uniform, fixed
bout <- quackNboots(temp, 
		    nboots = 100,
                    distribution = "beta",
                    type = "fixed",
                    uniform = 1, 
                    mixtures = c("diploid", "triploid"), 
                    samplename = samp)

write.csv(bout, file = paste0(outpath, samp, "-boots.csv"), row.names = FALSE)

##
## Analyze results
## 

# bootpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/bootstrap/"
# 
# samples <- list.files(path = bootpath, pattern = "*.csv" )
# samples <- lapply(samples, function(x) gsub(x = x, pattern = "-boots.csv", replacement ="")) %>% unlist()
# 
# # Leave commented out unless using
# 
# # Row 1 = best model for original dataset
# # Row 2 = bootstrap replicates best model
# boots <- lapply(list.files(path = bootpath, full.names = TRUE  ), read.csv)
# allboots <- do.call(rbind, boots)
# str(allboots)
