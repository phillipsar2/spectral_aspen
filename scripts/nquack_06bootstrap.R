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
