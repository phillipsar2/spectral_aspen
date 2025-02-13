#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")

getwd()

# Argument name
ap <- arg_parser("nQuack")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "bam", help = "Pre-processed bam")
ap <- add_argument(ap, "--tmp", help = "Path to the temp directory")

# parse arguments
argv <- parse_args(ap)

bam <- as.character(argv$bam)
tmp <- as.character(argv$tmp)
print(bam)

inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/prepared/"

# List files in the inpath and remove their ending
# filelist <- list.files(path = inpath, pattern = "*.rg.bam" )
# filelist <- gsub(".bam", "", filelist)

prepare_data(name = bam, inpath = inpath, outpath = outpath, tempfolder = tmp)
