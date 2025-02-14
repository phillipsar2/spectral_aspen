#devtools::install_github("mgaynor1/nQuack")
#devtools::install("/global/scratch/users/arphillips/toolz/nQuack")
library(nQuack)
library("argparser")

getwd()

# Argument name
ap <- arg_parser("nQuack")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "samp", help = "Text file from step 2")

# parse arguments
argv <- parse_args(ap)

samp <- as.character(argv$samp)
print(samp)

inpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/prepared/"
outpath <- "/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/"

## For loop version
# List files in the inpath and remove their ending
newfilelist <- list.files(path = inpath, pattern = "*.txt" )

# for(i in 1:length(newfilelist)){
for(i in 410:503){
  samp <- newfilelist[i]
  temp <- process_data(paste0(inpath, samp),
                       min.depth = 6,
                       max.depth.quantile.prob = 0.9,
                       error = 0.01,
                       trunc = c(0.1,0.9))


  write.csv(temp,
            file = paste0(outpath, gsub(".txt", "", samp), ".csv"),
            row.names = FALSE)
}

## By sample version
temp <- process_data(paste0(inpath, samp), 
                     min.depth = 6, 
                     max.depth.quantile.prob = 0.9, 
                     error = 0.01, 
                     trunc = c(0.1,0.9))


write.csv(temp, 
          file = paste0(outpath, gsub(".txt", "", samp), ".csv"),
          row.names = FALSE)


# xm <- temp
# hist(xm[,1])
# xi <- xm[,2]/xm[,1]
# hist(xi) # allele frequency plot
