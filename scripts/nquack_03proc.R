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
# for(i in 410:503){
#   samp <- newfilelist[i]
#   temp <- process_data(paste0(inpath, samp),
#                        min.depth = 6,
#                        max.depth.quantile.prob = 0.9,
#                        error = 0.01,
#                        trunc = c(0.1,0.9))
# 
# 
#   write.csv(temp,
#             file = paste0(outpath, gsub(".txt", "", samp), ".csv"),
#             row.names = FALSE)
# }

## By sample version
temp <- process_data(paste0(inpath, samp, ".txt"), 
                     min.depth = 6, 
                     max.depth.quantile.prob = 0.9, 
                     error = 0.01, 
                     trunc = c(0.1,0.9))


write.csv(temp, 
          file = paste0(outpath, gsub(".txt", "", samp), ".csv"),
          row.names = FALSE)


## Examine output
## Comment out when not using 

i = 1
# xm <- read.csv(paste0(outpath, gsub(".txt", "", newfilelist[i]), ".csv"))
xm <- read.csv(paste0(outpath, gsub(".txt", "", 	"SRR27496155.rg"), ".csv"))

# Plot coverage
hist(xm[,1], xlab = "Coverage", main = NULL, breaks = 15, xlim = c(5,25))
# abline(v = 6, col = "red")
# abline(v = quantile(xm[,1], probs=0.9, na.rm=TRUE), col = "red")

# Error cutoffs
### If I increase the sequence error rate, how many sites will likely be removed?
# new.e <- 0.02 # 2 sites out of every 100 
# removes <- c()
# for(i in 1:nrow(xm)){
#   if(xm[i,2] < (xm[i,1]*new.e) | xm[i,2] > (xm[i,1]*(1-new.e))){
#     removes[i] <- 1
#   }else{
#     removes[i] <- 0
#   }
# }
# sum(removes)

# Allele frequency plot
xi <- xm[,2]/xm[,1]
hist(xi, xlab = "Allele frequency", main = "SRR27496155.rg", breaks = 100)

# plot(x = xm[,1], y = xi, xlab = "coverage", ylab = "allele freq")
