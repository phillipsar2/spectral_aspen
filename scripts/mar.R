### Title: Aspen MAR
### Author: Alyssa Phillips
### Date: 2/18/25

# library(devtools)
# install_github("meixilin/mar")
library(mar)
library(stringr)
library(dplyr)
library(vcfR) #v1.15


# Prepare genotype data ----
vcf <- read.vcfR("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz",
                 verbose = FALSE)

## Filter for minor allele frequency
maf_est <- maf(vcf, element = 2) # 2 is the minor allele
hist(maf_est[,4])
sum(maf_est[,4] > 0.05)
sum(maf_est[,4] > 0.01)

# is.na(vcf@gt[,-1][maf_est[,4] <= 0.05,]) <- TRUE # mask low frequency variants
# vcf

vcf_maf <- vcf
vcf_maf@gt <- vcf@gt[maf_est[,4] > 0.01,]
vcf_maf@fix <- vcf@fix[maf_est[,4] > 0.01,] # need to subset both matrices
paste0("Number of SNPs after maf filter: ", dim(vcf_maf@gt)[1])

## Thin VCF to 1 SNP pre 500 bp
distance_thin <- function(vcfR,
                          min.distance=NULL){
  
  #if vcfR is not class vcfR, fail gracefully
  if (!inherits(vcfR, what = "vcfR")){
    stop("specified vcfR object must be of class 'vcfR'")
  }
  
  #function is useless if no distance is provided
  if (is.null(min.distance)){
    stop("filtering distance must be provided")
  }
  
  #logical test specifying that the minimum distance between SNPs for filtering must be at least 1 base pair or the logic doesn't work
  if (is.numeric(min.distance) != TRUE){
    stop("specified filtering distance must be numeric")
  }
  
  #logical test specifying that the minimum distance between SNPs for filtering must be at least 1 base pair or the logic doesn't work
  if (min.distance < 1){
    stop("filtering distance must be >= 1 base pair")
  }
  
  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[1] != "CHROM"){
    stop("vcfR incorrectly formatted. vcfR@fix column 1 must be 'CHROM'")
  }
  
  #logical test to ensure formatting of input vcf will allow accurate analysis of position in genome
  if (colnames(vcfR@fix)[2] != "POS"){
    stop("vcfR incorrectly formatted. vcfR@fix column 2 must be 'POS'")
  }
  
  #set min distance specified by user
  j=min.distance
  
  #generate dataframe containing information for chromosome and bp locality of each SNP
  df<-as.data.frame(vcfR@fix[,1:2])
  
  #write test to identify and remove duplicated SNPs in input vcf
  if (length(unique(paste(df$CHROM,df$POS))) < nrow(df)){
    #remove duplicated SNPs
    vcfR<-vcfR[!duplicated(paste(df$CHROM,df$POS)),]
    #report to user
    message(nrow(df) - length(unique(paste(df$CHROM,df$POS)))," duplicated SNPs removed from input vcf")
    #regenerate df without duplicate inputs
    df<-as.data.frame(vcfR@fix[,1:2])
  }
  
  #generate list of all unique chromosomes in alphabetical order
  chroms<-levels(as.factor(df$CHROM))
  
  #intialize empty df to hold filtering
  keep.df<-data.frame()
  
  #make progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(chroms), style = 3)
  
  #begin tracker
  pbtrack<-1
  
  #loop over each chromosome
  #for t in vector containing the name of each chromosome
  for (t in chroms){
    
    #isolate the SNP positions on the given chromosome
    fix.sub<-as.numeric(df$POS[df$CHROM == t])
    
    #order the positions numerically
    fix.sub<-fix.sub[order(fix.sub)]
    
    #set the first position
    prev<-fix.sub[1]
    
    #initialize empty vector
    k<-c()
    
    #always keep first SNP on the chromosome
    k[1]<-TRUE
    
    #loop to decide whether to keep each following SNP
    if (length(fix.sub) < 2){
      #if chrom only has 1 SNP, do nothing
    } else{
      
      #else, use a for loop to determine which SNPs to keep that satisfy our distance criteria
      for (i in 2:length(fix.sub)){
        
        #store logical indicating whether this SNP is greater than j base pairs from the previous SNP
        k[i]<- fix.sub[i] > prev+j
        
        #if it is, then we keep this SNP, making it the new 'previous' for assessing the next point.
        #If we don't keep the SNP, we don't update the closest point
        if (fix.sub[i] > prev+j){
          prev<-fix.sub[i]
        }
        
        #close for loop
      }
      #close else statement
    }
    
    #make a dataframe with the precise info for each SNP for this chromosome
    chrom.df<-data.frame(CHROM=rep(t, times=length(fix.sub)), POS=fix.sub, keep=k)
    
    #now we rbind in the information for this chromosome to the overall df
    keep.df<-rbind(keep.df,chrom.df)
    
    #empty df for this chrom to prepare for the next one
    chrom.df<-NULL
    
    #update progress bar
    utils::setTxtProgressBar(pb, pbtrack)
    
    #update tracker
    pbtrack<-pbtrack+1
    
  } #close for loop, start over on next chromosome
  
  #close progress bar
  close(pb)
  
  #order the dataframe to match the order of the input vcf file
  #remove scientific notation
  keep.df$POS<-format(keep.df$POS,scientific = FALSE)
  df$POS<-format(df$POS,scientific = FALSE)
  #make sure class matches between the columns you're trying to merge
  keep.df$POS<-as.numeric(as.character(keep.df$POS))
  df$POS<-as.numeric(as.character(df$POS))
  #make sure class matches between the columns you're trying to merge
  keep.df$CHROM<-as.character(keep.df$CHROM)
  df$CHROM<-as.character(df$CHROM)
  #add tracking column
  df$id<-c(1:nrow(df))
  #merge
  order.df<-merge(keep.df,df)
  #order based on tracking column
  order.df<-order.df[order(order.df$id),]
  
  #order.df<-keep.df[match(paste(df$CHROM,format(df$POS,scientific = FALSE)), paste(keep.df$CHROM,format(keep.df$POS,scientific = FALSE))),]
  #note: the position vector must be stripped of scientific notation otherwise identical numbers will not be recognized as matches, giving NA values
  #note: this old approach using match() has been deprecated because there were too many formatting issues causing
  #match to not be able to correctly match the order between 'df' and 'keep.df'. Now we use merge() which seems to be more robust
  
  #write a test to catch if this internal dataset is not able to merge correctly
  if (sum(is.na(order.df)) > .5){
    stop("internal error with the merge function. Please email a copy of your input vcf to devonderaad@gmail.com for a bug fix")
  }
  
  #write a test to catch if this internal dataset is not able to merge correctly
  if (sum(order.df$id != c(1:nrow(df))) > .5){
    stop("internal error with the merge function. Please email a copy of your input vcf to devonderaad@gmail.com for a bug fix")
  }
  
  #subset vcfR locus info based on the logical column from our dataframe
  #vcfR@fix<-vcfR@fix[order.df$keep,]
  #subset genotypes based on logical
  #vcfR@gt<-vcfR@gt[order.df$keep,]
  
  #realized there is no need to do this subsetting separately
  vcfR<-vcfR[order.df$keep,]
  
  #calculate number of total SNPs input
  z<-nrow(keep.df)
  
  #calculate total SNPs retained
  p<-nrow(vcfR@fix)
  
  #print info to screen
  message(p," out of ",z," input SNPs were not located within ",j," base-pairs of another SNP and were retained despite filtering")
  
  #return vcfR
  return(vcfR)
}

vcf.thin <- distance_thin(vcf_maf, min.distance = 500)
vcf.thin

## Subset to genotypes with ploidy data
ploidy_meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL.ploidycalls.2025-04-17.csv")
keep_col <- colnames(vcf.thin@gt) %in% ploidy_meta$sample
keep_col[1] <- TRUE # FORMAT column
vcf.thin@gt <- vcf.thin@gt[,keep_col]
dim(vcf.thin@gt)[2] == dim(ploidy_meta)[1] + 1

## Save VCF
write.vcf(vcf.thin, "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/filtered.RMBL.vcf")

## Extract genotypes
gt <- extract.gt(vcf.thin)
gt[1:5,1:5]

## Swap genotypes to 0,1,2
gt_cont <- NULL
gt_cont <- apply(gt, c(1,2), function(x) {
      if (is.na(x)){      # is.na needs to be the first test
        "0"               # assign to homozygous ref and remove the first alleles
        } else if (x == "0/0"){
          "0" 
          } else if (x == "0/1"){
            "1"
            } else if (x == "1/1"){
              "2"}})

gt[1:5,1:5]
gt_cont[1:5,1:5]
dim(gt_cont)

# Filter out invariant sites (~ 27)
## should accumulate slower on M x A plots
invar <- apply(gt_cont, c(2), as.numeric) %>%
  rowSums() 
sum(invar == 0)
gt_sub <- gt_cont[invar != 0,]
dim(gt_sub)

## Write to file without rownames or colnames. Save rownames and colnames to other files.
write.table(gt_sub,
            "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.genomat.all.thinned.tsv",
            row.names = F, col.names = F, sep = "\t", quote = F)
write.table(colnames(gt_sub),
            "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.genomat.colnames.tsv",
            row.names = F, col.names = F, sep = "\t")
write.table(gt_sub[,1],
            "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.genomat.pos.tsv",
            row.names = F, col.names = F, sep = "\t")

# Prepare Geographic data ----
geno_names <- read.table("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.genomat.colnames.tsv")
site_names <- str_split(geno_names$V1, pattern = "_", simplify =T)[,1]

meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv")
str(meta)

## subset metadata to those in genotype file
meta_sub <- meta %>%
  filter(Site_Code %in% site_names) %>%
  dplyr::select(Site_Code, Longitude, Latitude) # pick appropriate columns

## order metadata to match genotype file
meta_sub$Site_Code == site_names

# meta_sub_order <- meta_sub[match(site_names, meta_sub$Site_Code),] 

meta_sub$genotype <- geno_names$V1


## subset to necessary columns
meta_sub_order <- dplyr::select(meta_sub, genotype, Longitude, Latitude)
colnames(gt_cont) == meta_sub_order$genotype
write.csv(meta_sub_order, 
          file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.latlon.withnames.csv",
          row.names = F, quote = F)

colnames(meta_sub_order) <- c("ID", "LON", "LAT")
dim(meta_sub_order)
head(meta_sub_order)

## set ID to id number
meta_sub_order$ID <- seq(from = 1, to = length(meta_sub_order$ID), by = 1)

## export
write.csv(meta_sub_order, 
          file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.latlon.csv",
          row.names = F, quote = F)

# Run MAR pipeline ----
geno_file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.genomat.all.thinned.tsv"
geo_file = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.latlon.csv"
dir = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/", Sys.Date())

dir.create(dir)
MARPIPELINE(name = "RMBL_aspen", 
            workdir = dir, 
            genofile = geno_file, 
            lonlatfile = geo_file, 
            # option_marext = list(scheme = .MARsampling_schemes, nrep = 10, xfrac = 0.01, quorum = TRUE, animate = TRUE, myseed = NULL),
            saveobj = TRUE)

# geno_file_in <- read.table(geno_file)
# geo_file_in <- read.csv(geo_file)
# geno_file_in[1:5,1:5]
# geo_file_in$ID


# Get MAR output ----
library(stringr)
library(dplyr)

dir = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/")
list.files(dir)

load(paste0(dir, "gm_RMBL_aspen.rda" )) # genotype information
head(gm)

# load(paste0(dir, "extlist_RMBL_aspen.rda" ))
# head(extlist)


# Get sampling grids for Erin out of gm ----
## Load sampling pattern
load(paste0(dir, "mardflist_RMBL_aspen.rda"))
str(mardflist)

## make coordinate files with samples ----
samplingtype <- c("random","inwards","outwards","southnorth","northsouth")

rowcol_sample <- function(mm, bbox, revbbox = FALSE) {
  cellsnotna <- mar:::.rowcol_cellid(mm, bbox, revbbox)
  samples <- mar:::.cellid_sample(mm, cellsnotna)
  return(samples)
}

geno_names <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/mar.latlon.withnames.csv")
head(geno_names)

for (s in samplingtype){
  ## Extract sampling grids
  coor_list <- lapply(mardflist[[s]][,8], function(x){ 
    str_split(x, pattern = ";", simplify = T) # split by ';' into vectors in a list
  })
  coor_num <- lapply(coor_list, as.numeric)
  coor_df <- lapply(coor_num, function(x){ 
    mar:::.rowcol_extent(gm$maps, bbox = unlist(x)) # convert to lat long
  })
  coor_mat <- lapply(coor_df, function(x) { 
    c(x@xmin, x@xmax, x@ymin, x@ymax) }) # format into list
  
  coor_mat <- matrix(unlist(coor_mat), ncol = 4, byrow = TRUE) # convert to matrix
  colnames(coor_mat) <- c("xmin", "xmax", "ymin", "ymax")
  # head(coor_mat)
  
  # Get individuals within each grid 
  ## mm = gm$maps
  ## bbox = grid coordinates
  samples <- lapply(coor_num, function(x){rowcol_sample(gm$maps, bbox = unlist(x))})
  
  ## Pair geno IDs to sample names
  sample_names <- lapply(samples, function(x){ paste( geno_names$ID[unlist(x)], collapse = ',') })
  # str_split(names_stacked, pattern = ",") # how to undo
  
  mar_ouput_df <- data.frame(coor_mat,
                             trees_in_grid = unlist(sample_names))
  # head(mar_ouput_df)
  
  ## Write coordinates to file
  write.table(mar_ouput_df, paste0(dir, s,".", Sys.Date(),".mar_output_df.tsv"), sep = "\t", row.names = F)
}


# Create population file for pixy ----
output_files <- list.files(path = dir, pattern = "*.mar_output_df.tsv")
output_files

for (i in 1:length(output_files)){
  output_df <- read.table(paste0(dir, output_files[i]), header = T)
  output_list <- lapply(output_df$trees_in_grid, function(x){str_split(x,pattern = ",", simplify = T)})
  names(output_list) <- 1:length(output_list)
  
  pop_df <- data.frame(sample = unlist(output_list),
                       pop = paste0("pop",rep( names(output_list), lapply(output_list, length)))
  ) %>%
    arrange(pop) %>%
    filter(sample != "")
  # head(pop_df)
  # unique(pop_df[,1])
  colnames(pop_df) <- NULL # remove header
  
  sampling <- stringr::str_split(output_files[i], pattern = ".mar", simplify = T)[1,1]
  write.table(pop_df, paste0(dir, sampling, ".populations", ".txt"), quote = F, row.names = F, sep = "\t")
}




#################################

# Manually create bed file with window intervals
# sizes <- c(10, 50) # 10k, 50k
# fai <- read.table("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai")
# 
# bed <- c()
# for (s in sizes){
#   for (i in 1:dim(fai)[1]){
#     int_start <- seq(from = 1, to = fai$V2[i], by = s*1000)
#     int_end <- c( (int_start[-1] - 1),  fai$V2[i] )
#     chr_id <- rep(fai$V1[i], length(int_start) )
#     bed_chr <- cbind(chr_id, int_start, int_end)
#     bed <- rbind(bed, bed_chr)
#   }
#   write.table(bed, paste0(dir, "populations.", Sys.Date(), ".txt"), quote = F, row.names = F)
# }


