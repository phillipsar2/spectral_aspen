### Title: GUSrelate
### Author: Alyssa Phillips
### Date: 1/27/2025

# library("devtools")
# devtools::install_github("tpbilton/GUSbase")
# devtools::install_github("tpbilton/GUSrelate")
library("GUSrelate")
library(vcfR)
library(stringr)
# install.packages("pheatmap")
library(pheatmap)
# install.packages("egg")
library(egg)

# (1) Load data ----
file <- as.character("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.3dp30.nocall.vcf.gz", verbose = FALSE )

vcf <- read.vcfR(file, verbose = FALSE )
vcf

# (2) Filter VCF ----

## LD Thinning
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
  
vcf.thin <- distance_thin(vcf, min.distance = 500)
vcf.thin

## MAF
min_mac <- function(vcfR,
                    min.mac=NULL){
  
  #if vcfR is not class vcfR, fail gracefully
  if (!inherits(vcfR, what = "vcfR")){
    stop("specified vcfR object must be of class 'vcfR'")
  }
  
  #if all input SNPs are not bi-allelic, minor allele count can't be calculated accurately, let user know
  if (max(nchar(gsub(",","",vcfR@fix[,"ALT"])),na.rm=T) > 1){
    stop("Input vcf contains SNPs with > 2 alleles. MAC is calculated under a strict assumption that a single SNP can only possess two alleles. Please use 'filter_biallelic(vcfR)' to remove multi-allelic sites before implementing a MAC filter.")
  }
  
  if (is.null(min.mac)){
    
    #print message
    message("no filtering cutoff provided, vcf will be returned unfiltered")
    
    #convert vcfR to matrix and make numeric
    gt.matrix<-vcfR::extract.gt(vcfR)
    missingness.og<-sum(is.na(gt.matrix)) #store missingness
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/0"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"
    missingness.new<-sum(is.na(gt.matrix)) #store missingness after the conversion
    
    #if unrecognized genotype values were present throw an error
    if (missingness.og != missingness.new){
      stop("Unrecognized genotype values in input vcf. Only allowed non-missing genotype inputs are '0/0','0/1','1/0','1/1'.")
    }
    
    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)
    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }
    
    #hist folded mac with cutoff shown
    graphics::hist(sfs,
                   main="folded SFS",
                   xlab = "MAC")
    
    #return unfiltered vcf
    return(vcfR)
    
  }
  
  else{
    
    #if specified min.mac is not numeric, fail gracefully
    if (!inherits(min.mac, what = "numeric")){
      stop("specified min.mac must be numeric")
    }
    
    #convert vcfR to matrix and make numeric
    gt.matrix<-vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"
    
    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)
    
    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }
    
    #hist folded mac with cutoff shown
    graphics::hist(sfs,
                   main="folded SFS",
                   xlab = "MAC")
    graphics::abline(v=min.mac-1,
                     col="red")
    
    #calculate % of SNPs to be removed
    p<-round((sum(sfs < min.mac)/length(sfs))*100, 2)
    
    #print message to user
    message(p, "% of SNPs fell below a minor allele count of ", min.mac, " and were removed from the VCF")
    
    #filter vcfR
    vcfR <- vcfR[sfs >= min.mac,]
    
    #return vcf object
    return(vcfR)
    
  }
  
}

# vcf.maf <- min_mac(vcf.thin, min.mac = 2)
# vcf.maf

## Save filtered vcf
filt_thin_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/vcf/spectral_aspen.filt.LD.vcf")
write.vcf(vcf.thin, file = filt_thin_file)

# filt_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/obv_aspen/data/vcf/obv_aspen.filt.LD.MAC.vcf")
# write.vcf(vcf.maf, file = filt_file)

# (2) Load in metadata ----
meta <- read.csv("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv")
str(meta)
dim(meta)

# (3) Create required genotype metadata file ----
## Exclude samples with ploidy data
meta_ploidy <- meta[!is.na(meta$Ploidy_level),]

## Split colnames into site codes
site <- str_split(colnames(vcf.thin@gt)[-1], pattern = "_", simplify = T)[,1] %>%
  as.factor()
ID <- colnames(vcf.thin@gt)[-1] %>% 
  as.factor()
site_ID <- data.frame(site, ID)

## Subset to relevant columns and select sites with ploidy
meta_ploidy_sub <- meta_ploidy %>%
  dplyr::select(Site_Code, Ploidy_level, Genotype) %>%
  dplyr::filter(Site_Code %in% site_ID$site)

## Add ID
geno_info <- merge(site_ID, meta_ploidy_sub, by.x = "site", by.y = "Site_Code", all.x = T) %>%
  dplyr::select(ID, Ploidy_level, Genotype) %>%
  dplyr::filter(ID != "HSYDC_448_18")
dim(geno_info)

# ID <- ID[!ID == 'HSYDC_448_18']
# ID <- ID[site %in% meta_ploidy$Site_Code] ## is this order right? YES
# Group <- meta$Genotype[meta$Site_Code %in% site]
# ploidy_level <- meta$Ploidy_level[meta$Site_Code %in% site]

geno_info$Ploidy <- ifelse(geno_info$Ploidy_level == "Diploid", yes = 2, no = 3)
geno_info <- dplyr::select(geno_info, ID, Ploidy, Genotype)
dim(geno_info)
# geno_info <- data.frame(ID, Ploidy, Group) # required column names

geno_meta_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/spectral_aspen.filt.meta.csv")
write.csv(geno_info, geno_meta_file, # ID, Ploidy, Group
          row.names = F,
          quote = F)

## Remove samples without ploidy info
noploidy <- geno_info$ID[is.na(geno_info$Ploidy)]

vcf.sub <- vcf.thin[,!colnames(vcf.thin@gt) %in% noploidy]
vcf.sub

sub_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/vcf/spectral_aspen.subset_ploidy.vcf")
write.vcf(vcf.sub, file = sub_file)

geno_info_sub <- geno_info[!geno_info$ID %in% noploidy,]
dim(geno_info_sub)

geno_sub_file <- as.character("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/spectral_aspen.subset_ploidy.meta.csv")
write.csv(geno_info_sub, geno_sub_file, # ID, Ploidy, Group
          row.names = F,
          quote = F)


# (3) Create GRM object ----
ra <- VCFtoRA(sub_file) # convert VCF to RA format
radata <- readRA(ra) # process RA file to RA object
grm <- makeGRM(RAobj = radata, 
               samfile = geno_sub_file)

# (4) Compute genotype relatedness matrix (GRM) ----
# seq_error = 0.01
grm$computeGRM(name = "GRM_VR", method="VanRaden")

# (5) PCA of GRM ----
grm$PCA(name = "GRM_VR", colour="Ploidy", shape=NULL, )
ggsave(filename = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/grm.PCA.ploidy.",Sys.Date(),".jpeg"),
       width = 4, height = 5, unit = "in")

grm$PCA(name = "GRM_VR", colour="Genotype", shape=NULL,interactive = T) 
  # guides(colour="none")
ggsave(filename = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/grm.PCA.genotype.",Sys.Date(),".jpeg"),
       width = 4, height = 5, unit = "in")

# (6) Plot GRM ----
grm_mat <- grm$extractGRM(name = "GRM_VR")
dim(grm_mat)
heat <- pheatmap(grm_mat, fontsize = 4)
ggsave(heat, filename = paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/grm.heatmap.",Sys.Date(),".jpeg"), 
       height = 20, width = 20, units = "in")

kin <- which(grm_mat > 0.5, arr.ind = TRUE) %>% as.data.frame() # rownames are colnames

kin$colnames <- colnames(grm_mat)[kin$col] 
kin$rownames <- rownames(grm_mat)[kin$row]

kin_noself <- kin[kin$colnames != kin$rownames,]

kin$values <- grm_mat[cbind(kin$row, kin$col)]
head(kin)

# Expected self-relatedness ~ 1-1.25
# Parent-offspring ~ 0.5-0.75
# Full-siblings ~ 0.5
# Unrelated = 0

# (7) Export GRM ----
grm_file <- paste0("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gusrelate/grm.filt.",Sys.Date(),".csv")
grm$writeGRM(name = "GRM_VR", filename = grm_file, IDvar=NULL)
