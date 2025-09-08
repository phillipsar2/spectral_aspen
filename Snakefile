configfile: "profiles/config.yaml"
import pandas as pd
from random import randint
import datetime

# =================================================================================================
#    Wildcards
# =================================================================================================

# Sample names
sample_file = pd.read_csv("sequences_to_work_with.11132024.csv", header = 0)
SAMPLE = list(sample_file.filenames)
#print(SAMPLE)

## Diploid sample names (subset)
DIPLOIDS = ["ABEY_034_18","AQXKL_473_18","AWUSD_479_18","BFPMS_412_18","BJSRB01_488_18","BS_503_18","BYFM_035_18","CAJTV_468_18","CAKME_401_18","CCRWS01_372_18","CCTOR01_417_18","CFRID_474_18","CVHQZ_314_18","DDKJ_261_18","DEGLR_421_18","DRCAD01_374_18", "DRCRS01_378_18","DVNBJ_470_18"]

## Known ploidy sample names
ploidy_samp = pd.read_csv("metadata/known_flow_with_accnum.csv")
PLOID_SAMP = list(ploidy_samp["accession"])

## Additional SRA seq
srr_file = pd.read_csv("additional_sequences.txt", header = None)
SRR = list(srr_file[0])
#print(SRR) 

## Spectra genotypes
spectra_sampl = pd.read_csv("metadata/spectra_genotypes.csv")
SPEC_SAMP = list(spectra_sampl["Accession"])

# Chromosomes
fai =  pd.read_csv("/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai", header = None, sep = "\t")
CHROM = list(fai[0])
#print(CHROM)

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# SNP filters
## Depth
MAX_DP = ["30"]
MIN_DP = ["1"]

## Missing data per SNP (10%, 20%, and 25%) 
MISS = ["0.25", "0.2", "0.1"] 

## Which dataset I'm working with ["SWUS", or "RMBL"]
DATASET = ["RMBL"]


# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        ## Calling
#        ncbi = expand("/global/scratch/users/arphillips/spectral_aspen/raw/{srr}.fastq.gz", srr = SRR),
#        bamqc = expand("/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/{srr}_stats/genome_results.txt", srr = SPEC_SAMP),
#        bamqc = expand("/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE),
#        bamqc_stats = "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/stats.bamqc.txt"
#        vcf = expand("/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.{dataset}.raw.gvcf.gz", chr = CHROM, dataset = DATASET),
        vcftools_filt = expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.gvcf.gz", chr = CHROM, dataset = DATASET)
#        table = expand("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/rad_aspen.{chr}.table", chr = CHROM),
#        dp = expand("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/depth/rad_aspen.{chr}.filtered.nocall.table", chr = CHROM),
#        dp_filt = expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.nocall.{min_dp}dp{max_dp}.per{miss}.vcf.gz", chr = CHROM, min_dp = MIN_DP, max_dp = MAX_DP, miss = MISS),
#        comb_filt = expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.{min_dp}dp{max_dp}.per{miss}.vcf.gz", min_dp = MIN_DP, max_dp = MAX_DP, miss = MISS)
#        comb_raw =  "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/raw/rad_aspen.all.raw.vcf.gz"
        ## nQuack
#        prep = expand("/global/scratch/users/arphillips/spectral_aspen/data/nquack/prepared/{sample}.rg.txt", sample = PLOID_SAMP),
#        infer = expand("/global/scratch/users/arphillips/spectral_aspen/data/nquack/model_inference/{sample}.rg.csv", sample = PLOID_SAMP)
#        boot = expand("/global/scratch/users/arphillips/spectral_aspen/data/nquack/bootstrap/{sample}.rg-boots.csv", sample = PLOID_SAMP)
        ## Genotyping
#        gbs2ploidy = expand("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/{spec_samp}.propOut.csv", spec_samp = SAMPLE)
#        gbs2ploidy = expand("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/{spec_samp}.propOut.csv", spec_samp = SPEC_SAMP)
#        updog_dip = expand("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.{date}.txt", date = DATE),
#        updog_trip = expand("/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.{date}.txt", date = DATE),
        ## Subset for Obv study
#        merge_raw = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/raw/rad_aspen.all.raw.vcf.gz",
#        table = expand("/global/scratch/users/arphillips/obv_aspen/reports/filtering/rad_aspen.{chr}.table", chr = CHROM)
#        dp = expand("/global/scratch/users/arphillips/obv_aspen/reports/filtering/depth/rad_aspen.{chr}.filtered.nocall.table", chr = CHROM)
#         dp_filt = expand("/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.6dp30.nocall.vcf", chr = CHROM)

# =================================================================================================
#     Rule Modules
# =================================================================================================
#include: "rules/mapping.smk"
include: "rules/calling.smk"
#include: "rules/calling_obv.smk"
#include: "rules/genotyping.smk"
#include: "rules/nquack.smk"
#include: "rules/angsd.smk"
