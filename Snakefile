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

# Chromosomes
fai =  pd.read_csv("/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.fai", header = None, sep = "\t")
CHROM = list(fai[0])
#print(CHROM)

# Date
DATE = datetime.datetime.utcnow().strftime("%Y-%m-%d")

# SNP filters
MAX_DP = ["30"]
MIN_DP = ["1", "3"]



# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        ## Calling
#        fastqc = expand("/global/scratch/users/arphillips/spectral_aspen/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
#        fastp = expand("/global/scratch/users/arphillips/spectral_aspen/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
#        bamqc = expand("/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/{sample}_stats/genome_results.txt", sample = SAMPLE),
#        bamqc_stats = "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/stats.bamqc.txt"
#        vcf = expand("/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz", chr = CHROM),
#        table = expand("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/rad_aspen.{chr}.table", chr = CHROM)
        dp = expand("/global/scratch/users/arphillips/spectral_aspen/reports/filtering/depth/rad_aspen.{chr}.filtered.nocall.table", chr = CHROM)
        dp_filt = expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.{min_dp}dp{max_dp}.nocall.vcf", chr = CHROM)
#        filt_vcf = expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.{min_dp}dp{max_dp}.nocall.vcf.gz", min_dp = MIN_DP, max_dp = MAX_DP)
        ## nQuack
        filt = expand("/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/{sample}.txt", sample = SAMPLE)

# =================================================================================================
#     Rule Modules
# =================================================================================================
#include: "rules/mapping.smk"
include: "rules/calling.smk"
#include: "rules/genotyping.smk"
include: "rules/nquack.smk"