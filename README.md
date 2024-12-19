# Assessing the genetic diversity of RMBL aspen

Author: Alyssa Phillips

Snakemake pipeline pipeline for calling variants and estimating diversity across RMBL quaking aspen. This repository follows a CookieCutter directory structure.

Data source:

## Running the pipeline

`snakemake --executor slurm --profile profiles/ --use-conda --rerun-triggers input`

## Project organization
<pre>
├── README.md  
|  
├── rules  
|    ├── calling.smk  
|    ├── determine_ploidy.smk  
|    ├── filtering.smk  
|    └── mapping.smk  
|  
├─  environment.yml  
├─  scripts  
│    ├── README.md  
│    └── filtering       <- Custom scripts for variant filtering  
|  
├── data  
│    ├── raw 		 <- Original data dump  
│    ├── trimmed         <- Trimmed fastqs
│    ├── genome 	 <- Reference genome  
│    ├── interim  	 <- Intermediate files in read mapping and SNP calling  
│    └── processed	 <- Final vcfs for analysis  
|  
├── reports 		 <- Generated analyses as HTML, PDF, or .txt.  
├── qc 			 <- Quality check output for raw data  
├── Snakefile  
└── profiles/config.yaml  
</pre>

1. Pre-processing of reads
* Assess read quality with fastqc
* Trim reads with fastp and re-evaluate quality. Reads are trimmed via sliding windows (4 bp windows, min quality of 15) and automated detection of adapters.

2. Mapping
* Map reads to the reference with bwa-mem2
* Sort and add read groups to BAM files with samtools and GATK.
* Assess mapping quality with qualimap's bamqc (read duplication metric is incorrect)

3. Ploidy determination
* Ploidy was previously estimated for a majority of samples using [gbs2ploidy](https://doi.org/10.32614/CRAN.package.gbs2ploidy)
* Ploidy will be validated using [nQuack](https://github.com/mgaynor1/nQuack).
* For nQuack, BAMs are filtered for MQ > 10 and repetitive regions annotated by EDTA are excluded. 

4. Variant calling
* BCFtools is used to call variants.
	Raw SNPs: 3,825,061
* Variants were hard filtered to keep QUAL > 30, MQ > 30, and biallelic sites
        SNPs after hard filters:
* Multiple depth and genotype missingness filters were tested.
 	3 < DP < 30 & 10%: 
	1 < DP < 30 & 10%:
* Updog (https://dcgerard.github.io/updog/index.html) will be used to estimate genotypes. This software considers sequencing error, allele bias, and overdispersion. The input of updog is read counts for SNPs, so variats must be called and initially filtered prior to genotype calling.
* After genotyping, an LD filter will be applied 
