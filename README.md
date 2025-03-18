# Assessing the genetic diversity of RMBL aspen

Author: Alyssa Phillips

Snakemake pipeline pipeline for calling variants and estimating diversity across RMBL quaking aspen. This repository follows a CookieCutter directory structure.

Data source: [https://doi.org/10.26078/JNPN-4A13](https://doi.org/10.26078/JNPN-4A13) & SRA BioProject PRJNA1063753

## Running the pipeline

`rm -r .snakemake/log  .snakemake/metadata .snakemake/slurm_logs`
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
|    ├── backup_data	 <- Final VCFs and BAM files 
│    └── processed	 <- Intermediate VCF files in filtering
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
* Map reads to the reference with bwa-mem2. RMBL reads were mapped as pair-end (-p) and all other reads were mapped as single-end reads.
* Sort and add read groups to BAM files with samtools and GATK.
* Assess mapping quality with Qualimap's bamqc.

3. Ploidy determination
* Ploidy was previously estimated for a majority of samples using [gbs2ploidy](https://doi.org/10.32614/CRAN.package.gbs2ploidy)
* Ploidy will be validated using [nQuack](https://github.com/mgaynor1/nQuack).
* For nQuack, BAMs are filtered for MQ > 10 and repetitive regions annotated by EDTA are excluded. 

4. Variant calling
* BCFtools is used to call variants.
	RMBL raw SNPs: 3,942,970
        Spectral raw SNPs: 1,504,711
* Variants are seperately filtered for RMBL spatial dataset & spectral dataset. Variants were hard filtered to keep QUAL > 30, MQ > 30, and biallelic sites
        RMBL SNPs after hard filters: 1,413,730
        Spectral raw SNPs: 739,735
* Multiple depth and genotype missingness filters were tested.
 	RMBL 3 < DP < 30 & 10%: 26,475 
	RMBL 1 < DP < 30 & 10%: 111,953
        Spectral 3 < DP < 20 & 10%:
        Spectral 1 < DP < 20 & 10%: 
* Updog (https://dcgerard.github.io/updog/index.html) will be used to estimate genotypes. This software considers sequencing error, allele bias, and overdispersion. The input of updog is read counts for SNPs, so variats must be called and initially filtered prior to genotype calling.
* After genotyping, an LD filter will be applied 
