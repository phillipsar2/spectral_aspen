# Assessing the genetic diversity of RMBL aspen

Author: Alyssa Phillips

Snakemake pipeline pipeline for calling variants and estimating diversity across RMBL quaking aspen. This repository follows a CookieCutter directory structure.

Data source:

## Running the pipeline

`snakemake --executor slurm --profile profiles/ --use-conda`

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
└── config.yaml  
</pre>
