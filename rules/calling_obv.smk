# (7) Extract SNPs from each raw vcf and subsample to the correct samples
rule get_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz",
        samples = "/global/scratch/projects/fc_moilab/aphillips/obv_aspen/sample_list.txt" 
    output:
         "/global/scratch/users/arphillips/obv_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bcftools.yaml"
    shell:
        """
        bcftools view {input.vcf} --types snps --samples-file {input.samples} --force-samples -Oz -o {output} 
        """


# (8) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
# Extract alternate base quality (QUAL), mapping quality (MQ), depth at the site (DP), and allele depth for each genotype (AD)
rule diagnostics:
    input:
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/obv_aspen/reports/filtering/rad_aspen.{chr}.table"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F QUAL -F DP -F MQ -GF AD \
        -O {output}
        """

# (9) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz"
    output:
        "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.snps.vcf"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        """
        gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
        -O {output}
        """

# (10) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.snps.vcf"
    output:
        "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        """
        gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}
        """

# (11) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/obv_aspen/reports/filtering/depth/rad_aspen.{chr}.filtered.nocall.table"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output}
        """

# (12) Fitlter by depth of each genotype at each site
rule filter_depth:
    input:
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    output:
        dp = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.6dp30.vcf"
    shell:
        """
        gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < 6 || DP > 30 \" \
        -G-filter-name \"DP_6-30\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}
        """

# (13) Filter snps for genotype missingness
rule depth_nocall:
    input:
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.6dp30.vcf",
    output:
        vcf = "/global/scratch/users/arphillips/obv_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.6dp30.nocall.vcf"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered true --max-nocall-fraction 0.2 -O {output}"

