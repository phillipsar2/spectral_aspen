# (7) make bam list
rule bamlist:
    input:
       bams = expand("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam", sample = SAMPLE)
    output:
       bamlist = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{date}.bamlist.txt"
    params:
       path = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/"
    shell:
       "ls {params.path}*bam > {output}"

# (8) Call snps initially with bcftools to identify variable sites
# default only sites with max 250 reads considered at each positin, this is way above the max coverage
# -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
# -r output for only the given region
# --annotate FORMAT/AD,FORMAT/DP give allele and genotype depths
rule mpileup:
    input:
        ref = config["data"]["reference"]["genome"],
        bamlist = expand("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{date}.bamlist.txt", date = DATE) 
    output:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz"
    params:
        chr = "{chr}"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bcftools.yaml"
#    benchmark:
#         "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{chr}.mpileup.benchmark.txt"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -b {input.bamlist} -r {params.chr} \
        --annotate FORMAT/AD,FORMAT/DP --threads 8 | \
        bcftools call -mv -Oz -o {output.vcf}
        bcftools index -t {output}
        """

# (9) Extract SNPs from each vcf
rule get_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz"
    output:
         "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}
        """


# (10) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
# Extract alternate base quality (QUAL), mapping quality (MQ), depth at the site (DP), and allele depth for each genotype (AD)
rule diagnostics:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/spectral_aspen/reports/filtering/rad_aspen.{chr}.table"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F QUAL -F DP -F MQ -GF AD \
        -O {output}
        """

# (11) Concatenate files for analysis
#  head -n1 rad_aspen.Potre.1MX.Chr01.table > rad_aspen.all.table
# grep -v QUAL  <(cat rad_aspen.Potre.1MX.*.table) >> rad_aspen.all.table


# (12) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.snps.vcf"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 30.0\" --filter-name \"MQ30\" \
        -O {output}
        """

# (13) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.snps.vcf"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}
        """

# (14) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/spectral_aspen/reports/filtering/depth/rad_aspen.{chr}.filtered.nocall.table"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output}
        """

# (15) Concatenate files for analysis
#  head -n1 rad_aspen.Potre.1MX.Chr01.filtered.nocall.table > rad_aspen.all.filtered.nocall.table
# grep -v CHROM  <(cat rad_aspen.Potre.1MX.*.filtered.nocall.table) >> rad_aspen.all.filtered.nocall.table


# (16) Fitlter by depth of each genotype at each site
# 6 < DP < 30 with 10% missing is too strict - only 15 SNPs
rule filter_depth:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    output:
        dp = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.{min_dp}dp{max_dp}.vcf"
    params:
        min = "{min_dp}",
        max = "{max_dp}"
    shell:
        """
        gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < {params.min} || DP > {params.max} \" \
        -G-filter-name \"DP_{params.min}-{params.max}\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}
        """

# (17) Filter snps for genotype missingness
rule depth_nocall:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.{min_dp}dp{max_dp}.vcf",
    output:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.{min_dp}dp{max_dp}.nocall.vcf",
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered true --max-nocall-fraction 0.1 -O {output}"
