# (7) make bam list
rule bamlist:
    input:
       bams = expand("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.sorted.bam", sample = SAMPLE)
    output:
       bamlist = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{date}.bamlist.txt"
    params:
       path = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/"
    shell:
       "ls {params.path} > {output}"

# (8) Call snps initially with bcftools to identify variable sites
# default only sites with max 250 reads considered at each positin, this is way above the max coverage
# -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
# -r output for only the given region
# --annotate FORMAT/AD,FORMAT/DP give allele and genotype depths
rule mpileup:
    input:
        ref = config.ref,
        bamlist = expand("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{date}.bamlist.txt", date = DATE) 
    output:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz"
    params:
        chr = "{chr}"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bcftools.yaml"
    benchmark:
         "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{chr}.mpileup.benchmark.txt"
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
        ref = config.ref,
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz"
    output:
         "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.snps.vcf.gz"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    run:
        """
        gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}
        """


# (10) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
rule diagnostics:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/vcf/rad_aspen.{chr}.raw.vcf.gz",
        ref = config.ref
    output:
        "/global/scratch/users/arphillips/spectral_aspen/reports/filtering/rad_aspen.{chr}.table"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output}
        """

# (11) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config.ref,
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


# (12) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config.ref,
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.filtered.snps.vcf"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/{cov}/rad_aspen.{chr}.filtered.nocall.vcf"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
    shell: 
        """
        gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}
        """
