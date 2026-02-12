# Determine ploidy
rule gbs2ploidy:
    input:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/backup_data/vcf/RMBL_aspen.nocall.3dp30.vcf.gz"
#        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.1dp30.per0.25.vcf.gz"
    output:
         "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/{spec_samp}.propOut.csv"
    params:
        geno = "{spec_samp}",
        tmp_dir = "/global/scratch/users/arphillips/spectral_aspen/temp/gbs2ploidy/{spec_samp}",
        temp = "/global/scratch/users/arphillips/spectral_aspen/temp/gbs2ploidy/{spec_samp}/{spec_samp}.vcf.gz"
    conda: "/global/home/users/arphillips/.conda/envs/stuff_in_r"
    shell:
        """
        mkdir -p {params.tmp_dir}
        bcftools view -Oz -s {params.geno} {input} >> {params.temp}
        Rscript scripts/gbs2ploidy.R {params.temp} --out {output}
        rm -rf {params.tmp_dir}
        """

# Subset SNPs for updog
rule get_snps_fromgvcf:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.gvcf"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.snps.vcf.gz"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk SelectVariants --java-options ""-Xmx4G"" \
        -R {input.ref} \
        -V {input.vcf} \
        --select-type-to-include SNP \
        -O {output}
        """

# Subset invariant sites for later
rule get_invar:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.gvcf"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.invar.vcf.gz"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk SelectVariants --java-options ""-Xmx4G"" \
        -R {input.ref} \
        -V {input.vcf} \
        --select-type-to-include NO_VARIATION \
        -O {output}
        """

# Additional filters applied that exclude poor quality genotypes.
rule updog_dips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.snps.vcf.gz"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.{chr}.{dataset}.txt" 
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "diploid",
        chr = "{chr}"
    conda: "stuff_in_r"
    shell:
        "Rscript scripts/updog.R {input.vcf} --ploidy {params.ploidy} --cores 4 --outdir {params.outdir} --chr {params.chr}"

rule updog_trips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.snps.vcf.gz"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.{chr}.{dataset}.txt"
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "triploid",
        chr = "{chr}"
    conda: "stuff_in_r"
    shell:
        "Rscript scripts/updog.R {input.vcf} --ploidy {params.ploidy} --cores 8 --outdir {params.outdir} --chr {params.chr}"

# Merge diploids and triploid vcfs output from updog
rule updog_merge:
    input:
        dip = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.{chr}.vcf.gz",
        trip = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.{chr}.vcf.gz"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.{chr}.vcf.gz"
    shell:
        """
        bcftools index -f {input.dip}
        bcftools index -f {input.trip}
        bcftools merge {input.dip} {input.trip} -Oz -o {output}
        """

# Concatonate SNPs and invariant sites
rule comb_snp_inv:
    input:
        snp = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.{chr}.vcf.gz",
        invar = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.invar.vcf.gz",
        keep = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/gbs2ploidy/RMBL.passedsamples.2025-04-17.csv"
    params:
        tmp = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.{dataset}.filtered.invar.ploidy.vcf.gz"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.{chr}.{dataset}.gvcf.gz"
    shell:
        """
        bcftools index -f {input.invar}
        bcftools view -S {input.keep} --force-samples {input.invar} -Oz -o {params.tmp}
        bcftools concat -n {params.tmp} {input.snp} > {output}
        rm {params.tmp}
        """
