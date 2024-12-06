# (18) Combine vcfs with bcftools
rule combine_vcfs:
    input:
        expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.depth.{min_dp}dp{max_dp}.nocall.vcf", chr = CHROM, min_dp = MIN_DP, max_dp = MAX_DP)
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.{min_dp}dp{max_dp}.nocall.vcf.gz"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bcftools.yaml"
    shell:
        "bcftools concat {input} -Oz -o {output}"
