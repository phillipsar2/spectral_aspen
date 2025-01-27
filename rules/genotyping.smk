# (18) Combine vcfs with bcftools
rule combine_vcfs:
    input:
        expand("/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.nocall.{min_dp}dp{max_dp}.vcf", chr = CHROM, min_dp = MIN_DP, max_dp = MAX_DP)
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.{min_dp}dp{max_dp}.nocall.vcf.gz"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bcftools.yaml"
    shell:
        "bcftools concat {input} -Oz -o {output}"


# (19) Genotype with updog
# Additional filters applied that exclude poor quality genotypes.
rule updog_dips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.3dp30.nocall.vcf.gz",
        meta = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv" 
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.{date}.txt" 
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "diploid"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/updog.yaml"
    shell:
        "Rscript scripts/updog.R {input.vcf} --meta {input.meta} --ploidy {params.ploidy} --cores 4 --outdir {params.outdir}"

rule updog_trips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.3dp30.nocall.vcf.gz",
        meta = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.{date}.txt"
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "triploid"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/updog.yaml"
    shell:
        "Rscript scripts/updog.R {input.vcf} --meta {input.meta} --ploidy {params.ploidy} --cores 8 --outdir {params.outdir}"
