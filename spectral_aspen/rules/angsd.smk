# (1) Create the sites file for ANGSD
rule create_sites:
    input:
#        sites = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.0.99_0.2.txt",
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.{min_dp}dp{max_dp}.nocall.vcf.gz"
    output:
        sites = "/global/scratch/users/arphillips/spectral_aspen/data/angsd/rad_aspen.{min_dp}dp{max_dp}.positions"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bcftools.yaml"
    shell:
        """
        bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT{{0}}\n' {input.vcf} \
        > {output.sites}
        """

# (2) Extract single reads positions that have all individuals present (no missing data)
# sites file must contain the major and minor (ref and alt) allele
rule ibs:
    input:
        ref = config["data"]["reference"]["genome"],
        bamlist = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/2024-11-26.bamlist.txt",
        sites = "/global/scratch/users/arphillips/spectral_aspen/data/angsd/rad_aspen.{min_dp}dp{max_dp}.positions"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/angsd/rad_aspen.{min_dp}dp{max_dp}.ibs.gz"
    params:
        prefix = "/global/scratch/users/arphillips/spectral_aspen/data/angsd/rad_aspen.{min_dp}dp{max_dp}"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/angsd.yaml"
    shell:
        """
        angsd sites index {input.sites}
        angsd \
        -sites {input.sites} \
        -bam {input.bamlist} \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}
        """
