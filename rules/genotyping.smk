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

# Additional filters applied that exclude poor quality genotypes.
rule updog_dips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.3dp30.nocall.vcf.gz",
#        meta = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv" 
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.diploid.{date}.txt" 
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "diploid"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/updog.yaml"
    shell:
        "Rscript scripts/updog.R {input.vcf} --ploidy {params.ploidy} --cores 4 --outdir {params.outdir}"

rule updog_trips:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.all.depth.3dp30.nocall.vcf.gz",
#        meta = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv"
    output:
        "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/updog.genomat.triploid.{date}.txt"
    params:
        outdir = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog",
        ploidy = "triploid"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/updog.yaml"
    shell:
        "Rscript scripts/updog.R {input.vcf} --ploidy {params.ploidy} --cores 8 --outdir {params.outdir}"
