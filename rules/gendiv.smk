rule pixy:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.RMBL.filtered.gvcf.gz",
        popfile =
    output:
    params:
        winsize = "{winsize}",
        dir = "/global/scratch/users/arphillips/spectral_aspen/data/pixy" 
    shell:
        """
        mkdir -p {params.dir}
        pixy --stats pi watterson_theta \
        --vcf {input.vcf} \
        --populations {input.popfile} \
        --window_size {params.winsize} \
        --output_folder {params.dir} \
        --output_prefix {params.winsize} \
        --n_cores 4
        """ 
