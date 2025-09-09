rule pixy:
    input:
        vcf = "/global/scratch/users/arphillips/spectral_aspen/data/processed/filtered_snps/rad_aspen.{chr}.RMBL.filtered.gvcf.gz",
        popfile = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/{mar}.populations.2025-09-08.txt"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/pixy/{mar}.{chr}.w{winsize}_pi.txt"
    params:
        winsize = "{winsize}",
        dir = "/global/scratch/users/arphillips/spectral_aspen/data/pixy",
        pre = "{mar}.{chr}.w{winsize}"
    conda: "/global/scratch/users/arphillips/toolz/.conda/envs/pixy"
    shell:
        """
        mkdir -p {params.dir}
        #gunzip {input.vcf}
        #bgzip {input.vcf}
        #tabix {input.vcf}
        pixy --stats pi watterson_theta \
        --vcf {input.vcf} \
        --populations {input.popfile} \
        --window_size {params.winsize} \
        --output_folder {params.dir} \
        --output_prefix {params.pre} \
        --n_cores 1
        """ 
