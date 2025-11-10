rule pixy:
    input:
        vcf = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.{chr}.RMBL.gvcf.gz",
        popfile = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/mar/2025-07-02/{mar}.2025-09-29.populations.txt"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/data/pixy/09292025/{mar}.{chr}.w{winsize}_pi.txt"
    params:
        tmp = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/data/updog/vcf/updog.genomat.{chr}.RMBL.sorted.gvcf.gz",
        winsize = "{winsize}",
        dir = "/global/scratch/users/arphillips/spectral_aspen/data/pixy/09292025",
        pre = "{mar}.{chr}.w{winsize}"
    conda: "/global/scratch/users/arphillips/toolz/.conda/envs/pixy"
    shell:
        """
        mkdir -p {params.dir}
        #gunzip {input.vcf}
        #bgzip {input.vcf}
        #bcftools sort {input.vcf} -Oz -o {params.tmp}
        #tabix {params.tmp}
        pixy --stats pi watterson_theta \
        --vcf {params.tmp} \
        --populations {input.popfile} \
        --window_size {params.winsize} \
        --output_folder {params.dir} \
        --output_prefix {params.pre} \
        --n_cores 1
        """ 
