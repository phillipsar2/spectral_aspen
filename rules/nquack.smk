# Determining ploidy with nQuack


# (1) Data preparation! 
# Running on each bam in paralell,rather than using a for loop
#
# -U output regions that don't meet filtering criteria
# -L select regions in repeat bed file
# Have to skip the repetitive region filtering for now as we don't have an annotation:
# #  samtools view {input.bam} -h -U -L {input.bed} | \
# #      samtools view {input.bam} -b -q 10 | \
#        samtools mpileup --no-BAQ --ff UNMAP,DUP -A -Q 0 -q 0 /dev/stdin > {output.txt}
rule prepare:
    input:
        bam = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam",
    output:
        txt = "/global/scratch/users/arphillips/spectral_aspen/data/nquack/prepared/{sample}.rg.txt"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/nquack.yaml"
    params:
        bam = "{sample}.rg",
        tmp = "/global/scratch/users/arphillips/tmp/nquack_prep/{sample}/"
    shell:
        """
#        mkdir -p {params.tmp}
        Rscript scripts/nquack_02prep.R {params.bam} --tmp {params.tmp}
#        rm -r {params.tmp}
        """
