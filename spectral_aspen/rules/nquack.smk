# Determining ploidy with nQuack


# (1) Data preparationi
# -U output regions that don't meet filtering criteria
# -L select regions in repeat bed file
rule rm_rpts:
    input:
        bam = "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam",
        bed = "/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/mex_genome/genome.1MX.modifiedchrnames.fasta.mod.EDTA.intact.gff.bed"
    output:
#        bam = "/global/scratch/users/arphillips/spectral_aspen/data/nquack/filtered/{sample}.filt.bam"
        txt = "/global/scratch/users/arphillips/spectral_aspen/data/nquack/processed/{sample}.txt"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/samtools.yaml"
    shell:
        """
        samtools view {input.bam} -h -U -L {input.bed} | \
        samtools view -b -q 10 | \
        samtools mpileup --no-BAQ --ff UNMAP,DUP -A -Q 0 -q 0 /dev/stdin > {output.txt}
        """
