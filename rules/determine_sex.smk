rule toz_cov:
    input:
        bam = 
    output:
    conda: samtools, bedtools
    params:
        subbam = ,
        bed = 
    shell:
        """
        samtools view {input.bam} ""Potre.1MX.sc0049:320490-324352"" > {params.subbam}
        bedtools bamtobed -i {params.subbam} | \
        sort -k1,1 -k2,2n > {params.bed}
        bedtools coverage -d {input.bam}
        """
