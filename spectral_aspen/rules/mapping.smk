# (1) Evaluate quality of raw reads with fastqc
# scripts/fastqc.sh

rule fastqc:
    input:
        "/global/scratch/projects/fc_moilab/PROJECTS/aspen/moi_aspen/fastq/samples/{sample}.fastq.gz"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/qc/fastqc/{sample}_fastqc.zip"
    params:
        tmp = "/global/scratch/users/arphillips/tmp/spectral_aspen/fastqc/{sample}",
        outdir = "/global/scratch/users/arphillips/spectral_aspen/qc/fastqc"
    conda:
        "/global/home/users/arphillips/aspen/spectral_aspen/envs/fastqc.yaml"
    resources:
        cpus_per_task=2
    shell:
        """
        mkdir -p {params.tmp}
        fastqc -o {params.outdir} -t {resources.cpus_per_task} -d {params.tmp} -f fastq {input}
        rm -rf {params.tmp}
        """

# (2) Trim reads sequenced
# threads (-w)
# Minimum length is 36 (-l 36)
# Don't filter for quality (-Q)
# Adapter trimming is enabled by default -- don't need to specify adapter seq
# Default detects and trims polyG tails for NovaSeq data
# --cut_front is sliding window trimming from 5' ot 3'
rule fastp_trim:
    input:
        fastq = "/global/scratch/projects/fc_moilab/PROJECTS/aspen/moi_aspen/fastq/samples/{sample}.fastq.gz"
    output:
        trim = temp("/global/scratch/users/arphillips/spectral_aspen/data/trimmed/{sample}.trim.fastq.gz")
    conda:
        "/global/home/users/arphillips/aspen/spectral_aspen/envs/fastp.yaml"
#    benchmark:
#        "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{sample}.trim.benchmark.txt"
    shell:
        """
        fastp -w 1 \
        -l 36 -Q \
        -i {input.fastq} \
        -o {output.trim} \
        --cut_front --cut_front_window_size 4 --cut_front_mean_quality 15 
        """

# (3a) Prepare reference file
# 9.6 GB per core available (core == threads)
rule bwa_prep:
    input: 
        config["data"]["reference"]["genome"]
    output:
        index = "/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bwa-mem2.yaml"
    shell:
        """
        ~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index {input}
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input}
        """


# (3b) Align reads to the reference genome
# Paired-end reads in single file
rule bwa_map:
    input:
        ref = config["data"]["reference"]["genome"],
        index = "/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta.0123",
        trim = "/global/scratch/users/arphillips/spectral_aspen/data/trimmed/{sample}.trim.fastq.gz"
    output:
        temp("/global/scratch/users/arphillips/spectral_aspen/data/interm/mapped_bam/{sample}.mapped.bam")
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/bwa_map.yaml"
#    benchmark:
#         "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{sample}.bwa.benchmark.txt"
    shell:
        "~/toolz/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 1 -p {input.ref} {input.trim} |"
        "samtools view -Sb > {output}"

# (4) Sort bams
rule samtools_sort:
    input:
        "/global/scratch/users/arphillips/spectral_aspen/data/interm/mapped_bam/{sample}.mapped.bam"
    output:
        temp("/global/scratch/users/arphillips/spectral_aspen/data/interm/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/global/scratch/users/arphillips/tmp/spectral_aspen/sort_bam/{sample}"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/samtools.yaml"
#    benchmark:
#       "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{sample}.sort.benchmark.txt"
    shell:
        """ 
        mkdir -p {params.tmp}
        samtools sort -T {params.tmp} -@ 1 {input} > {output}
        rm -rf {params.tmp}
        """

# (5) Add read groups
rule add_rg:
    input:
        "/global/scratch/users/arphillips/spectral_aspen/data/interm/sorted_bam/{sample}.sorted.bam"
    output:
        bam = touch("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam")
    params:
        tmp = "/global/scratch/users/arphillips/tmp/spectral_aspen/addrg/{sample}",
        sample = "{sample}",
        rg = randint(1,1000)
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/gatk.yaml"
#    benchmark:
#       "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{sample}.add_rg.benchmark.txt"
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -RGID {params.rg} \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX true
        rm -rf {params.tmp}
        """

# (6) Assess alignment quality metrics with qualimap
# nr is normally 100000 and -nt is normally 8, java mem size = 48
# nw is normally 400
# for higher cov, make nr 1000 and -nt 12, java mem size = 64
rule bamqc:
    input:
        "/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam"
    output:
        "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/{sample}_stats/qualimapReport.html"
    params:
        dir = "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/{sample}_stats",
        stats_dir = "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc"
    conda: "/global/home/users/arphillips/aspen/spectral_aspen/envs/qualimap.yaml"
#    benchmark:
#         "/global/scratch/users/arphillips/spectral_aspen/benchmarks/{sample}.bamqc.benchmark.txt"
    shell:
        """
        qualimap bamqc \
        -bam {input} \
        -nt 12 \
        -nr 1000 \
        -nw 400 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=24G
        grep -h "bam file =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/bamlist.bamqc.txt
        grep -h "number of reads =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/numreads.bamqc.txt
        grep -h "number of duplicated reads (flagged) =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/numdups.bamqc.txt
        grep -h "median insert size =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/medianinsertsize.bamqc.txt
        grep -h "GC percentage =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/GC.bamqc.txt
        grep -h "mean coverageData =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/meancoverage.bamqc.txt
        grep -h "mean mapping quality =" {params.dir}/genome_results.txt | cut -d"=" -f2 >> {params.stats_dir}/meanMQ.bamqc.txt
        grep -h "number of mapped reads =" {params.dir}/genome_results.txt | cut -d "(" -f2 >> {params.stats_dir}/perreadsmapped.bamqc.txt 
        """

rule bamqc_stats:
    input:
        expand("/global/scratch/users/arphillips/spectral_aspen/data/interm/addrg/{sample}.rg.bam", sample = SAMPLE)
    output:
        "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc/stats.bamqc.txt"
    params:
        stats_dir = "/global/scratch/users/arphillips/spectral_aspen/reports/bamqc"
    shell:
        """
        paste {params.stats_dir}/bamlist.bamqc.txt \
        {params.stats_dir}/numreads.bamqc.txt \
        {params.stats_dir}/numdups.bamqc.txt \
        {params.stats_dir}/medianinsertsize.bamqc.txt \
        {params.stats_dir}/GC.bamqc.txt \
        {params.stats_dir}/meancoverage.bamqc.txt \
        {params.stats_dir}/meanMQ.bamqc.txt \
        {params.stats_dir}/perreadsmapped.bamqc.txt | \
        column -s $'\t' -t >> {output}
        """
