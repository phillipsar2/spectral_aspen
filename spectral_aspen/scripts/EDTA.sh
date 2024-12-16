#!/bin/sh
#SBATCH --job-name=edta
#SBATCH --account=fc_mel
#SBATCH --partition=savio3
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-00:00:00
#SBATCH --output /global/home/users/arphillips/aspen/spectral_aspen/slurm_log/edta_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/spectral_aspen/slurm_log/edta_%j.err
#SBATCH --chdir /global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/mex_genome/

# Run EDTA

#module load anaconda3
#module load repeatmasker/4.1.3
#module load repeatmodeler/2.0.3
#source activate /global/home/users/arphillips/.conda/envs/EDTA2

# Had to create a ref file with modified chromosome name as the original was too long
ref="genome.1MX.modifiedchrnames.fasta"
cds="Potre.1MX.cds.fasta"
threads=32

# Run EDTA
apptainer run /global/scratch/users/arphillips/toolz/singularity/EDTA.sif EDTA.pl --genome $ref \
--cds $cds \
-t $threads \
--overwrite 1

# Prep for nQuack
cut -f1,4,5 '$ref'.mod.EDTA.intact.gff3 | perl -pi -e 's/^#.*\n//g' > '$ref'.mod.EDTA.intact.gff.bed
