#!/bin/sh
#SBATCH --job-name=edta
#SBATCH --account=fc_mel
#SBATCH --partition=savio3
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-00:00:00
#SBATCH --output /global/home/users/arphillips/aspen/spectral_aspen/slurm_log/edta_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/spectral_aspen/slurm_log/edta_%j.err
#SBATCH --chdir /global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/mex_genome

# Run EDTA

module load anaconda3
module load repeatmasker/4.1.3
module load repeatmodeler/2.0.3
source activate /global/home/users/arphillips/.conda/envs/EDTA2

ref="genome.1MX.fasta"
cds="Potre.1MX.cds.fasta"
threads=32

perl /global/home/users/arphillips/toolz/EDTA/EDTA.pl --genome $ref \
--cds $cds \
-t $threads \
--overwrite 1
