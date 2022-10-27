#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH -t 12:00:00
#SBATCH -J script_spectrogram
#SBATCH -o JOB%j.out # File to which STDOUT will be written
#SBATCH -e JOB%j.out # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.bower@yale.edu

module load R/4.1.0-foss-2020b

Rscript ./script_spectrogram.R

cat JOB${SLURM_JOB_ID}.out | mail -s "$SLURM_JOB_NAME Ended id=$SLURM_JOB_ID" mark.bower@yale.edu
