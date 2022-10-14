#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH -t 12:00:00
#SBATCH -J script_spectrogram
#SBATCH -o log.%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo '-------------------------------'
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'

module load R/4.1.0-foss-2020b

Rscript ./script_spectrogram.R
