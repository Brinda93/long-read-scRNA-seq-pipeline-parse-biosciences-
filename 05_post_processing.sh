#runs split-pipe --mode post

#!/bin/bash
#SBATCH --job-name=Post-Processing_split-pipe
#SBATCH --output=PP_splite_pipe_%j.out
#SBATCH --error=PP_splite_pipe_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-user=bxp5423@psu.edu
#SBATCH --mail-type=ALL

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spipe
module load parallel
module load fastqc/0.12.1

MAIN=/gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/10of18/20250320_1510_3B_PAY54190_959d6892/fastq_pass/
REF=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/newvolume/genomes/GRCm39/


split-pipe \
    --mode post \
    --chemistry v3 \
    --nthreads 16 \
    --genome_dir $REF \
    --parfile ${MAIN}parfile.txt \
    --output_dir ${MAIN}S1-out
