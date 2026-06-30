#runs minimap2 alignment + stats + bam conversion 
#!/bin/bash
#SBATCH --job-name=alignment_minimap2
#SBATCH --output=alignment_minimap2%j.out
#SBATCH --error=alignment_minimap2%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate spipe
module load parallel
module load fastqc/0.12.1

PROCESS=/gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/10of18/20250320_1510_3B_PAY54190_959d6892/fastq_pass/S1-out/process/
MMAP2=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/minimap2-2.24_x64-linux/minimap2
IDX=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/GRCm39.pre.mmi
BED=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/GRCm39.pre.bed
FQ=/gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/10of18/20250320_1510_3B_PAY54190_959d6892/fastq_pass/S1-out/process/barcode_head.fastq.gz


# Takes approximately 1.5 hours for 60M reads

# ONT parameters
"$MMAP2" --MD -a -u f -x splice -t 16 --junc-bed "$BED" "$IDX" "$FQ" > ${PROCESS}aligned.sam
