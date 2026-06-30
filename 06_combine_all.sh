# runs split-pipe --mode comb

#!/bin/bash
#SBATCH --job-name=combine_all
#SBATCH --output=combine_all%j.out
#SBATCH --error=combine_all%j.err
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

split-pipe \
  --mode comb \
  --sublibraries \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/10of18/20250320_1510_3B_PAY54190_959d6892/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/11of18/20250320_1510_3D_PAY57422_5f985bf4/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/12of18/20250320_1510_3F_PAY54219_3f157c43/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/13of18/20250320_1510_3H_PAY54157_791424d6/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/14of18/20250320_1526_2B_PAY54152_9bab7d9b/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/15of18/20250320_1526_2C_PAY92136_7349ee5c/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/16of18/20250320_1526_2G_PAY54192_be60dbfe/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/17of18/20250320_1526_3D_PBA39721_7cfe665e/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/18of18/20250320_1526_3F_PBA54480_ac1ddb96/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/1of18/20250319_1348_2C_PAY88348_f3b27a2c/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/2of18/20250319_1348_2E_PAY87960_a1c75287/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/3of18/20250320_1417_1A_PAY91968_0215e244/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/4of18/20250320_1430_3B_PBA54423_d977081c/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/5of18/20250320_1510_1B_PAY92859_24d7ae01/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/6of18/20250320_1510_1D_PAY91919_a024ebc2/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/7of18/20250320_1510_1E_PAW82933_d5f43f7b/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/8of18/20250320_1510_1F_PAY91861_0ee821ab/fastq_pass/S1-out \
    /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/9of18/20250320_1510_1H_PAW80554_554f628a/fastq_pass/S1-out \
  --output_dir /gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/combined_analysis
