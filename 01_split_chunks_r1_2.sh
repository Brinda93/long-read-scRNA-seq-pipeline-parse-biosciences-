 #splits long reads into paired-end chunks
 #!/bin/bash
#SBATCH --job-name=split_chunks_r1_2
#SBATCH --output=fastq_r1_2_%j.out
#SBATCH --error=fastq_r1_2_%j.err
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

SCRIPT="/gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/10of18/20250320_1510_3B_PAY54190_959d6892/fastq_pass/LR_generate_pairs_1.0.0.py"  # Replace with actual path to the script
INPUT_DIR="/gpfs/Scratch/DovatLab/Derm_Single_Cell/Unaligned/8of18/20250320_1510_1F_PAY91861_0ee821ab/fastq_pass"
OUTPUT_DIR="$INPUT_DIR"  # Same directory for output
CHEMISTRY="v3"
NCPU=16

#split long reads to chunks of paired-end
# Find all fastq.gz files in input folder matching your chunk pattern
FASTQ_FILES=($INPUT_DIR/PAY91861_pass_0ee821ab_b9e3bbf7_*.fastq.gz)

# Run LR_generate_pairs_1.0.0.py on all chunk files in parallel
parallel -j $NCPU python $SCRIPT \
  --out_dir $OUTPUT_DIR \
  --chemistry $CHEMISTRY \
  --fastq {} \
  --multiple_fq \
  --new_fname split \
  --l1dist 1 \
  --l2dist 1 ::: "${FASTQ_FILES[@]}"
