# builds GTF/BED/minimap2 index (this one has no #SBATCH header — you may want to add one)
awk '{ if ($0 ~ /transcript_id/) print $0; else print $0 " transcript_id \"" NR "\";"; }' newvolume/genomes/Mus_musculus.GRCm39.109.gtf | ./bin/convert2bed --input=gtf --output=bed --do-not-sort > GRCm39.bed

sed 's/^[^\t]*/GRCm39_&/' GRCm39.bed > GRCm39.pre.bed

sed 's/^>/\>GRCm39_/g' newvolume/genomes/Mus_musculus.GRCm39.dna.primary_assembly.fa > GRCm39.pre.fa

pigz -k GRCm39.pre.fa




MMAP2=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/minimap2-2.24_x64-linux/minimap2
FA=/gpfs/Labs/Dovat/bxp5423/Nanopore_ONT_scRNAseq/GRCm39.pre.fa.gz
$MMAP2 -x map-ont -d GRCm39.pre.mmi $FA
