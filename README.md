# long-read-scRNA-seq-pipeline-parse-biosciences

# Long-Read scRNA-seq Pipeline (Parse Biosciences + ONT)

A SLURM-based pipeline for processing single-cell RNA-seq data generated with **Parse Biosciences (Evercode WT)** combinatorial barcoding chemistry, sequenced on **Oxford Nanopore (ONT)** long-read instruments. The pipeline splits long reads into pseudo paired-end format, demultiplexes cell barcodes with `split-pipe`, aligns to a reference genome with `minimap2`, and combines multiple sublibraries into a single analysis-ready dataset.

## Pipeline Overview

```
Raw ONT FASTQ (long reads)
        │
        ▼
1. Split long reads → paired-end chunks (LR_generate_pairs script)
        │
        ▼
2. Build reference genome index (GTF → BED, FASTA prep, minimap2 index)
        │
        ▼
3. Demultiplex cell barcodes (split-pipe --mode pre)
        │
        ▼
4. Align reads to genome (minimap2) + alignment QC stats
        │
        ▼
5. Post-process split-pipe output (split-pipe --mode post)
        │
        ▼
6. Combine sublibraries into final dataset (split-pipe --mode comb)
```

## Requirements

- SLURM workload manager
- [Conda](https://docs.conda.io/) environment named `spipe` with:
  - [`split-pipe`](https://support.parsebiosciences.com/) (Parse Biosciences pipeline, v3 chemistry)
  - Python 3 (for `LR_generate_pairs_1.0.0.py`)
- [`minimap2`](https://github.com/lh3/minimap2) (v2.24 used here)
- [`bedops`](https://bedops.readthedocs.io/) (`convert2bed`) for GTF → BED conversion
- `GNU parallel` (loaded as an environment module)
- `fastqc` (loaded as an environment module)
- `samtools`
- `pigz`

## Repository Structure

```
.
├── scripts/
│   ├── 01_split_chunks_r1_2.sh       # Split ONT long reads into paired-end chunks
│   ├── 02_build_reference.sh         # Build GTF/BED + minimap2 index
│   ├── 03_demultiplexing.sh          # split-pipe --mode pre (barcode demultiplexing)
│   ├── 04_alignment_minimap2.sh      # minimap2 alignment + QC stats + BAM conversion
│   ├── 05_post_processing.sh         # split-pipe --mode post
│   └── 06_combine_all.sh             # split-pipe --mode comb (merge sublibraries)
└── README.md
```

## Usage

> **Note:** Scripts contain hardcoded paths specific to our cluster/project (e.g. `/gpfs/Scratch/DovatLab/...`). Update `INPUT_DIR`, `OUTPUT_DIR`, `MAIN`, `REF`, and similar variables at the top of each script before running on your own data.

### 1. Split long reads into chunks
Splits raw ONT FASTQ files into pseudo paired-end format required by `split-pipe`.
```bash
sbatch scripts/01_split_chunks_r1_2.sh
```

### 2. Build the reference genome
One-time step per reference genome. Converts a GTF to BED format, prepends a prefix to chromosome/contig names, and builds a `minimap2` splice-aware index.
```bash
bash scripts/02_build_reference.sh
```

### 3. Demultiplex cell barcodes
Runs `split-pipe` in `pre` mode to identify and demultiplex Parse Biosciences cell barcodes from the split FASTQ.
```bash
sbatch scripts/03_demultiplexing.sh
```

### 4. Align reads
Aligns demultiplexed reads to the reference using `minimap2` with ONT splice-aware presets (`-x splice`), then generates basic alignment QC stats (unique vs. multimapped reads) and converts SAM → BAM.
```bash
sbatch scripts/04_alignment_minimap2.sh
```

### 5. Post-process
Runs `split-pipe` in `post` mode to generate the cell-by-gene count matrix and QC reports for the sublibrary.
```bash
sbatch scripts/05_post_processing.sh
```

### 6. Combine sublibraries
Merges all processed sublibraries (e.g. 18 flow cells/runs) into a single combined output directory for downstream analysis.
```bash
sbatch scripts/06_combine_all.sh
```

## Notes

- Chemistry version used: `v3` (Parse Biosciences Evercode WT kit)
- Each SLURM job requests 32 CPUs, 128 GB RAM, and a 48-hour walltime; adjust based on your dataset size and cluster availability.
- Alignment of ~60M reads takes approximately 1.5 hours with these settings.
- Email notifications (`--mail-user`) are configured for job start/end/fail; update to your own address.

## Citation / Acknowledgments

This pipeline integrates:
- [Parse Biosciences `split-pipe`](https://support.parsebiosciences.com/)
- [`minimap2`](https://github.com/lh3/minimap2) (Li, H. 2018)
- [`bedops`](https://bedops.readthedocs.io/) (Neph et al. 2012)

## License

Add a license here (e.g. MIT) if you intend others to reuse this pipeline.
