#!/bin/bash

# For each MinION run, this script will run a Nextflow workflow which:
#   - Runs Guppy basecalling in sup mode
#   - Concatenates all FASTQ files into 1 file
#   - Concatenates the "sequencing_summary" files and run PycoQC for read QC
#   - Removes organelle-derived reads
# The working dir is /fs/project/PAS2380/assembly/jelmer

# General inputs and variables
genome_id=sus145
ref_fa_ncbi=data/ref/gmax/ncbi/GCF_000004515.6_Glycine_max_v4.0_genomic.fna
organel_contigs="NC_020455.1,NC_007942.1"

# ==============================================================================
#                         PREP SEQS FROM RUN 2022-10-05
#       This run was interrupted - MinKNOW reported an error and it stopped
# ==============================================================================
# Define settings and input files
run_id=GM145_2022-10-05
fast5_dir=../minion/TR145_2022-10-05/fast5_pass
guppy_config=run/config/dna_r9.4.1_450bps_sup.cfg
pyco_minq=9
outdir=results/"$genome_id"/nf_ontreadprep/"$run_id"

# Run the ontreadprep workflow
sbatch mcic-scripts/ont/nf_ontreadprep.sh \
    --fast5_dir "$fast5_dir" \
    --outdir "$outdir" \
    --guppy_config "$guppy_config" \
    --ref_assembly "$ref_fa_ncbi" \
    --contig_blacklist "$organel_contigs" \
    --more_args "--pyco_minq $pyco_minq"

# ==============================================================================
#                       PREP SEQS FROM RUN 2022-12-15
# ==============================================================================
# Define settings and input files
run_id=GM145_2022-12-15
fast5_dir=../minion/TR145_2022-12-15/fast5
guppy_config=run/config/dna_r9.4.1_450bps_sup.cfg
pyco_minq=9
outdir=results/"$genome_id"/nf_ontreadprep/"$run_id"

# Run the ontreadprep workflow
sbatch mcic-scripts/ont/nf_ontreadprep.sh \
    --fast5_dir "$fast5_dir" \
    --outdir "$outdir" \
    --guppy_config "$guppy_config" \
    --ref_assembly "$ref_fa_ncbi" \
    --contig_blacklist "$organel_contigs" \
    --more_args "--pyco_minq $pyco_minq"

# ==============================================================================
#                       PREP SEQS FROM RUN 2023-02-16
#               Flow cell type: FLO-MIN114, Kit type: SQK-LSK114
# ==============================================================================
# Find out default basecall method (e.g., to know the 'bps')
# json=../minion/TR145_2023-02-16/report_FAW32537_20230216_1139_779bff98.json
# tr , "\n" < "$json" | grep "default basecall"

# Define settings and input files
run_id=GM145_2023-02-16
fast5_dir=../minion/TR145_2023-02-16/fast5
guppy_config=run/config/dna_r10.4.1_e8.2_400bps_sup.cfg
pyco_minq=10
outdir=results/"$genome_id"/nf_ontreadprep/"$run_id"

# Run the ontreadprep workflow
sbatch mcic-scripts/ont/nf_ontreadprep.sh \
    --fast5_dir "$fast5_dir" \
    --outdir "$outdir" \
    --guppy_config "$guppy_config" \
    --ref_assembly "$ref_fa_ncbi" \
    --contig_blacklist "$organel_contigs" \
    --more_args "--pyco_minq $pyco_minq"

#! LOST READS
#../minion/TR145_2023-02-16/fast5/FAW32537_779bff98_c64957a2_1157.fast5
#../minion/TR145_2023-02-16/fast5/FAW32537_779bff98_c64957a2_1752.fast5

# ==============================================================================
#                               SUMMARY
# ==============================================================================

# Run1
#> All/pass/1k+ reads: 2.56/2.54/1.93 M // All/pass/1k+ bases: 12.37/12.32/11.95 Gbp
#> Pass reads length: median=4,031 // 90%=10,436
#> Quality: median=13.27
#> Organel removal: Nr of input / output sequences in FASTQ files: 2,541,653  / 2,280,804

# Run2
#> All/pass/1k+ reads: 3.32/2.83/1.46 M // All/pass/1k+ bases: 14.63/12.92/12.42 Gbp
#> Pass reads length: median=2,201 / 90%=11,440
#> Quality: median=13.44
#> Organel removal: Nr of input / output sequences in FASTQ files: 2,831,581  /  2,478,929

# Run3
#> All/pass/1k+ reads: 8.38/7.32/5.16 M // All/pass/1k+ bases: 27.47/24.15/22.8 Gbp
#> Pass reads length: median=1,870 / 90%=8,786
#> Quality: median=15.40
#> Organel removal: Nr of input / output sequences in FASTQ files: 7,320,390 / 6,687,355

# ==============================================================================
#                       CONCATENATE FINAL FASTQ FILES
# ==============================================================================

##TODO - THE BELOW IS JUST COPIED FROM SUS44
fq_concat_all=results/sus44/rm_organel_reads/"$genome_id"/"$genome_id".fastq.gz
fa_concat_all=results/sus44/rm_organel_reads/"$genome_id"/"$genome_id".fasta

# Concatenate
#fq_filt=results/sus145/rm_organel_reads/"$genome_id"/"$run_id".fastq.gz
sbatch mcic-scripts/utils/fqcat.sh -i results/sus44/rm_organel_reads/"$genome_id" -o "$fq_concat_all"
# Number of sequences in output file: ...

# Create a FASTA file
sbatch -A PAS2380 --wrap="
    ml python && source activate /fs/ess/PAS0471/jelmer/conda/seqtk
    seqtk seq -a $fq_concat_all > $fa_concat_all
    ls -lh $fa_concat_all"

# Stats
sbatch -A PAS2380 -c 2 -t 120 --wrap="
    ml python && source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    seqkit stats --all $fq_concat_all"
#>file           format  type   num_seqs         sum_len  min_len  avg_len  max_len     Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)  GC(%)
#>GM44.fastq.gz  FASTQ   DNA   5,817,594  33,528,794,819       17  5,763.3  104,628  1,850  5,530  8,498        0  8,458   60.32   19.57  35.19
