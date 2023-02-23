#!/bin/bash

# For each MinION run, this script will:
# - Run Guppy basecalling in sup mode
# - Concatenate all FASTQ files into 1 file
# - Concatenate the "sequencing_summary" files and run PycoQC for read QC
# - Remove organelle-derived reads

# The working dir is /fs/project/PAS2380/assembly/jelmer

# Inputs and variables
genome_id=GM44
ref_fa_ncbi=data/ref/gmax/ncbi/GCF_000004515.6_Glycine_max_v4.0_genomic.fna
organel_ids="NC_020455.1,NC_007942.1"

# ==============================================================================
#                 PREP SEQS FROM RUN 2022-09-21 (G. max line 44)
# ==============================================================================
# We ran high-accuracy (hac) basecalling on the PC connected to the MinION,
# but this was so slow that <10% of reads were basecalled after about a week.
# We stopped the basecalling and then 'recovered' the reads on that PC.

# Define input files
run_id=GM44_2022-09-21
indir_pass=../minion/TR44_2022-09-21/MCIC_basecalled/fast5_pass/  # These reads were basecalled by the PC connected to the MinION
indir_recov=../minion/TR44_2022-09-21/fast5_recovered/            # These reads were not basecalled by the PC
guppy_config=run/config/dna_r9.4.1_450bps_sup.cfg

# Define output files
guppy_basedir=results/sus44/guppy/"$run_id"
seqsum="$guppy_basedir"/"$run_id"_seqsum_all.txt                    # Concatenated sequencing summary, input for pycoqc
fq_concat_dups=results/sus44/guppy/concat/wdups/"$run_id"_wdups.fastq.gz  # Single FASTQ file with all passed reads, seems to contain some dups due to rescue
fq_concat=results/sus44/guppy/concat/"$run_id".fastq.gz                   # Single FASTQ file with all passed reads, no dups
fq_filt=results/sus44/rm_organel/"$genome_id"/"$run_id".fastq.gz          # Single FASTQ file with all passed reads, no organelle-derived reads

# Guppy (re)basecalling for already basecalled files
for infile in "$indir_pass"/*fast5; do
    outdir="$guppy_basedir"/$(basename "$infile" .fast5)
    sbatch mcic-scripts/ont/guppy_gpu.sh -i "$infile" -o "$outdir" --config "$guppy_config"
done

# Guppy basecalling for 'recovered' FAST5 files
for infile in "$indir_recov"/*.fast5; do
    outdir="$guppy_basedir"/$(basename "$infile" .fast5)
    sbatch mcic-scripts/ont/guppy_gpu.sh -i "$infile" -o "$outdir" \
        --config "$guppy_config" --min_qscore 9
done

# Concatenate all PASS FASTQ files
sbatch mcic-scripts/utils/fqcat.sh -i "$guppy_basedir" -o "$fq_concat_dups" --subdir "pass"

# Remove duplicate reads (with duplicate IDs - incidentally present due to 'rescued' reads ordeal)
sbatch scripts/rmdups.sh "$fq_concat_dups" "$fq_concat"

# Run PycoQC
pyco_out=results/sus44/pycoqc/"$run_id"
sbatch mcic-scripts/ont/seqsum_cat.sh -i "$guppy_basedir" -o "$seqsum"
sbatch mcic-scripts/ont/pycoqc.sh -i "$seqsum" -o "$pyco_out" --min_pass_qual 9
sbatch mcic-scripts/ont/pycoqc.sh -i "$seqsum" -o "$pyco_out"_len1k \
    --min_pass_qual 9 --min_pass_len 1000
#> All/pass/1k+ reads: 4.19/3.40 million/2.82 // All/pass/1k+ bases: 18.45/15.69/15.33 Gbp
#> Pass reads length: median=4,159 / 90%=9,356 // Quality: median=12.49

# Remove organelle-derived reads
sbatch mcic-scripts/ont/rm_organel.sh -i "$fq_concat" -o "$fq_filt" \
    --ref "$ref_fa_ncbi" --seqids "$organel_ids"
#> Nr of input / output sequences in FASTQ files: 3,379,718  /  2,964,870


# ==============================================================================
#             PREP SEQS FROM RUN 2022-12-13 (G. max line 44)
# ==============================================================================
# Define input files
run_id=GM44_2022-12-13
indir=../minion/TR44_2022-12-13/fast5
guppy_config=run/config/dna_r9.4.1_450bps_sup.cfg

# Define output files
guppy_basedir=results/sus44/guppy/"$run_id"
seqsum="$guppy_basedir"/"$run_id"_seqsum_all.txt
fq_concat=results/sus44/guppy/concat/"$run_id".fastq.gz
fq_filt=results/sus44/rm_organel_reads/"$genome_id"/"$run_id".fastq.gz        # Single FASTQ file with all passed reads, no organelle-derived reads

# 1. Guppy basecalling
for infile in "$indir"/*.fast5; do
    outdir="$guppy_basedir"/$(basename "$infile" .fast5)
    sbatch mcic-scripts/ont/guppy_gpu.sh -i "$infile" -o "$outdir" --config "$guppy_config" --min_qscore 9
done

# 2. Concatenate all PASS FASTQ files
sbatch mcic-scripts/utils/fqcat.sh -i "$guppy_basedir" -o test.fastq.gz --subdir "pass"

# 3. Run PycoQC
pyco_out=results/sus44/pycoqc/"$run_id"
sbatch mcic-scripts/ont/seqsum_cat.sh -i "$guppy_basedir" -o "$seqsum"
sbatch mcic-scripts/ont/pycoqc.sh -i "$seqsum" -o "$pyco_out" --min_pass_qual 9
sbatch mcic-scripts/ont/pycoqc.sh -i "$seqsum" -o "$pyco_out"_len1k \
    --min_pass_qual 9 --min_pass_len 1000
#> All/pass/1k+ reads: 3.80/3.26/2.93 million // All/pass/1k+ bases: 24.66/21.73/21.52 Gbp
#> Pass reads length: median=6,680 // 90%=12,434 // Quality: median=13.39

# 4. Remove organelle-derived reads
sbatch mcic-scripts/ont/rm_organel.sh -i "$fq_concat" -o "$fq_filt" \
    --ref "$ref_fa_ncbi" --seqids "$organel_ids"
#> Nr of input / output sequences in FASTQ files: 3,258,336 / 2,852,724


# ==============================================================================
#                   CONCATENATE FINAL FASTQ FILES
# ==============================================================================
fq_concat_all=results/sus44/rm_organel_reads/"$genome_id"/"$genome_id".fastq.gz
fa_concat_all=results/sus44/rm_organel_reads/"$genome_id"/"$genome_id".fasta

# Concatenate
sbatch mcic-scripts/utils/fqcat.sh -i results/sus44/rm_organel_reads/"$genome_id" -o "$fq_concat_all"
# Number of sequences in output file: 5,817,594

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
