#!/bin/bash

# Variables
fqdir=results/rnaseq/fastq_concat_all
ref_fasta_pathogen=data/ref/sclerotinia/GCA_001857865.1_ASM185786v1_genomic.fna
ref_gtf_pathogen=data/ref/sclerotinia/GCA_001857865.1_ASM185786v1_genomic.gtf

# ==============================================================================
#                         SET UP AND QC
# ==============================================================================
# Get reference genomes
mkdir -p data/ref/sclerotinia data/ref/gmax
wget -P data/ref/sclerotinia https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/865/GCA_001857865.1_ASM185786v1/GCA_001857865.1_ASM185786v1_genomic.fna.gz
wget -P data/ref/sclerotinia https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/857/865/GCA_001857865.1_ASM185786v1/GCA_001857865.1_ASM185786v1_genomic.gtf.gz
gunzip data/ref/sclerotinia/GCA_001857865.1_ASM185786v1_genomic.fna.gz
gunzip data/ref/sclerotinia/GCA_001857865.1_ASM185786v1_genomic.gtf.gz
cp -r ../2022-08_linda/data/ref/gmax_jgi_Wm82.a4.v1 data/ref/gmax/jgi_Wm82.a4.v1

# Run FastQC & MultiQC on non-concat files
for fq in data/rnaseq/fastq_raw_all/*gz; do
    sbatch -A PAS2380 mcic-scripts/qc/fastqc.sh -i "$fq" -o results/rnaseq/fastqc 
done
sbatch -A PAS2380 mcic-scripts/qc/multiqc.sh -i results/rnaseq/fastqc -o results/rnaseq/multiqc/fastqc

# Concatenate files from different lanes
mkdir -p data/meta && samplesheet=data/meta/samplesheet_forconcat.tsv
paste <(ls -1 data/rnaseq/fastq_raw_all/*gz) \
    <(ls -1 data/rnaseq/fastq_raw_all/*gz | xargs -n 1 basename | sed -E 's/run.*([PST][0-9][^_]*)_.*/\1/') \
    > "$samplesheet"
sbatch -A PAS2380 mcic-scripts/utils/fqconcat_smpsheet.sh -i "$samplesheet" -o "$fqdir"

# Make links to FASTQ files in separate dirs for 44 (T) and 145 (S) (`P` => Sclerotinia culture)
mkdir -p results/rnaseq/fastq_concat_44 results/rnaseq/fastq_concat_145
ln -s -t results/rnaseq/fastq_concat_44 "$PWD"/"$fqdir"/T*fastq.gz
ln -s -t results/rnaseq/fastq_concat_145 "$PWD"/"$fqdir"/S*fastq.gz


# ==============================================================================
#                         REMOVE SCLEROTINIA READS
# ==============================================================================
# Run TrimGalore - I am clipping the first 3 bases since they were consistently of poor qual
for fq in "$fqdir"/[ST]*fastq.gz; do
    sbatch -A PAS2380 mcic-scripts/trim/trimgalore.sh -i "$fq" -o results/rnaseq/trimgalore \
        --length 50 --single_end --more_args "--clip_R1 3"
done

# Map to Sclerotinia with STAR
star_idx_dir=results/rnaseq/star_index/sclerotinia
star_outdir=results/rnaseq/star/sclerotinia

sbatch -A PAS2380 mcic-scripts/rnaseq/star_index.sh \
    -i "$ref_fasta_pathogen" -o "$star_idx_dir" -a "$ref_gtf_pathogen" 

for fq in results/rnaseq/trimgalore/trimmed/*fastq.gz; do
    sbatch -A PAS2380 mcic-scripts/rnaseq/star_align.sh -i "$fq" -r "$star_idx_dir" \
        -o "$star_outdir" -a "$ref_gtf_pathogen" -u -P
done

# Check STAR mapping rates with MultiQC
sbatch -A PAS2380 mcic-scripts/qc/multiqc.sh -i "$star_outdir" -o results/rnaseq/multiqc/star

# Stats
sbatch -A PAS2380 -c 15 --wrap="
    ml python && source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    seqkit stats -j 15 --tabular results/rnaseq/star/sclerotinia/unmapped/T*q.gz > results/rnaseq/star/sclerotinia/unmapped/seqkit_stats.tsv
    cat results/rnaseq/star/sclerotinia/unmapped/seqkit_stats.tsv"
tail -n +2 results/rnaseq/star/sclerotinia/unmapped/seqkit_stats.tsv |
    awk -F"\t" '{sumlen+=$5}END{print sumlen}'

sbatch -A PAS2380 -c 15 --wrap="
    ml python && source activate /fs/ess/PAS0471/jelmer/conda/seqkit
    seqkit stats -j 15 --tabular results/rnaseq/fastq_concat_all/T*q.gz > results/rnaseq/fastq_concat_all/seqkit_stats.tsv
    cat results/rnaseq/fastq_concat_all/seqkit_stats.tsv"
tail -n +2 results/rnaseq/fastq_concat_all/seqkit_stats.tsv |
    awk -F"\t" '{sumlen+=$5; sumreads+=$4}END{print sumlen, sumreads}'

fastqc_smr=results/rnaseq/multiqc/fastqc/multiqc_data/multiqc_fastqc.txt
tail -n +2 "$fastqc_smr" | awk -F"\t" '{nreads+=$5}END{print nreads}'
