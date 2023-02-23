#!/bin/bash

#TODO - Check Canu rerun output

# ==============================================================================
#                        DEFINE VARIABLES FOR THE RUN
# ==============================================================================
# Input files
ref_fa=data/ref/gmax/jgi_Wm82.a4.v1/assembly/Gmax_508_v4.0.fa
ref_gtf=data/ref/gmax/jgi_Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gtf
fq=results/rm_organel/GM44/GM44.fastq.gz
fa_reads=results/rm_organel/GM44/GM44.fasta

# Parameters
genome_id=sus44                            # Susceptible line (resistant line is 145)
basecall_model=r941_min_sup_g507           # We used Guppy 6.4.2, but there are no v6 models
gsize="1g"                                 # Genome size
busco_db=embryophyta_odb10
readlen_smart=1000merged_sc
kraken_db=/fs/ess/PAS0471/jelmer/refdata/kraken/std-plus-fungi
asm_id=merged_scf

# Output files
rata_dir=results/"$genome_id"/ratatosk
pilon_dir=results/"$genome_id"/pilon/"$asm_id"
bbstats_dir="results/$genome_id/bbstats" && mkdir -p "$bbstats_dir"
busco_dir=results/"$genome_id"/busco
kraken_dir=results/"$genome_id"/kraken/"$asm_id"

fq_corr="$rata_dir"/merged/"$genome_id".fastq.gz

asm_flye_raw=results/"$genome_id"/flye/"$genome_id"_flye.fasta
asm_smart_raw=results/"$genome_id"/smartdenovo/smartdenovo.fasta
asm_canu_raw=results/"$genome_id"/canu/canu.fasta
assemblies_raw=( "$asm_smart_raw" "$asm_flye_raw" "$asm_canu_raw" )
asm_smart_medaka=results/"$genome_id"/medaka/smartdenovo/smartdenovo.fasta
asm_canu_medaka=results/"$genome_id"/medaka/canu/canu.fasta
asm_flye_medaka=results/"$genome_id"/medaka/flye/flye.fasta
asm_merged_sc=results/"$genome_id"/quickmerge/q-smart_r-canu/merged_q-smart_r-canu.fasta
asm_merged_scf=results/"$genome_id"/quickmerge/round2_scf/merged_scf.fasta
assemblies_merged=("$asm_merged_sc" "$asm_merged_scf")
asm_pilon="$pilon_dir"/pilon/"$asm_id".fasta
asm_ragtag_c=results/"$genome_id"/ragtag_correct/"$asm_id"/"$asm_id".fasta
asm_ragtag_s=results/"$genome_id"/ragtag_scaffold/"$asm_id"/"$asm_id".fasta
asm_inspector=results/"$genome_id"/inspector_correct/"$asm_id"/"$asm_id".fasta
asm_kraken="$kraken_dir"/unclassified/"$asm_id".fa
gapcloser_dir=results/"$genome_id"/tgs_gapcloser/"$asm_id"
asm_gapcl="$gapcloser_dir"/"$asm_id".fasta

# ==============================================================================
#                   CORRECT LONG READS WITH SHORT READS
# ==============================================================================
# Split FASTQ file
micromamba activate /fs/ess/PAS0471/jelmer/conda/seqkit
seqkit split2 -s 100000 --out-dir results/"$genome_id"/fq_split "$fq"

# Run Ratatosk 
mkdir -p "$rata_dir"
fq_short_list="$rata_dir"/fq_short_list.txt
ls -1 results/rnaseq/star/sclerotinia/unmapped/T*fastq.gz > "$fq_short_list"
for fq_long in results/"$genome_id"/fq_split/*fastq.gz; do
    sbatch -A PAS2380 mcic-scripts/ont/ratatosk.sh -o "$rata_dir" \
        --fq_short_list "$fq_short_list" --fq_long "$fq_long" --insert_size 100
done

# Merge FASTQ files back into one & create a FASTA
sbatch -A PAS2380 mcic-scripts/utils/fqcat.sh -i "$rata_dir" -o "$fq_corr"

# ==============================================================================
#                                  ASSEMBLY
# ==============================================================================
# Run Flye
sbatch -A PAS2380 mcic-scripts/assembly/flye.sh -i "$fq" -o "$asm_flye_raw" \
    --genome-size "$gsize" --iterations 3 

# Run SmartDenovo 
sbatch -A PAS2380 mcic-scripts/assembly/smartdenovo.sh -i "$fq" -o "$asm_smart_raw" \
    --readlen "$readlen_smart"

# Run Canu
canu_workdir=/fs/scratch/PAS0471/"$USER"/tiff/canu/"$genome_id"
sbatch -A PAS2380 mcic-scripts/assembly/canu.sh \
    -i "$fq" --work-dir "$canu_workdir" --out-prefix "$genome_id" \
    --genome-size "$gsize" --time "30:00:00" --more-args "correctedErrorRate=0.16"
bash mcic-scripts/assembly/canu_cp.sh -i "$canu_workdir" -o "$asm_canu_raw" # Need to copy output files afterwards

#! =============================================================================
asm_canu_corr=results/"$genome_id"/canu_corr/"$genome_id"_canu-pr.fasta
canu_workdir=/fs/scratch/PAS0471/"$USER"/tiff/canu/"$genome_id"_corr_err0.105
sbatch -A PAS2380 mcic-scripts/assembly/canu.sh \
    -i "$fq_corr" --work-dir "$canu_workdir" --out-prefix "$genome_id" \
    --genome-size "$gsize" --time "30:00:00" --more-args "correctedErrorRate=0.105"
bash mcic-scripts/assembly/canu_cp.sh -i "$canu_workdir" -o "$asm_canu_corr" # Need to copy output files afterwards

#TODO Try with uncorrected reads and lower error rate?
#canu_workdir=/fs/scratch/PAS0471/"$USER"/tiff/canu/"$genome_id"_corr_err0.16
#sbatch -A PAS2380 mcic-scripts/assembly/canu.sh \
#    -i "$fq_corr" --work-dir "$canu_workdir" --out-prefix "$genome_id" \
#    --genome-size "$gsize" --time "30:00:00" --more-args "correctedErrorRate=0.16"
#bash mcic-scripts/assembly/canu_cp.sh -i "$canu_workdir" -o "$asm_canu_raw" # Need to copy output files afterwards
#! =============================================================================

# ==============================================================================
#                           ASSEMBLY QC 1
# ==============================================================================
# Quick stats with BBtools
micromamba activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96
bbstats_dir="results/$genome_id/bbstats" && mkdir -p "$bbstats_dir"
for asm in "${assemblies_raw[@]}"; do
    asm_id=$(basename "$asm" | sed 's/\..*//') 
    stats.sh "$asm" > "$bbstats_dir"/"$asm_id".txt
    cat "$bbstats_dir"/"$asm_id".txt
    echo -e "-----------------------------------------------\n"
done
#statswrapper.sh in=$merged_scf,$merged_csf format=6 2>/dev/null | column -t tmp.txt | less -S
#>                  size       N50             % in contigs > 50kb      total contigs
#> canu             1062 MB    1292/190 KB     81.32%                   15,095
#> flye             886 MB     1200/200 KB     87.09%                   11,235
#> smartdenovo      957 MB     1105/222 KB     90.99%                   7,506                              

# Busco
for asm in "${assemblies_raw[@]}"; do
    assembler=$(basename "$asm" .fasta)
    outdir=results/"$genome_id"/busco/"$assembler"
    sbatch -A PAS2380 mcic-scripts/assembly/busco.sh -i "$asm" -o "$outdir" -d "$busco_db"
done
#> reference    C:99.8%[S:42.1%,D:57.7%],F:0.1%,M:0.1%,n:1614 
#> canu         C:98.5%[S:45.0%,D:53.5%],F:0.6%,M:0.9%,n:1614
#> flye         C:98.6%[S:46.2%,D:52.4%],F:1.1%,M:0.3%,n:1614
#> smartdenovo  C:99.0%[S:47.2%,D:51.8%],F:0.4%,M:0.6%,n:1614

# Quast
quast_dir=results/"$genome_id"/quast/all_raw
sbatch -A PAS2380 mcic-scripts/assembly/quast.sh -o "$quast_dir" \
    --ref_fa "$ref_fa" --ref_annot "$ref_gtf" --fragmented --large --kmer_stats \
    "${assemblies_raw[@]}"
#>              genome-%        kmer-compl. #missass.   #mismatch   #genes
#> canu         89.07           56.96%      8,903       1358        219k     
#> flye         83.12           61.88       4,403       1030        162k                 
#> smartdenovo  86.57           61.62       4,734       968         179k         

# Inspector
for asm in "${assemblies_raw[@]}"; do
    outdir=results/"$genome_id"/inspector/"$(basename "$asm" .fasta)"
    sbatch -A PAS2380 mcic-scripts/assembly/inspector.sh --assembly "$asm" --reads "$fq" -o "$outdir"
done
#>              QV      split-read  struct.err  sc-ass.err/Mbp  base sub.   maprate
#> canu         29.60   6.9%        9           2058            1.8M        98.97%      
#> flye         40.03   14.31%      1           95              47k         98.89%    
#> smart        31.77   8.37%       12          644             164k        98.96%

# ==============================================================================
#                   LONG-READ POLISHING WITH RACON AND MEDAKA
# ==============================================================================
# Long-read polishing with Racon (2 rounds) for non-Flye assemblies
for asm in "$asm_smart_raw" "$asm_canu_raw"; do
    asm_id=$(basename "$asm" .fasta)
    racon_outdir=results/"$genome_id"/racon/"$asm_id"
    sbatch -A PAS2380 mcic-scripts/assembly/racon.sh --assembly "$asm" --reads "$fq" -o "$racon_outdir"
done

# Long-read polishing with Medaka
for asm in "${assemblies_raw[@]}"; do
    asm_id=$(basename "$asm" .fasta)
    [[ ! "$asm" =~ flye ]] && asm=results/"$genome_id"/racon/"$asm_id"/"$asm_id"_racon2.fasta
    asm_out=results/"$genome_id"/medaka/"$asm_id"/"$asm_id".fasta
    sbatch -A PAS2380 mcic-scripts/ont/medaka.sh -o "$asm_out" \
        --reads "$fq" --assembly "$asm" --model "$basecall_model"
done

# ==============================================================================
#                     MERGE ASSEMBLIES WITH QUICKMERGE
# ==============================================================================
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$asm_smart_medaka" --ref "$asm_canu_medaka" -o "$asm_merged_sc"
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$asm_merged_sc" --ref "$asm_flye_medaka" -o "$asm_merged_scf"

# ==============================================================================
#                      POST-MERGING ASSEMBLY QC
# ==============================================================================
# Quick stats with BBtools
micromamba activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96
for asm in "${assemblies_merged[@]}"; do
    asm_id=$(basename "$asm" | sed 's/\..*//') 
    stats.sh "$asm" > "$bbstats_dir"/"$asm_id".txt
    cat "$bbstats_dir"/"$asm_id".txt
done
#> assembly      size        N50             % in contigs > 50kb      total contigs
#> sc            968 MB      740/315 KB      91.80%                   6,592
#> scf           975 MB      561/435 KB      92.49%                   5,916

# Busco
for asm in "${assemblies_merged[@]}"; do
    outdir="$busco_dir"/$(basename "$asm" .fasta)
    sbatch -A PAS2380 mcic-scripts/assembly/busco.sh -i "$asm" -o "$outdir" -d "$busco_db"
done
#> reference    C:99.8%[S:42.1%,D:57.7%],F:0.1%,M:0.1%,n:1614 
#> sc           C:99.4%[S:44.7%,D:54.7%],F:0.2%,M:0.4%,n:1614
#> scf          C:99.2%[S:43.6%,D:55.6%],F:0.4%,M:0.4%,n:1614 

# Quast
quast_dir=results/"$genome_id"/quast/post_merge3
sbatch -A PAS2380 mcic-scripts/assembly/quast.sh -o "$quast_dir" \
    --ref_fa "$ref_fa" --ref_annot "$ref_gtf" --fragmented --large --kmer_stats \
    "${assemblies_merged[@]}"
#               genome-%        kmer-compl. #missass.   #mismatch   #genes
#> sc           88.1%           64.56%      5,542       1040        195.2k
#> scf          87.38%          62.99%      5,630       1085        195.7k               

# Inspector
for asm in "${assemblies_merged[@]}"; do
    outdir=results/"$genome_id"/inspector/"$(basename "$asm" .fasta)"
    sbatch -A PAS2380 mcic-scripts/assembly/inspector.sh --assembly "$asm" --reads "$fq" -o "$outdir"
done
#>          QV      split-read  struct.err  sc-ass.err/Mbp  base sub.   maprate
# sc        35.95   8.87%       66          479             114k        98.89%
# scf       32.39   9.22%       99          503             147k        98.87%

# ==============================================================================
#                 POST-MERGING LONG-READ POLISHING
# ==============================================================================

#! =============================================================================
#TODO Another round of Racon and/or Medaka? Use corrected reads?
#! Submitted Jan 17
# Long-read polishing with Racon (2 rounds) for non-Flye assemblies
racon_outdir=results/"$genome_id"/racon/$(basename "$asm_merged_scf" .fasta)
sbatch -A PAS2380 mcic-scripts/assembly/racon.sh --assembly "$asm_merged_scf" --reads "$fq" -o "$racon_outdir"

# Long-read polishing with Medaka
asm_id=$(basename "$asm_merged_scf" .fasta)
asm_racon=results/"$genome_id"/racon/"$asm_id"/"$asm_id"_racon2.fasta
asm_medaka=results/"$genome_id"/medaka/"$asm_id"/"$asm_id".fasta
sbatch -A PAS2380 mcic-scripts/ont/medaka.sh -o "$asm_medaka" \
    --reads "$fq" --assembly "$asm_racon" --model "$basecall_model"
#! =============================================================================

# Inspector-correct 
inspector_dir=results/"$genome_id"/inspector/"$asm_id"
sbatch -A PAS2380 mcic-scripts/assembly/inspector.sh --assembly "$asm" --reads "$fq" -o "$inspector_dir"
sbatch -A PAS2380 mcic-scripts/assembly/inspector_correct.sh --inspector_dir "$inspector_dir"  -o "$asm_inspector"

# ==============================================================================
#                   SHORT-READ POLISHING WITH PILON
# ==============================================================================
# Mapping RNAseq reads to the genome before running Pilon
sbatch -A PAS2380 mcic-scripts/rnaseq/star_index.sh -i "$asm_inspector" -o "$pilon_dir"/star_idx
for fq_rna in results/rnaseq/trimgalore/trimmed/T*fastq.gz; do
    sbatch -A PAS2380 mcic-scripts/rnaseq/star_align.sh \
        --single_end --index_bam \
        --reads "$fq_rna" --index_dir "$pilon_dir"/star_idx -o "$pilon_dir"/star
done

# Run Pilon #TODO - redo with inspector-corrected - submitted Jan 18
sbatch -A PAS2380 mcic-scripts/assembly/pilon.sh -o "$asm_pilon" \
    --assembly "$asm_inspector" --bam_dir "$pilon_dir"/star/bam --fix bases

# ==============================================================================
#                      FILTER THE ASSEMBLY
# ==============================================================================
# Run purge_dups - #! Seems to remove too much - change settings or omit
#TODO- read https://academic.oup.com/bioinformatics/article/36/9/2896/5714742
#TODO - Alt software, HapSolo, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7788845/
sbatch -A PAS2380 mcic-scripts/assembly/purge_dups.sh \
    -o results/sus44/purge_dups/"$asm_id" --assembly "$asm_pilon" --reads "$fq"
# 2: scf        975 MB     561/435 KB      92.49%       5,916
# scf_purged    820 MB     423/494 KB      95.97%       3,866

# Run Kraken2 to remove contaminant contigs
sbatch -A PAS2380 mcic-scripts/meta/kraken.sh --unclassified-out --classified-out \
    -i "$asm_pilon" -o "$kraken_dir" --db-dir "$kraken_db"
#> scf: 5 sequences classified, 5,911 unclassified
#>grep -v "^>" results/sus44/kraken/consensus/classified/merged_scf.fa | while read l; do echo "${#l}"; done

#TODO - Remove low-coverage contigs? Blobtools2?

# ==============================================================================
#                     REFERENCE-GUIDED SCAFFOLDING
# ==============================================================================
# Reference-guided misassembly correction with Ragtag #! Do this, or not?
sbatch -A PAS2380 mcic-scripts/assembly/ragtag_correct.sh -o "$asm_ragtag_c" \
    --assembly "$asm_kraken" --reference "$ref_fa" --reads "$fq" 

# Reference-guided scaffolding with Ragtag #TODO - Currently done without Kraken (or purge_dups)
sbatch -A PAS2380 mcic-scripts/assembly/ragtag_scaffold.sh -o "$asm_ragtag_s" \
    --assembly "$asm_ragtag_c" --reference "$ref_fa"

# ==============================================================================
#                     GAP CLOSING
# ==============================================================================
# TGS-GapCloser #TODO redo with Kraken+insp-correct output
sbatch -A PAS2380 mcic-scripts/assembly/tgs_gapcloser.sh \
    --outfile "$asm_gapcl" --assembly "$asm_ragtag_s" --reads "$fa_reads"

# ==============================================================================
#                      FINAL QC
# ==============================================================================
asm="$asm_gapcl"

# Lengths and assignment to chromosomes
ma /fs/ess/PAS0471/jelmer/conda/seqkit
seqkit fx2tab --length --name --header-line "$asm" | sort -k2,2nr # Shortest is just over 10 kb
seqkit fx2tab --length --name --header-line "$asm" | grep "^Gm" | awk '{s+=$2} END{print s}'
seqkit fx2tab --length --name --header-line "$asm" | grep -v "^Gm" | awk '{s+=$2} END{print s}'
#> 929,606,913 assigned to Gm chromosomes
#> 46,309,759 not assigned to chromosomes

# Quick stats with BBtools
micromamba activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96
stats.sh -Xmx4G "$asm" > "$bbstats_dir"/"$asm_id".txt
# STATS WITHOUT KRAKEN #TODO Replace with after-Kraken
# assembly          n_contigs   n_scaffolds scaffold-N50    contig-N50  %-in-scaf>50kb
# scf_pilon         5,916                                   435 KB      92.49%        
# scf_ragtag_c      6,531                                   372 KB      91.59%   
# scf_ragtag_s      6,531       712         46.72 MB        372 KB      98.81%
# scf_gapcloser     2,969       712         47.21 MB        755 KB      98.92%

# Busco
sbatch -A PAS2380 mcic-scripts/assembly/busco.sh -i "$asm" -o "$busco_dir"/"$asm_id" -d "$busco_db"

# Quast
sbatch -A PAS2380 mcic-scripts/assembly/quast.sh -o results/"$genome_id"/quast \
    --nanopore "$fq" --ref_fa "$ref_fa" --ref_annot "$ref_gtf" \
    --fragmented --large --kmer_stats "$asm"

# Inspector
sbatch -A PAS2380 mcic-scripts/assembly/inspector.sh \
    --assembly "$asm" --reads "$fq" -o results/"$genome_id"/inspector/"$asm_id"


#===============================================================================
#                               TRY MUMMER
#===============================================================================
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/mummer4

# Get the Soybase reference
mkdir -p data/ref/soybase
wget https://www.soybase.org/data/v2/Glycine/max/genomes/Wm82.gnm4.4PTR/glyma.Wm82.gnm4.4PTR.genome_main.fna.gz
gunzip data/ref/soybase/glyma.Wm82.gnm4.4PTR.genome_main.fna.gz # Unzip the gzipped file

ref=data/ref/soybase/glyma.Wm82.gnm4.4PTR.genome_main.fna
asm=results/sus44/assembly/sus44_2023-01-19.fasta

# Only get chromosome 16
ref_c16=results/mummer/Wm82_Gm16.fa
awk '/^>/ {P=index($0,"Gm16")>0} {if(P) print}' "$ref" > "$ref_c16"
#grep "Gm16" "$ref"
# Gene: 16g193900

sbatch mcic-scripts/misc/mummer.sh "$ref_c16" "$asm" results/mummer


sbatch -A PAS2380 -t 1:00:00 --wrap="
    rsync -av --progress /fs/project/PAS0471/jelmer/assist/2022-10_tiffana/ /fs/ess/PAS2380/assembly/jelmer/;
    echo 'Done'; date"
