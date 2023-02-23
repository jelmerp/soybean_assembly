

## 2023-02-22 -- Find out what the correct Guppy config is

- Find out default basecall method (e.g., to know the 'bps', bases-per-second rate)

```bash
json=../minion/TR145_2023-02-16/report_FAW32537_20230216_1139_779bff98.json
tr , "\n" < "$json" | grep "default basecall"

#> "default basecall model":{"string_value":"dna_r10.4.1_e8.2_400bps_fast.cfg"}
```

## 2023-01-19 -- Send an assembly to Tiffanna

- Copied results/sus44/tgs_gapcloser/merged_scf/merged_scf.fasta to results/sus44/assembly/sus44_2023-01-19.fasta
- This assembly does not include the inspector-correct step

## 2023-01-11 -- If wanting to do Inspector-correct before merging

```bash
# Inspector-correct 
for asm in "${assemblies_medaka[@]}"; do
    asm_id=$(basename "$asm" .fasta)
    outdir=results/"$genome_id"/inspector/"$asm_id"
    sbatch -A PAS2380 mcic-scripts/assembly/inspector.sh --assembly "$asm" --reads "$fq_corr" -o "$outdir"
done
for asm in "${assemblies_medaka[@]}"; do
    asm_id=$(basename "$asm" .fasta)
    asm_out=results/"$genome_id"/inspector_correct/"$asm_id"/"$asm_id".fasta
    inspector_dir=results/"$genome_id"/inspector/"$asm_id"
    ls -lhd "$inspector_dir"
    sbatch -A PAS2380 mcic-scripts/assembly/inspector_correct.sh --inspector_dir "$inspector_dir"  -o "$asm_out"
done
```

## 2023-01-11 -- Assemblies with Ratatosk-corrected reads are a bit worse

```bash
asm_flye_corr=results/"$genome_id"/flye_postrata/"$genome_id"_flye-pr.fasta
asm_smart_corr=results/"$genome_id"/smartdenovo_postrata/"$genome_id"_smartdenovo-pr.fasta
asm_canu_corr=results/"$genome_id"/canu_corr/"$genome_id"_canu-pr.fasta
assemblies_corr=( "$asm_smart_corr" "$asm_flye_corr" "$asm_canu_corr" )
assemblies_corr_noflye=( "$asm_smart_corr" "$asm_canu_corr" )

sbatch -A PAS2380 mcic-scripts/assembly/flye.sh -i "$fq_corr" -o "$asm_flye_corr" --genome-size "$gsize" --iterations 3
sbatch -A PAS2380 mcic-scripts/assembly/smartdenovo.sh -i "$fq_corr" -o "$asm_smart_corr" --readlen "$readlen_smartdenovo"

fa_corr="$outdir_rata"/merged/"$genome_id".fasta
sbatch -A PAS2380 --wrap="
    ml python && source activate /fs/ess/PAS0471/jelmer/conda/seqtk
    seqtk seq -a $fq_corr > $fa_corr
    ls -lh $fa_corr"

#> assembly         size       N50             % in contigs > 50kb      total contigs
#> canu             1062 MB    1292/190 KB     81.32%                   15,095
#> canu-corr
#> flye             886 MB     1200/200 KB     87.09%                   11,235
#> flye-corr        895 MB     1342/177 KB     85.33%                   11,967
#> smartdenovo      957 MB     1105/222 KB     90.99%                   7,506                              
#> smart-corr       956 MB     1117/222 KB     91.27%                   7,463

#> reference    C:99.8%[S:42.1%,D:57.7%],F:0.1%,M:0.1%,n:1614 
#> canu         C:98.5%[S:45.0%,D:53.5%],F:0.6%,M:0.9%,n:1614
#> canu-corr
#> flye         C:98.6%[S:46.2%,D:52.4%],F:1.1%,M:0.3%,n:1614
#> flye-corr    C:99.0%[S:44.7%,D:54.3%],F:0.6%,M:0.4%,n:1614
#> smartdenovo  C:99.0%[S:47.2%,D:51.8%],F:0.4%,M:0.6%,n:1614
#> smart-corr   C:99.5%[S:46.6%,D:52.9%],F:0.1%,M:0.4%,n:1614

#               genome-%        kmer-compl. #missass.   #mismatch   #genes
#> canu         #TODO
#> canu-corr
#> flye         83.119          61.88       4,403       8.16M       162k                 
#> flye-corr    83.761          58.94       5,700       9.00M       170k
#> smartdenovo  86.57           61.62       4,734       7.99M       179k         
#> smart-corr   86.48           58.72       5,426       9.30M       184k

# Inspector
#               QV      split-read  struct.err  sc-ass.err/Mbp  Base sub.   Maprate
#> canu         #TODO
#> canu-corr
#> flye         40.03   14.31%      1           95              47k         98.89%    
#> flye-corr    38.37   15.39%      2           137             72k         99.0%
#> smart        31.77   8.37%       12          644             164k        98.96%
#> smart-corr   31.42   8.83%       38          693             193k        98.87%

# STATS WITHOUT KRAKEN
# assembly          n_contigs   n_scaffolds scaffold-N50    contig-N50  %-in-scaf>50kb
# scf_pilon         5,916                                   435 KB      92.49%        
# scf_ragtag_c      6,531                                   372 KB      91.59%   
# scf_ragtag_s      6,531       712         975 MB          372 KB      98.81%
# csf_pilon         11,894                                  372 KB      84.39%
# csf_ragtag_c      12,407                                  311 KB      83.97%
# csf_ragtag_s      12,407      3,616       49 MB           311 KB      95.11% 
```


## 2023-01-11 -- 'scf' assembly is better than 'csf'

- scf: Merge of smartdenovo-as-query with (canu-query + flye-ref)-as-ref
- csf: Merge of canu-as-query with (smartdenovo-query + flye-ref)-as-ref

```bash
merged_cs=results/"$genome_id"/quickmerge/q-canu_r-smart/merged_q-canu_r-smart.fasta
merged_csf=results/"$genome_id"/quickmerge/round2_csf/merged_csf.fasta

sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$asm_canu_medaka" --ref "$asm_smart_medaka"  -o "$merged_cs"
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$merged_cs" --ref "$asm_flye_medaka" -o "$merged_csf"
```

```bash
#> stats.sh
# assembly      size        N50             % in contigs > 50kb      total contigs
# sc            968 MB      740/315 KB      91.80%                   6,592
# cs            1060 MB     891/272 KB      83.62%                   12,706
# 2: scf        975 MB      561/435 KB      92.49%                   5,916
# 2: csf        1059 MB     668/372 KB      84.39%                   11,894
# 3: c-scf      1079 MB     658/385 KB      83.97%                   12,338
# scf_inspector 974         562/432         92.49%                   5,907

#> BUSCO
# reference             C:99.8%[S:42.1%,D:57.7%],F:0.1%,M:0.1%,n:1614 
# sc                    C:99.4%[S:44.7%,D:54.7%],F:0.2%,M:0.4%,n:1614
# cs                    C:99.3%[S:42.8%,D:56.5%],F:0.2%,M:0.5%,n:1614 
# 2: scf                C:99.2%[S:43.6%,D:55.6%],F:0.4%,M:0.4%,n:1614 
# 2: csf                C:99.0%[S:43.1%,D:55.9%],F:0.4%,M:0.6%,n:1614
# scf_dedup             C:98.1%[S:49.8%,D:48.3%],F:0.6%,M:1.3%,n:1614  
# csf_kraken            C:99.0%[S:43.1%,D:55.9%],F:0.4%,M:0.6%,n:1614

#> Inspector
#           QV      split-read  struct.err
# csf       32.05   7.94%       247
# scf       32.39   9.22%       99
```

## 2023-01-07 -- Quality scores usage

- Flye does not use quality scores: https://github.com/fenderglass/Flye/issues/357
- Looks like Canu doesn't either: https://github.com/marbl/canu/issues/1009

## 2023-01-07 -- Duplicate contigs in merged assemblies

In the end, I'm not using the `smartdenovo-flye_canu-fly`-merged assembly,
so this is not necessary:

```bash
#! =============================================================================
# Duplicate entries in sf-cf
# For now, remove duplicates, but note that in 2 of 3 cases, the sequences are not the same!!
# Though note that Pilon will relabel contigs anyway...
grep "^>" "$merged_sfcf" | sort | uniq -c | sort -nr | head
grep -A1 ">contig_7216" "$merged_sfcf" | less -S
grep -A1 ">contig_3826" "$merged_sfcf" | less -S
grep -A1 ">contig_10260" "$merged_sfcf" | less -S

merged_sfcf_nodups=results/sus44/quickmerge/round2_sf-cf/merged_sf-cf_nodups.fasta
micromamba activate /fs/ess/PAS0471/jelmer/conda/seqkit
seqkit rmdup -n "$merged_sfcf" > "$merged_sfcf_nodups"

awk 'BEGIN{RS=">"}{if(NR>1)print ">"$1"_"(NR-1)"\n"$2}' "$merged_sfcf" > test2.fasta
#! =============================================================================
```

## 2023-01-07 -- Scaffolding

- Long-read scaffolding seems to make things worse? Nr of contigs increases dramatically.
  When doing reference-guided scaffolding for the Pilon assembly versus Pilon+Longstitch,
  the Pilon one also looks a lot better!

```bash
# assembly          n_contigs   n_scaffolds scaffold-N50    contig-N50  %-in-scaf>50kb
# scf_longstitch    11,251      9,001       451 KB          284 KB      93.77%
# scf_ragtag_c      11,559      9,463       404 KB          277 KB      93.19%            
# scf_ragtag_s      11,559      4,168       47.2 MB         277 KB      98.11%
```

```bash
# ==============================================================================
#                     LONG-READ SCAFFOLDING
# ==============================================================================
# Longstitch
gsize_int=1e9                              # Genome size as an integer
sbatch -A PAS2380 mcic-scripts/assembly/longstitch.sh --genome_size $gsize_int \
    --assembly "$asm_pilon" --fastq "$fq_corr" -o "$asm_longstitch"

# Links #TODO Should this (still) be run after longstitch?
asm_links=X
sbatch -A PAS2380 mcic-scripts/assembly/links.sh \
    --assembly "$asm_longstitch" --fastq "$fq" -o "$links_dir" #TODO Use corr. reads!
```

## 2023-01-07 -- Quickmerge

Old quickmerge runs:

```bash
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$asm_smart_medaka" --ref "$asm_flye_medaka" -o "$merged_sf" 
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$asm_canu_medaka" --ref "$asm_flye_medaka"  -o "$merged_cf" 
sbatch -A PAS2380 mcic-scripts/assembly/quickmerge.sh --minlen_anchor 200000 \
    --query "$merged_sf" --ref "$merged_cf"  -o "$merged_sfcf"
```

The strategy of merging each of Smartdenovo and Canu with Flye first did lead
to the longest N50, but the Busco completeness score and Inspector QV value was
actually lower.

## 2023-01-05 -- Pilon

Thought I needed to split the assembly for Pilon into many smaller files
due to memory usage, but this was just due to not specifying memory using `java xmx...`.

Code with split assembly:

```bash
# Split assembly into parts -- Pilon uses too much memory to run it all at once
micromamba activate /fs/project/PAS0471/jelmer/conda/seqkit
seqkit split2 --by-length 100M "$assembly" --out-dir "$pilon_dir"/split_assembly

for asm_part in  "$pilon_dir"/split_assembly/*fasta; do
    sbatch -A PAS2380 mcic-scripts/assembly/pilon.sh \
        --genome "$asm_part" --bam_dir "$pilon_dir"/star/bam \
        -o "$pilon_dir"/pilon --out_prefix "$(basename "$asm_part" .fasta)" \
        --fix bases
done
```

## 2022-12-31 -- QuickMerge

```bash
#! 3rd round is not improving anything, remove
merged_round3a=results/"$genome_id"/quickmerge/round3/c-scf/merged_c-scf.fasta
sbatch mcic-scripts/assembly/quickmerge.sh --minlen_anchor 300000 \
    --query "$asm_canu" --ref "$merged_round2" -o "$merged_round3a"

merged_round3b=results/"$genome_id"/quickmerge/round3/scf-c/merged_scf-c.fasta
sbatch mcic-scripts/assembly/quickmerge.sh --minlen_anchor 300000 \
    --query "$merged_round2" --ref "$asm_canu" -o "$merged_round3b"
```

## 2022-12-24 -- QuickMerge

- The better assembly should be the 'query' in Quickmerge,
  but I also tried it the other way around.
  That indeed resulted in a worse assembly than the 'right way around' merger'.

```bash
genome_id=sus44
asm_smart=results/sus44/medaka/sus44_smartdenovo/consensus.fasta
asm_canu=results/sus44/medaka/sus44_canu/consensus.fasta
merged_round1_alt=results/"$genome_id"/quickmerge/q-canu_r-smart/merged_q-canu_r-smart.fasta
sbatch mcic-scripts/assembly/quickmerge.sh --query "$asm_canu" --ref "$asm_smart" \
    -o "$merged_round1_alt" --minlen_anchor 200000

## BBstats output:
# assembly          size           N50                 % in contigs > 50kb      total contigs
# q-smart_r-canu    968.325 MB     740/315.836 KB      91.80%                   6592
# q-canu_r-smart    1060.448 MB    891/272.631 KB      83.62%                   12706
```

## 2022-12-15 -- Racon & Medaka Polishing

- Polishing doesn't have a big effect on basic assembly stats:

```bash
# ==============================================================================
#                      POST-POLISHING ASSEMBLY QC
# ==============================================================================
mapfile -t assemblies < <(find results/"$genome_id"/medaka -name "consensus.fasta")

# Quick stats with BBtools
micromamba activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96
bbstats_dir="results/$genome_id/bbstats/post_medaka" && mkdir -p "$bbstats_dir"
for asm in "${assemblies[@]}"; do
    ls -lh "$asm"
    asm_id=$(basename "$asm" | sed 's/\..*//') 
    stats.sh "$asm" > "$bbstats_dir"/"$asm_id".txt
    cat "$bbstats_dir"/"$asm_id".txt
    echo -e "-----------------------------------------------\n"
done
#               size            N50
# canu          1055.028 MB     1263/195.422 KB
# flye          886.877 MB      1200/200.703 KB     
# smartdenovo   960.590 MB      1106/222.795 KB
```

## 2022-11-24 -- Guppy

- In the past week, updated Guppy to v 6.4.2 and figured out how to use GPUs for basecalling,
  and because the SUP accuracy mode barely took longer than the HAC mode with GPUs,
  I switched to that.

- **Accidentally still used a min. qual score of 9, whereas the default for the updated**
  **Guppy version + SUP config is 10**  

## 2022-11-09 -- Remove organelle reads

- https://www.ncbi.nlm.nih.gov/assembly/GCF_000004515.6#/def_asm_non-nuclear

## 2022-11-08 -- Tested Nanolyse

Tried NanoLyse on the FASTQ subset -- 0 matches, should be able to skip:

```bash
micromamba activate /fs/ess/PAS0471/jelmer/conda/nanolyse-1.2.1
zcat "$fqsub" | NanoLyse | gzip > test.fastq.gz
#> NanoLyse: removed 0 reads
```

## 2022-11-08 -- Use Pilon instead of FMLRC2

Considered polishing with FMLRC2, but not sure if suitable for RNAseq data

```bash
## Short-read polishing with FMLRC2 #TODO - Is this OK?
fq_short=$(ls results/rnaseq/star/sclerotinia/unmapped/S1*fastq.gz)
sbatch mcic-scripts/nanopore/fmlrc2.sh --infile "$asm" --fq_short "$fq_short" \
    -o results/"$genome_id"/fmlrc2/"$asm_id" --eukaryote
#> results/sus44_run1/fmlrc2/smartdenovo/smartdenovo.fasta

grep -v "^>" results/sus44_run1/fmlrc2/smartdenovo/smartdenovo.fasta | wc -c   # 880,524,361
grep -v "^>" results/sus44_run1/all_assemblies/smartdenovo.fasta | wc -c       # 889,633,694
```

## 2022-11-08 -- No Ratatosk

```bash
## Run Ratatosk     # Output: 3,316,426 reads ##TODO - skip this
fq_short_list="$outdir"/fq_short_list.txt
outdir_rata=results/"$genome_id"/ratatosk && mkdir -p "$outdir"
fq_rata="$outdir_rata"/$(basename "$fq")

ls -1 results/rnaseq/star/sclerotinia/unmapped/S1*fastq.gz > "$fq_short_list"
sbatch mcic-scripts/nanopore/ratatosk.sh --fq_short_list "$fq_short_list" \
    --fq_long "$fq" -o "$outdir_rata" --insert_size 100
```

## 2022-11-07 -- No NextDenovo

Skipping `NextDenovo` for now, since it behaves poorly when running,
and seems more sensitive to low coverage than other assemblers.

```bash
## Run NextDenovo => Skip due to lower performance with lower coverage? #TODO retry
outdir=results/"$genome_id"/nextdenovo && mkdir -p "$outdir"
fofn="$outdir"/fq.fofn 
echo "$(realpath $fq)" > "$fofn" # Needs to contain absolute paths!
sbatch mcic-scripts/assembly/nextdenovo.sh -i "$fofn" -o "$outdir" -s "$gsize"
#TODO - check /fs/scratch/PAS0471/jelmer/nextdenovo/sus44_run2022-09-21/03.ctg_graph/03.ctg_cns.sh.work/ctg_cns08/nextdenovo.sh.o
```

## 2022-11-07 -- Misc

- Don't polish `Flye` assemblies with `Racon` before `Medaka`, see https://github.com/nanoporetech/medaka#origin-of-the-draft-sequence
- Perhaps not necessary to run `Racon` at all, see https://github.com/rrwick/August-2019-consensus-accuracy-update

## 2022-11-07 -- Assemblies

Can run Busco with either `fabales_odb10.2019-11-20` (5,366 genes, Order)
or `embryophyta_odb10` (1,614 genes, Subkingdom).
I think `embryophyta_odb10` is more widely used.

`Busco` (`embryophyta_odb10`) results for raw assemblies:

- `Smartdenovo`:   |C:97.0%[S:61.7%,D:35.3%],F:1.5%,M:1.5%,n:1614
- `Canu`:          |C:95.5%[S:53.2%,D:42.3%],F:2.5%,M:2.0%,n:1614
- `Flye`:          |C:96.4%[S:52.5%,D:43.9%],F:2.9%,M:0.7%,n:1614
- `Redbean`:       |C:79.8%[S:73.4%,D:6.4%],F:6.7%,M:13.5%,n:1614

## 2022-11-07 -- Read numbers

Summary files: 4,187,769 total reads / 3,678,949 passed reads

```bash
fq_concat_dups=results/guppy/concat/"$run_id"_wdups.fastq.gz     # Single FASTQ file with all passed reads, seems to contain some dups due to rescue
fq_concat=results/guppy/concat/"$run_id".fastq.gz                # Single FASTQ file with all passed reads, no dups

zcat $fq_concat_dups | wc -l
#13347816 => 3,336,954

zcat $fq_concat | wc -l
#13265704 => 3,316,426
```

## 2022-11-05 -- Raw long-read error correction

First tried `FMLRC2`, which came best out of benchmark papers,
but this only outputs FASTA files.
Especially in out case, where only a fraction of the reads will be corrected
because the short-read data is from RNAseq, this makes no sense --
we need to retain the quality scores.

Tried `Ratatosk` instead.

```bash
fq_short=$(ls results/rnaseq/star/sclerotinia/unmapped/S1*fastq.gz)
sbatch mcic-scripts/nanopore/fmlrc2.sh --fq_short "$fq_short" --fq_long "$fq" -o results/"$genome_id"/fmlrc2
```

## 2022-11-05 -- Mapping Illumina short reads to the genome

See:
- https://denbi-nanopore-training-course.readthedocs.io/en/latest/artic/pilon/Mapping_1.html
- https://timkahlke.github.io/LongRead_tutorials/ECR_P.html

```bash
source activate /fs/ess/PAS0471/jelmer/conda/bwa-mem-2.2.1
bwa index "$genome"
bwa mem -t "$n_threads" "$genome" "$fastq" | \
    samtools view -b - | \
    samtools sort - > $bam
```

## 2022-11-05 -- Flye detected duplicate read IDs

These duplicate read IDs were there due to mixing "recovered" reads and already
basecalled reads -- must have been an overlapping file.

```bash
#![2022-11-05 10:47:17] ERROR: The input contain reads with duplicated IDs. Make sure all reads have unique IDs and restart. The first problematic ID was: b8dc0aa0-cb79-4294-aff2-8821d0269bfc
zgrep "b8dc0aa0-cb79-4294-aff2-8821d0269bfc" "$fq"
#@b8dc0aa0-cb79-4294-aff2-8821d0269bfc runid=61d5b40425e502fa287d79365a779b6fcc6e2609 sampleid=44 read=45557 ch=193 start_time=2022-09-23T16:11:47Z
#@b8dc0aa0-cb79-4294-aff2-8821d0269bfc runid=61d5b40425e502fa287d79365a779b6fcc6e2609 sampleid=44 read=45557 ch=193 start_time=2022-09-23T16:11:47Z

grep "b8dc0aa0-cb79-4294-aff2-8821d0269bfc" results/guppy/run2022-09-21/sequencing_summary_all.txt
zgrep -r "b8dc0aa0-cb79-4294-aff2-8821d0269bfc" results/guppy/run2022-09-21

seqsum=results/guppy/run2022-09-21/sequencing_summary_all.txt
tail -n +2 "$seqsum" | cut -f 2 | sort | uniq -c | sort -nr | awk '$1 > 1' > dupreads_counts.txt
tail -n +2 dupreads_counts.txt | awk '{print $NF}' > dupreads.txt

grep -f dupreads.txt "$seqsum" > seqsum_dupreads.txt # This takes forever
```

## 2022-11-04 -- Run Kraken after mapping to Sclerotinia

This resulted in the detection of _very_ few fungal reads
(and some bacteria, viruses, human) -- don't think it's worth specifically removing these.

```sh
## Run Kraken2 on unmapped reads
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-plant-fungi-protozoa
kraken_outdir=results/rnaseq/kraken/post-sclero-map/
for fq in "$star_outdir"/unmapped/S1*fastq.gz; do
    sbatch --mem=130G -t30 mcic-scripts/meta/kraken.sh -i "$fq" -o "$kraken_outdir" -d "$kraken_db" -s
done
```

## 2022-11-02 -- Percentage Sclerotinia reads in samples

Below are the results from nf-core RNAseq (using BBSplit to map to both Sclerotinia and Soybean).
(See `sandbox/nfcore_rnaseq.sh` for running this workflow -- not included in main.)
Mapping of reads to Sclerotinia genome with STAR results in similar (but slightly lower)
mapping rates.

```sh
cat results/rnaseq/nfc_rnaseq/bbsplit/summary.txt
#> P1      98.51596
#> P2      98.35859
#> P3      98.76110
#> S1-III  0.10855
#> S1-II   0.06309
#> S1-I    0.05304
#> S2-III  3.07190
#> S2-II   3.12651
#> S2-I    4.76603
#> S3-III  2.48674
#> S3-II   6.06794
#> S3-I    3.17456
#> S4-III  2.75851
#> S4-II   6.87724
#> S4-I    12.03746
#> S5-III  4.52513
#> S5-II   10.93413
#> S5-I    7.73383
#> T1-III  0.07412
#> T1-II   0.05406
#> T1-I    0.06117
#> T2-III  5.23400
#> T2-II   4.54485
#> T2-I    3.96657
#> T3-III  3.88660
#> T3-II   4.04851
#> T3-I    4.61344
#> T4-III  6.37922
#> T4-II   28.18617
#> T4-I    13.18633
#> T5-III  27.23632
#> T5-II   23.41938
#> T5-I    35.73355
```

## 2022-11-02 - Kraken on raw RNAseq data

```bash
## Initial check of fungal contents with Kraken
for R1 in "$fqdir"/P*gz "$fqdir"/S1*gz; do
    kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
    outdir=results/rnaseq/kraken/"$(basename "$kraken_db")"
    sbatch mcic-scripts/meta/kraken.sh -i "$R1" -o "$outdir" -d "$kraken_db"
    
    kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-plant-fungi-protozoa
    outdir=results/rnaseq/kraken/"$(basename "$kraken_db")"
    sbatch --mem=130G mcic-scripts/meta/kraken.sh -i "$R1" -o "$outdir" -d "$kraken_db"
done
```

- Ran only P samples (fungal culture) and S1 (plant timepoint 1, no fungal seqs expected)
- Some 'Sclerotinia sclerotiorum mitovirus 2' is detected
- Comparing Kraken the std-plus-fungi and std-plus-plant-fungi-protozoa database,
  less fungal sequences are detected in the latter, but the difference is smaller at
  the `Sclerotiniaceae` level, so some of the others may be false positives?

grep "Fungi" results/rnaseq/kraken/std-plus-fungi/*report.txt 
> results/rnaseq/kraken/P1_report.txt:  8.68      1151494 476     36563910        149009  K       4751            Fungi
> results/rnaseq/kraken/P2_report.txt: 13.96      1792560 496     49217948        146803  K       4751            Fungi
> results/rnaseq/kraken/P3_report.txt:  4.68      538651  521     21240270        143138  K       4751            Fungi
> results/rnaseq/kraken/S1-III_report.txt:  0.01  2748    19      911423  27860   K       4751            Fungi
> results/rnaseq/kraken/S1-II_report.txt:  0.00   519     21      698889  18961   K       4751            Fungi
> results/rnaseq/kraken/S1-I_report.txt:  0.00    765     30      1233059 24680   K       4751            Fungi

grep "Sclerotiniaceae" results/rnaseq/kraken/std-plus-fungi/*report.txt 
> results/rnaseq/kraken/std-plus-fungi/P1_report.txt:  4.83       640676  0       26601875        120172  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-fungi/P2_report.txt:  7.12       913950  0       32442241        118428  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-fungi/P3_report.txt:  2.93       337518  0       17116150        118985  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-fungi/S1-III_report.txt:  0.00   1331    0       59854   15785   F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-fungi/S1-II_report.txt:  0.00    286     0       19881   9312    F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-fungi/S1-I_report.txt:  0.00     387     0       30363   12288   F       28983                             Sclerotiniaceae

grep "Fungi" results/rnaseq/kraken/std-plus-plant-fungi-protozoa/*report.txt
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P1_report.txt:  5.44        721775  13      28614103        143501  K       4751            Fungi
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P2_report.txt:  6.63        851098  34      31693291        141738  K       4751            Fungi
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P3_report.txt:  3.94        453538  2       19794523        138811  K       4751            Fungi
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-III_report.txt:  0.01    1412    0       106122  21432   K       4751            Fungi
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-II_report.txt:  0.00     374     0       54333   13162   K       4751            Fungi
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-I_report.txt:  0.00      537     0       90065   17748   K       4751            Fungi

grep "Sclerotiniaceae" results/rnaseq/kraken/std-plus-plant-fungi-protozoa/*report.txt 
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P1_report.txt:  4.06        538921  0       24496081        119832  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P2_report.txt:  5.33        684081  0       27732372        118279  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/P3_report.txt:  2.77        319429  0       16744728        118721  F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-III_report.txt:  0.00    1083    0       50479   15444   F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-II_report.txt:  0.00     259     0       16312   9025    F       28983                             Sclerotiniaceae
> results/rnaseq/kraken/std-plus-plant-fungi-protozoa/S1-I_report.txt:  0.00      357     0       24188   11951   F       28983                             Sclerotiniaceae

## 2022-10-26

Trying `Canu` with option `correctedErrorRate=0.16`.
See <https://canu.readthedocs.io/en/latest/faq.html#what-parameters-can-i-tweak>

> For less than 30X raw PacBio or Nanopore coverage, increase the allowed difference in
> overlaps by a few percent (from 4.5% to 8.5% (or more) with correctedErrorRate=0.105 for PacBio
> and from 14.4% to 16% (or more) with correctedErrorRate=0.16 for Nanopore), to adjust for inferior read correction.

## 2022-10-15 -- Initial runs with sus44 run2022-09-21

- `SmartDenovo` finished in 15:30 hours
- `Flye` finished in 04:32 hours
- `Redbean/Wtdbg2` finished in 1:07 hours
