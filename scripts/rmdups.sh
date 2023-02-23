#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=180
#SBATCH --mem=8G
#SBATCH --output=slurm-rmdup-%j.out

#? This script will remove duplicate reads that were introduced
#? due to mixing "recovered" reads and already
#? basecalled reads -- must have been an overlapping file.

## Command-line args
infile=$1       # FASTQ with duplicates
outfile=$2      # FASTQ without duplicates

## Input and output
dupreads=results/guppy/concat/dupreads.txt

## Report
echo "Starting script rmdups.sh"
date

## Bash script settings
set -euo pipefail

## Load software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/seqkit

## Only select unique reads
zcat "$infile" | \
    seqkit rmdup \
        --dup-num-file "$dupreads" \
        -o "$outfile"

## Report
echo "=========================================================================="
n_in=$(zcat "$infile" | awk '{ s++ } END{ print s/4 }')
n_out=$(zcat "$outfile" | awk '{ s++ } END{ print s/4 }')
echo "Nr of reads in input/output file: $n_in $n_out"
echo
echo "Listing the input and output file:"
ls -lh "$infile" "$outfile"
echo "Done with script"
date
