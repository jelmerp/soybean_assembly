#!/bin/bash

# 2022-02-21 -- Copy Guppy config
cp /fs/ess/PAS0471/jelmer/software/guppy-6.4.2/data/dna_r10.4.1_e8.2_260bps_sup.cfg run/config

# Get NCBI G. max genome
mkdir -p data/ref/gmax/ncbi
wget -P data/ref/gmax/ncbi https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/515/GCF_000004515.6_Glycine_max_v4.0/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz
gunzip data/ref/gmax/ncbi/GCF_000004515.6_Glycine_max_v4.0_genomic.fna.gz

# Install Inspector
# - First, I tried regular conda install of https://github.com/ChongLab/Inspector,
#     but got a NumPy attribute error.
# - Then, I installed the older version from https://github.com/ChongLab/Inspector
micromamba create -y -p /fs/ess/PAS0471/jelmer/conda/inspector-1.0.2 -c bioconda python=2.7 minimap2=2.15 samtools=1.9 pysam=0.16.0.1 statsmodels=0.10.1 flye=2.8.3
cd software || exit
git clone https://github.com/Maggi-Chen/Inspector.git
chmod +x Inspector/*py
cp Inspector/*py /fs/ess/PAS0471/jelmer/conda/inspector-1.0.2/bin/