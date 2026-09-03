#!/bin/bash


#SBATCH --job-name=guppy_basecalling
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --gres=gpu:A100:2  #Request 1 A100 GPUs
#SBATCH --time=240
#SBATCH --mem=100G
#SBATCH --partition=gpu
#SBATCH --error=/private/groups/brookslab/smehrete/RNA_modification/basecalling_output/sbatch/w1d1_log/lsf_%j_%x.err
#SBATCH --output=/private/groups/brookslab/smehrete/RNA_modification/basecalling_output/sbatch/w1d1_log/lsf_%j_%x.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smehrete@ucsc.edu

# Define paths
GUPPY_PATH=/private/groups/brookslab/bin/ont-guppy_v6.4.6/bin
INPUT_DIR=/private/groups/brookslab/smehrete/RNA_modification/clone1cscdmsorep1/wt1dmso_3/20231011_1648_MN18795_FAX27637_d643273f/fast5
# Input directory with FAST5 files
OUTPUT_DIR=/private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt1dmso1_output  # Output directory for basecalled data
CONFIG_FILE=/private/groups/brookslab/bin/ont-guppy_v6.4.6/data/rna_r9.4.1_70bps_hac.cfg  # RNA basecalling config for R9.4.1 flow cell


# Create the output directory if it doesn't exist
set -x
mkdir -p $OUTPUT_DIR

# Run Guppy basecaller using the local installation
$GUPPY_PATH/guppy_basecaller -i $INPUT_DIR -s $OUTPUT_DIR \
--config $CONFIG_FILE \
-x "cuda:all" \
--compress_fastq \
--num_callers 2 \
--min_qscore 9  # Set minimum Q-score to 9 (you can change this to 10 if needed)                                                  
