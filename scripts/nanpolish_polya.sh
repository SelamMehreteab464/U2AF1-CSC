#!/bin/bash

set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)
head -n 1 /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt1csc1_output/sequencing_summary.txt > all_samples_sequencing_summary.txt
for f in /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt1csc1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt1csc2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt2csc1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt2csc2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt1dmso1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt1dmso2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt2dmso1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/mt2dmso2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt1csc1_output/basecalled_samples/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt1csc2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt2csc1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt2csc2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt1dmso1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt1dmso2_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt2dmso1_output/sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/wt2dmso2_output/sequencing_summary.txt; do     tail -n +2 "$f" >> all_samples_sequencing_summary.txt 
done


./nanopolish index --directory=/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample --sequencing-summary=/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample/all_samples_sequencing_summary.txt /private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample/all_samples.fastq

minimap2 -a -x map-ont -t 12 --secondary=no /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.fa /private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample/all_samples.fastq | samtools view -b - -o /private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/polyAlength.bam
cd /private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/ && samtools sort -T tmp -o polyAlength.sorted.bam polyAlength.bam && samtools index polyAlength.sorted.bam && cd ..

export HDF5_PLUGIN_PATH=/private/groups/brookslab/smehrete/RNA_modification/nanopolish/usr/local/hdf5/lib/plugin

./nanopolish polya --threads=8 --reads=/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/all_sample/all_samples.fastq --bam=/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/polyAlength.sorted.bam --genome=/private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.fa > /private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/polya_results.tsv

