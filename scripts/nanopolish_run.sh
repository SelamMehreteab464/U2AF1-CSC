!/bin/bash
set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)

#move all fast5 files into one directory
./nanopolish index \
	-d /scratch/smehrete/all_sample \
	/scratch/smehrete/all_sample/all_samples.fastq


minimap2 -a -t 12 --secondary no /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.fa /scratch/smehrete/all_sample/all_samples.fastq    | samtools sort -o /scratch/smehrete/bams_m6a_minimap/07172026_combined_m6anet.flair.aligned.bam -T /scratch/smehrete/all_samples.tmp

export HDF5_PLUGIN_PATH=/private/groups/brookslab/smehrete/RNA_modification/nanopolish/usr/local/hdf5/lib/plugin


#script to generate eventaling file was modified to include both read name and read index (nanpolish_eventalign.cpp)
./nanopolish eventalign \
	--reads /scratch/smehrete/all_sample/all_samples.fastq \
	--bam /scratch/smehrete/bams_m6a_minimap/07172026_combined_sorted_m6anet.flair.aligned.bam \
	--genome /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.fa  \
	--scale-events \
	--signal-index  \
	--threads 50 > /scratch/smehrete/07.17.2026.all_samples_cdna_drna.eventalign.txt


