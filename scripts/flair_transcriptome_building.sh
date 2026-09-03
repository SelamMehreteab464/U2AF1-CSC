#!/bin/bash
set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)

#building ONT cdna transcriptome
flair align --threads 12 \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	-r /private/groups/brookslab/abehera/Nanopore/HBEC3kt_U2AF1/fastq_guppy4.2.2_gpu/guppy4.2.2_basecalled/all.samples.fastq.gz \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/092525_cdna_FLAIR \
	--quality 0

flair correct --threads 12 \
	-q /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07072026_cdna_FLAIR/flair_align/092525_cdna_FLAIR.bed \
	-f /private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf \
	-j /private/groups/brookslab/smehrete/csc_dmso_illumina_alignments/all_SJ.out.tab \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07072026_cdna_FLAIR/07072026_cdna_FLAIR 

flair collapse --threads 12 \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa --gtf /private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf \
	-q /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07072026_cdna_FLAIR/flair_correct/07072026_cdna_FLAIR_all_corrected.bed \
	-r /private/groups/brookslab/abehera/Nanopore/HBEC3kt_U2AF1/fastq_guppy4.2.2_gpu/guppy4.2.2_basecalled/all.samples.fastq.gz \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07072026_cdna_FLAIR/07072026_cdna_FLAIR \
	--annotation_reliant generate \
	--quality 0 \
	--stringent \
	--check_splice \
	--isoformtss \
	--no_redundant best_only

#making transcriptome for ONT dRNA 
flair align --threads 12 \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	-r /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/all_samples.fastq.gz \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR \
	--nvrn

flair correct --threads 12 \
	-q /private/groups/brookslab/smehrete/RNA_modification/FLAIR/FLAIR.bed \
	-f /private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf \
	-j /private/groups/brookslab/smehrete/csc_dmso_illumina_alignments/all_SJ.out.tab \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/042126 \
	--nvrn 


flair collapse --threads 12 \
	-g /private/groups/brookslab/smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	--gtf /private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf \
	-q /private/groups/brookslab/smehrete/RNA_modification/FLAIR/042126_all_corrected.bed \
	-r /private/groups/brookslab/smehrete/RNA_modification/basecalling_output/all_samples.fastq.gz \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/06042026.FLAIR_corrected_run \
	--generate_map \
	--annotation_reliant generate \
	--stringent \
	--keep_intermediate

#combining the two transcriptomes
flair combine \
	-m /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07072026_ont_drna_FLAIR/flair_collapse/flair_combine_manifest.txt \
	-o /private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome \
	-p 0 \
	-f usageonly



