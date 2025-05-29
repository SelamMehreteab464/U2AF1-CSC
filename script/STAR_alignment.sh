#Script and analysis was written and performed by Selam Mehreteab(smehrete@ucsc.edu)


#########################################################################################This is the script used to do STAR alignment on our samples. We have 4 samples per condition(2 clones and 2 replicates per clone) and we have 4 conditions(WT DMSO, WT CSC, S34F DMSO, S34F CSC)
#########################################################################################

#Sample Wt clone 1 DMSO rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1D2_S24_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1D2_S24_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate 
	--twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W1D2_; echo "w1d2 alignment done :)" | mail -s "w1d2 alignment done :)" smehrete@ucsc.edu

#Sample Wt clone 1 DMSO rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1D3_S32_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1D3_S32_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W1D3_; echo "W1D3 alignment done :)" | mail -s "W1D3 alignment done :)" smehrete@ucsc.edu


#Sample Wt clone 1 CSC rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1C2_S25_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1C2_S25_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W1C2_; echo "W1C2 alignment done :)" | mail -s "W1C2 alignment done :)" smehrete@ucsc.edu

#Sample Wt clone 1 CSC rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1C3_S33_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W1C3_S33_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W1C3_; echo "W1C3 alignment done :)" | mail -s "W1C3 alignment done :)" smehrete@ucsc.edu

#Sample Wt clone 2 DMSO rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2D2_S22_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2D2_S22_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W2D2_; echo "W2D2 alignment done :)" | mail -s "W2D2 alignment done :)" smehrete@ucsc.edu

#Sample Wt clone 2 DMSO rep 2
private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2D3_S30_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2D3_S30_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W2D3_; echo "W2D3 alignment done :)" | mail -s "W2D3 alignment done :)" smehrete@ucsc.edu


#Sample Wt clone 2 CSC Rep 1
private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2C2_S23_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2C2_S23_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W2C2_; echo "W2C2 alignment done :)" | mail -s "W2C2 alignment done :)" smehrete@ucsc.edu

#Sample Wt clone 2 CSC rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2C3_S31_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/W2C3_S31_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/W2C3_; echo "W2C3 alignment done :)" | mail -s "W2C3 alignment done :)" smehrete@ucsc.edu

#Sample Mt Clone 1 DMSO rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1D2_S28_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1D2_S28_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M1D2_; echo "M1D2 alignment done :)" | mail -s "M1D2 alignment done :)" smehrete@ucsc.edu

#Sample Mt Clone 1 DMSO rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
    --runThreadN 14 
    --genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/
    --readFilesCommand zcat
    --readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1D3_S36_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1D3_S36_L002_R2_001.fastq.gz 
    --outSAMtype BAM SortedByCoordinate
    --twopassMode Basic --quantMode GeneCounts --bamRemoveDuplicatesType UniqueIdentical     --sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf
    --outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M1D3_; echo "M1D3 alignment done :)" | mail -s "M1D3 alignment done :)" smehrete@ucsc.edu


#Sample Mt Clone 1 CSC rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1C2_S29_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1C2_S29_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M1C2_; echo "M1C2 alignment done :)" | mail -s "M1C2 alignment done :)" smehrete@ucsc.edu


#Sample Mt clone 1 CSC rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1C3_S37_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M1C3_S37_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M1C3_; echo "M1C3 alignment done :)" | mail -s "M1C3 alignment done :)" smehrete@ucsc.edu

#Sample Mt clone 2 DMSO rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2D2_S26_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2D2_S26_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M2D2_; echo "M2D2 alignment done :)" | mail -s "M2D2 alignment done :)" smehrete@ucsc.edu

#Sample Mt clone 2 DMSO rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2D3_S34_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2D3_S34_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M2D3_; echo "M2D3 alignment done :)" | mail -s "M2D3 alignment done :)" smehrete@ucsc.edu

#sample Mt clone 2 csc rep 1
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2C2_S27_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2C2_S27_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M2C2_; echo "M2C2 alignment done :)" | mail -s "M2C2 alignment done :)" smehrete@ucsc.edu

#sample mt clone 2 csc rep 2
/private/home/abehera/Applications/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR 
	--runThreadN 14 
	--genomeDir /private/groups/brookslab/reference_indices/STARIndex_v3_GRCh38.u2af1.fix.v1_gencode.v33_Amit/ 
	--readFilesCommand zcat 
	--readFilesIn /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2C3_S35_L002_R1_001.fastq.gz /private/groups/brookslab/data.rep/U2AF1-CSC-DMSO-RNAseq/M2C3_S35_L002_R2_001.fastq.gz 
	--outSAMtype BAM SortedByCoordinate --twopassMode Basic --quantMode GeneCounts 
	--bamRemoveDuplicatesType UniqueIdentical 
	--sjdbGTFfile /private/groups/brookslab/reference_annotations/gencode.v33.primary_assembly.annotation.gtf 
	--outFileNamePrefix /scratch/smehrete/csc_dmso_illumina_alignments/M2C3_; echo "M2C3 alignment done :)" | mail -s "M2C3 alignment done :)" smehrete@ucsc.edu


