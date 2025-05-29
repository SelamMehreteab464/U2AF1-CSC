#Analysis written and performed by Selam Mehreteab(smehrete@ucsc.edu)


#########################################################################################HT-seq count
#########################################################################################

#M1C2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M1C2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone1_cscrep1.tsv

#M1C3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M1C3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone1_cscrep2.tsv

#M1D2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M1D2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone1_dmsorep1.tsv

#M1D3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M1D3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone1_dmsorep2.tsv

#M2C2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M2C2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone2_cscrep1.tsv

#M2C3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M2C3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone2_cscrep2.tsv

#M2D2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M2D2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone2_dmsorep1.tsv

#M2D3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/M2D3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > s34fclone2_dmsorep2.tsv

#W1C2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W1C2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone1_cscrep1.tsv

#W1C3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W1C3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone1_cscrep2.tsv

#W1D2
ython /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W1D2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone1_dmsorep1.tsv

#W1D3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W1D3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone1_dmsorep2.tsv


#W2C2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W2C2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone2_cscrep1.tsv

#W2C3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W2C3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone2_cscrep2.tsv

#W2D2
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W2D2_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone2_dmsorep1.tsv


#W2D3
python /private/home/smehrete/.local/lib/python3.10/site-packages/HTSeq/scripts/count.py -r pos --stranded=reverse -t exon -i gene_name -m union --nonunique=none /scratch/smehrete/csc_dmso_illumina_alignments/W2D3_Aligned.sortedByCoord.out.bam /private/groups/brookslab/U2AF1-KRAS-RNA-seq/gencode.v33.primary_assembly.annotation.gtf > Wtclone2_dmsorep2.tsv
