library('DESeq2')

csc_u2af1_illumina_samplefiles<-grep('clone',list.files('Desktop/htseq_count/'),value=TRUE)

csc_u2af1_illumina_samplecondition <- c('U2AF1_S34F_csc', 'U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_S34F_dmso', 'U2AF1_S34F_csc', 'U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_S34F_dmso','U2AF1_WT_csc', 'U2AF1_WT_csc', 'U2AF1_WT_dmso', 'U2AF1_WT_dmso', 'U2AF1_WT_csc', 'U2AF1_WT_csc', 'U2AF1_WT_dmso', 'U2AF1_WT_dmso')

csc_u2af1_illumina_Batch<-c("clone 1", "clone 1", "clone 1", "clone 1","clone 2", "clone 2", "clone 2", "clone 2", "clone 1", "clone 1", "clone 1", "clone 1", "clone 2", "clone 2", "clone 2", "clone 2")


csc_u2af1_illumina_sampleTable<-data.frame(sampleName=csc_u2af1_illumina_samplefiles, fileName=csc_u2af1_illumina_samplefiles, condition=csc_u2af1_illumina_samplecondition, batch=csc_u2af1_illumina_Batch)

csc_u2af1_illumina_ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=csc_u2af1_illumina_sampleTable, directory="Desktop/htseq_count/", design=~batch + condition)


colData(csc_u2af1_illumina_ddsHTSeq)$condition<- factor(colData(csc_u2af1_illumina_ddsHTSeq)$condition, levels=c('U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_WT_csc', 'U2AF1_WT_dmso'))

csc_u2af1_illumina_ddsHTSeq$condition <- relevel(csc_u2af1_illumina_ddsHTSeq$condition, ref = 'U2AF1_WT_dmso')

csc_u2af1_illumina_dds<-DESeq(csc_u2af1_illumina_ddsHTSeq)
