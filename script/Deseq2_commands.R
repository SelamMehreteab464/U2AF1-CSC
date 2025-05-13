library('DESeq2')

csc_u2af1_illumina_samplefiles<-grep('clone',list.files('Desktop/htseq_count/'),value=TRUE)

csc_u2af1_illumina_samplecondition <- c('U2AF1_S34F_csc', 'U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_S34F_dmso', 'U2AF1_S34F_csc', 'U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_S34F_dmso','U2AF1_WT_csc', 'U2AF1_WT_csc', 'U2AF1_WT_dmso', 'U2AF1_WT_dmso', 'U2AF1_WT_csc', 'U2AF1_WT_csc', 'U2AF1_WT_dmso', 'U2AF1_WT_dmso')

csc_u2af1_illumina_Batch<-c("clone 1", "clone 1", "clone 1", "clone 1","clone 2", "clone 2", "clone 2", "clone 2", "clone 1", "clone 1", "clone 1", "clone 1", "clone 2", "clone 2", "clone 2", "clone 2")


csc_u2af1_illumina_sampleTable<-data.frame(sampleName=csc_u2af1_illumina_samplefiles, fileName=csc_u2af1_illumina_samplefiles, condition=csc_u2af1_illumina_samplecondition, batch=csc_u2af1_illumina_Batch)

csc_u2af1_illumina_ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=csc_u2af1_illumina_sampleTable, directory="Desktop/htseq_count/", design=~batch + condition)


colData(csc_u2af1_illumina_ddsHTSeq)$condition<- factor(colData(csc_u2af1_illumina_ddsHTSeq)$condition, levels=c('U2AF1_S34F_csc', 'U2AF1_S34F_dmso', 'U2AF1_WT_csc', 'U2AF1_WT_dmso'))

csc_u2af1_illumina_ddsHTSeq$condition <- relevel(csc_u2af1_illumina_ddsHTSeq$condition, ref = 'U2AF1_WT_dmso')

csc_u2af1_illumina_dds<-DESeq(csc_u2af1_illumina_ddsHTSeq)

#script to plot gene expresssion pca plot
csc_u2af1_illumina_vst <- vst(csc_u2af1_illumina_dds, blind=FALSE)
csc_u2af1_illumina_mat <- assay(csc_u2af1_illumina_vst)
assay(csc_u2af1_illumina_vst) <- csc_u2af1_illumina_mat
csc_u2af1_illumina_mat<- assay(csc_u2af1_illumina_vst)
plotPCA(csc_u2af1_illumina_vst, intgroup=c("condition", "batch"))
plotPCA(csc_u2af1_illumina_vst, intgroup=c("condition", "batch"), returnData=TRUE)

plotCounts(csc_u2af1_illumina_dds, gene='U2AF1', intgroup=c('condition', 'batch'), returnData=TRUE)

#script to plot gene expression volcano plot
LRTbatch_csc_u2af1_illumina_dds <- DESeq(csc_u2af1_illumina_dds, test="LRT", reduced=~batch)

u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso<- results(LRTbatch_csc_u2af1_illumina_dds, contrast=c("condition", "U2AF1_S34F_dmso", "U2AF1_WT_dmso"))

u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso <- results(LRTbatch_csc_u2af1_illumina_dds, contrast=c("condition", "U2AF1_WT_csc", "U2AF1_WT_dmso"))

u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso <- results(LRTbatch_csc_u2af1_illumina_dds, contrast=c("condition", "U2AF1_S34F_csc", "U2AF1_WT_dmso"))

library(EnhancedVolcano)

EnhancedVolcano(u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso, lab=rownames(u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso), x='log2FoldChange', y='pvalue', title="u2af1 wt dmso vs. u2af1 s34f csc", xlim = c(-15, 15), ylim = c(0,300))

EnhancedVolcano(u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso, lab=rownames(u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso), x='log2FoldChange', y='pvalue', title="u2af1 illumina wt dmso vs. u2af1 wt csc", xlim = c(-15, 15), ylim = c(0,300))

EnhancedVolcano(u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso, lab=rownames(u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso), x='log2FoldChange', y='pvalue', title="u2af1 wt dmso vs. u2af1 s34f csc", xlim = c(-15, 15), ylim = c(0,300))

#filtering for pval 
>u2af1_illumina_filtered_s34f_dmso_vs_u2af1_illumina_wt_dmso <- u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso[!is.na(u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso$padj) & u2af1_illumina_s34f_dmso_vs_u2af1_illumina_wt_dmso$padj<0.05,]


u2af1_illumina_filtered_u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso <- u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso[!is.na(u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso$padj) & u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso$padj<0.05,]


u2af1_illumina_filtered_s34f_csc_vs_u2af1_illumina_wt_dmso <- u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso[!is.na(u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso$padj) & u2af1_illumina_s34f_csc_vs_u2af1_illumina_wt_dmso$padj<0.05,]


write.csv(as.data.frame(u2af1_illumina_filtered_s34f_dmso_vs_u2af1_illumina_wt_dmso), file='~/Desktop/u2af1_illumina_filtered_s34f_dmso_vs_u2af1_illumina_wt_dmso')

write.csv(as.data.frame(u2af1_illumina_filtered_u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso), file='~/Desktop/u2af1_illumina_filtered_u2af1_illumina_wt_csc_vs_u2af1_illumina_wt_dmso')

write.csv(as.data.frame(u2af1_illumina_filtered_s34f_csc_vs_u2af1_illumina_wt_dmso), file='~/Desktop/u2af1_illumina_filtered_s34f_csc_vs_u2af1_illumina_wt_dmso')
