import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import os
import csv
import logging
from typing import Dict, List, Any

bed_file = '/private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.bed'

#parsing through the bed file
bed_intervals = defaultdict(list)
with open(bed_file, 'r') as bed:
    for line in bed:
        if line.startswith('#') or line.strip() == "":
            continue
        fields = line.strip().split('\t')
        if len(fields) < 12:
            logging.warning(f"Skipping line due to insufficient columns: {line.strip()}")
            continue
        chrom = fields[0]
        chrom_start = int(fields[1])
        chrom_end = int(fields[2])
        feature_id = fields[3]
        split_feature = feature_id.split('_')
        #print(split_feature)
        if len(split_feature) == 3:
            transcript_id = split_feature[1]
            gene_id = split_feature[2]
        elif len(split_feature) ==2:
            transcript_id = split_feature[0]
            gene_id = split_feature[1]
        else:
            transcript_id = np.nan
            gene_id = np.nan
        #print(transcript_id)
       #print(gene_id)
        strand = fields[5]
        bed_intervals[chrom].append((chrom_start,chrom_end,transcript_id, gene_id))
#        print(bed_intervals)

#parse through the input file

m6a_output = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/fishers_results/wtdmso_wtcsc_pairwise_fisher_unfiltered.tsv'

df = pd.read_csv(m6a_output, sep = '\t')

def find_overlap(chrom, pos, bed_intervals):
    transcripts = set()
    genes = set()

    if chrom in bed_intervals:
        for chrom_start, chrom_end, transcript_id, gene_id in bed_intervals[chrom]:
            if chrom_start <= pos < chrom_end:
                if isinstance(transcript_id, str):
                    transcripts.add(transcript_id)
                if isinstance(gene_id, str):
                    genes.add(gene_id)


    transcript_str = ",".join(sorted(transcripts)) if transcripts else np.nan
    genes_str =  ",".join(sorted(genes)) if genes else np.nan
    return transcript_str, genes_str

output = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/07302026_wtdmso_wtcsc_pairwise_fisher_unfiltered.tsv'
gene_name_file = '/private/groups/brookslab/smehrete/gencode.v33.primary_assembly.annotation.gtf'
read_m6a = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/data.indiv_genomic_samplenames_renamed.csv'

read_m6a = pd.read_csv(read_m6a, sep = ',')

#read_m6a['chrom:genomepos'] = read_m6a['chrom'] + ':' + read_m6a['genomepos']

df[['chrom', 'genome_pos']] = (df['site'].str.split(':', expand=True))

df['genome_pos'] = df['genome_pos'].astype(int)
read_m6a['genomepos'] = read_m6a['genomepos'].astype(int)

df[['transcript_id', 'gene_ids']] = df.apply(lambda row: pd.Series(find_overlap(row['chrom'], int(row['genome_pos']), bed_intervals)), axis = 1)

m6a_set = set(
    zip(
        read_m6a['chrom'],
        read_m6a['genomepos'],
        read_m6a['transcript_id'].apply(lambda name: name.split("_")[1] if len(name.split("_")) >= 2 and name.split("_")[1].upper().startswith("ENST") else name.split("_")[0])
    )
)


#df['transcript_id'] = df.apply(
 #           lambda row: [tx for tx in str(row['transcript_id']).split(",") 
  #                               if (row['chrom'], row['genome_pos'], tx) in m6a_set]
   #             if pd.notna(row['transcript_id']) else [],
    #                axis=1
     #)

df['transcript_id'] = df.apply(
            lambda row: ",".join(
                        [tx for tx in str(row['transcript_id']).split(",")
                                     if (row['chrom'], row['genome_pos'], tx) in m6a_set]
                            ) if pd.notna(row['transcript_id']) else np.nan,
                axis=1
                )


#df = df.explode('transcript_id').reset_index(drop=True)

#read_m6a['transcript_id_only'] = df['transcript_id'].apply(lambda x: [part for part in x.replace('-', '_').split('_') if part.startswith("ENST")][0] if any(part.startswith("ENST") for part in x.replace('-', '_').split('_')) else x)

#df['tx_match'] = df.apply(lambda row: find_tx(row['chrom'], row['genomepos'], tx_intervals), axis=1)

#df_filtered = df[df.apply(lambda row: row['transcript_id_only'] in row['tx_match'].split(','), axis=1)]

#df = df.drop(columns=['transcript_id'])
#Now i want to map all the ensg ids to gene_ids
#I want to print all the rows that have more than 2 values for gene_id
def cleanup_gene(gene_ids):
    if pd.isna(gene_ids):
        return []
    return re.findall(r'ENSG\d+\.\d+', gene_ids)

df['gene_ids'] = df['gene_ids'].apply(cleanup_gene)

df = df.explode('gene_ids').reset_index(drop=True)


with open('/private/groups/brookslab/smehrete/RNA_modification/10272025.m6anet.results/gencode.v33.primary_assembly.annotation.gtf') as f:
    gtf = list(f)

gtf = [x for x in gtf if not x.startswith('#')]
gtf = [x for x in gtf if 'gene_name "' in x]

gtf_list = list(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf))

gtfset = set(gtf_list)
gtfdict = dict(gtfset)


df['gene_name'] = df['gene_ids'].map(gtfdict)

#find delta mod ratio
df['total_reads_cond1'] = df['mod_cond1'] + df['unmod_cond1']
df['total_reads_cond2'] = df['mod_cond2'] + df['unmod_cond2']

df['mod_ratio_cond1'] = df['mod_cond1'] / df['total_reads_cond1'].replace(0, np.nan)
df['mod_ratio_cond2'] = df['mod_cond2'] / df['total_reads_cond2'].replace(0, np.nan)

df['delta_modratio'] = df['mod_ratio_cond2'] - df['mod_ratio_cond1']

#uncomment if doing analysis with signficant genes this is only for headmap
#df = df[df['delta_modratio'].abs() >= 0.1].copy()
#df = df[df['pvalue'] <= 0.05].copy()

df = df.drop(columns=['total_reads_cond1', 'total_reads_cond2'])
#print(df)

#first match the genomepos and also conditions and extract the probability modified values for each condition then calculate the delta mod prob
wtcsc_mod_prob = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/wtcsc_delta_site_prob_vals.csv'
mtdmso_mod_prob = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/mtdmso_delta_site_prob_vals.csv'
mtcsc_mod_prob = '/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/mtcsc_delta_site_prob_vals.csv'


wtcsc_mod_prob = pd.read_csv(wtcsc_mod_prob, sep = ',')
mtdmso_mod_prob = pd.read_csv(mtdmso_mod_prob, sep=',')
mtcsc_mod_prob = pd.read_csv(mtcsc_mod_prob, sep = ',')


wtcsc_mod_prob = wtcsc_mod_prob.rename(columns={'chrom:genomepos':'site'})
mtcsc_mod_prob = mtdmso_mod_prob.rename(columns={'chrom:genomepos':'site'})
mtdmso_mod_prob = mtcsc_mod_prob.rename(columns={'chrom:genomepos':'site'})
#print(mtdmso_mod_prob )

#wtcsc
#df = df.merge(wtcsc_mod_prob,on="site",how="left")
#df = df.dropna(subset=['delta_prob'])
#print(df)
#mtdmso
#df = df.merge(mtdmso_mod_prob,on="site",how="inner")
#df = df.dropna(subset=['delta_prob'])
#print(df)
#mtcsc
#df = df.merge(mtcsc_mod_prob,on="site",how="inner")
#df = df.dropna(subset=['delta_prob'])
print(df)


df.to_csv(output, sep='\t', index=False)

#df["base_site"] = df.index.str
