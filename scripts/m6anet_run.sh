#!/bin/bash
set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)


m6anet dataprep --eventalign /scratch/smehrete/07.17.2026.all_samples_cdna_drna.eventalign.noreadname.txt \
	--out_dir /scratch/smehrete/07212026_m6anet_results \
	--n_processes 4
# 1.run convert_trancriptome_coord_to_genomic.py to convert transcript positions into genomic positions
# 2.map the read indexes in data.indiv file to the read names from the eventalign file and the read names to the sample ids by running optimized_mapping_samplenames.py so when we run m6anet on our aggregated samples we can tell it to group by genomic coordinates and sample ids
# 3.after mapping the sample names then run the combining reps script to format the sample_ids to combine reps and standardize the naming of the sample names


#After preprocessing using above scripts then run this code to run m6anet inference
python calculate_site_proba.py \
	-i /scratch/smehrete/07212026_m6anet_result/data.indiv_genomic_samplenames.csv \
	-o /scratch/smehrete/07212026_m6anet_result/data.site_withsampleids.proba.csv


