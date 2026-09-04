# order to run m6anet steps and all the pre and post processing steps

##pre m6anet inference scripts to run
### 1.run convert_trancriptome_coord_to_genomic.py to convert transcript positions into genomic positions
### 2.map the read indexes in data.indiv file to the read names from the eventalign file and the read names to the sample ids by running optimized_mapping_samplenames.py so when we run m6anet on our aggregated samples we can tell it to group by genomic coordinates and sample ids

##to run m6anet inference
###run calculate_site_proba.py

##post processing to do statistics
### 1.run mod_ratio_stats.py - helps extract the num of mod and unmod reads for stats
### 2.run fisherspair_input_filtering.py - helps to format files for fishers test
### 3.run m6anet_modratio_fishers_test.py - run fishers test between wtdmso vs wtcsc, wtdmso vs mtdmso, wtdmso vs mtcsc
### 4.mapping_m6apos_togene.py - maps which gene the m6a site overlaps with and annotates, also calculates the delta mod ratio



