#!/bin/bash
set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)

##########################################################################################################################################JuncBASE command
#######################################################################################################################################

python /JuncBASE/run_preProcess_by_chr_step1.py \
	-i /smehrete/Juncbase2nd_try/u2af1csc_samp2bam_docker.txt \
	-o /smehrete/Juncbase2nd_try/sample_folders \
	--preProcess_options "--unique -j /smehrete/gencode.v33.primary_assembly.annotation.txt" \
	-p 20 

find /smehrete/Juncbase2nd_try/sample_folders -name "*_chrGL*" -exec rm -r {} \;

find /smehrete/Juncbase2nd_try/sample_folders -name "*_chrKI*" -exec rm -r {} \;

find /smehrete/Juncbase2nd_try/sample_folders -name "*_chrM*" -exec rm -r {} \;

python /JuncBASE/disambiguate_junctions.py \
	-i /smehrete/Juncbase2nd_try/sample_folders \
	-g /smehrete/GRCh38.u2af1_fix.v1.2020_04_01.fa \
	--by_chr --majority_rules 

python /JuncBASE/preProcess_getASEventReadCounts_by_chr_step2.py \
	-i /smehrete/Juncbase2nd_try/sample_folders \
	--by_chr 

python /JuncBASE/run_preProcess_step3_by_chr.py \
	--input_dir /smehrete/Juncbase2nd_try/sample_folders \
	--min_overhang 6 \
	--num_processes 20

python /JuncBASE/createPseudoSample.py \
	-i /smehrete/Juncbase2nd_try/sample_folders \
	-s W1D2_S24 \
       	--by_chr 

python /JuncBASE/run_getASEventReadCounts_multiSample.py \
	-s /smehrete/Juncbase2nd_try/samplemanifest.txt \
	-i /smehrete/Juncbase2nd_try/sample_folders \
	-o /smehrete/Juncbase2nd_try/getASEventReadCount \
	--sqlite_db_dir /smehrete/ \
	--txt_db1 gencode.v33 \
	--txt_db2 gencode.v33.basic \
	--jcn_seq_len 290 \
	-p 20 \
	--by_chr 


python /JuncBASE/run_createAS_CountTables.py \
	-d /smehrete/Juncbase2nd_try/getASEventReadCount/ \
	-i /smehrete/Juncbase2nd_try/sample_folders/ \
	--jcn_seq_len 290 \
	-s W1D2_S24,W1D3_S32,W2D2_S22,W2D3_S30,W1C2_S25,W1C3_S33,W2C2_S23,W2C3_S31,M1D2_S28,M1D3_S36,M2D2_S26,M2D3_S34,M1C2_S29,M1C3_S37,M2C2_S27,M2C3_S35 \
	--num_processes 20 

python /JuncBASE/combine_createAS_CountTables_by_chr.py \
	-d /smehrete/Juncbase2nd_try/getASEventReadCount/ \
	-o /smehrete/Juncbase2nd_try/JuncBaseEventReadCounts_020425





