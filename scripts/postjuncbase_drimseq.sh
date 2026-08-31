#!/bin/bash
set -euo pipefail

#This script and analysis was written and performed by Selam Mehreteab (smehrete@ucsc.edu/selammeh2004@gmail.com)

##########################################################################################################################################post-JuncBASE and pre-Drimseq commands
#######################################################################################################################

python3 /smehrete/Juncbase2nd_try/juncbase_filteroutput_0925-2.py \
	-j /smehrete/Juncbase2nd_try/JuncBaseEventReadCounts_020425_AS_exclusion_inclusion_counts_lenNorm.txt \
	-o /smehrete/Juncbase2nd_try/JuncBaseEventReadCounts_092625_AS_exclusion_inclusion_counts_lenNorm \
	--gtf /smehrete/gencode.v33.primary_assembly.annotation.gtf

python /private/groups/brookslab/smehrete/Juncbase2nd_try/filtering.drimseq.input.tables.py JuncBaseEventReadCounts_092625_AS_exclusion_inclusion_counts_lenNorm.drimformat.tsv

python3 /smehrete/Juncbase2nd_try/juncbase_drim_prefilter.py \
	-d /smehrete/Juncbase2nd_try/wtdmso_wtcsc.tsv \
	--cond1 wtdmso \
	--cond2 wtcsc \
	-o /smehrete/Juncbase2nd_try/092925.wtdmsovswtcsc.drimformat.mysamples.filtered.tsv

python3 /smehrete/Juncbase2nd_try/juncbase_drim_prefilter.py \
	-d /smehrete/Juncbase2nd_try/wtdmso_mtdmso.tsv \
	--cond1 wtdmso \
	--cond2 mtdmso \
	-o /smehrete/Juncbase2nd_try/092925.wtdmsovsmtdmso.drimformat.mysamples.filtered.tsv


python3 /smehrete/Juncbase2nd_try/juncbase_drim_prefilter.py \
	-d /smehrete/Juncbase2nd_try/wtdmso_mtcsc.tsv \
	--cond1 wtdmso \
	--cond2 mtcsc \
	-o /smehrete/Juncbase2nd_try/092925.wtdmsovsmtcsc.drimformat.mysamples.filtered.tsv
